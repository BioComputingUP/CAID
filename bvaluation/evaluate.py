# module imports
import os
import numpy as np
import logging
import warnings
import pandas as pd
import argparse
# function imports
from itertools import groupby
# relative imports
from bvaluation.assessment import set_logger, load_names
from bvaluation.assessment.dataset import ReferencePool, PredictionEntry
from bvaluation.assessment.evaluation import Evaluation
from bvaluation.assessment.dataset import Reference, Prediction

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


def parse_args(arglist=None):
    parser = argparse.ArgumentParser(
        prog='evaluate.py', description="Binary evaluation package",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('reference',
                        help='reference file to which predictions are to be compared')

    parser.add_argument('prediction', nargs='+',
                        help='file(s) containing binary predictions for reference targets.')

    parser.add_argument('-s', '--strictNegatives', action='store_true', default=False,
                        help="filter out (of ref and pred) all unlabeled positions "
                             "in reference ('-')")
    parser.add_argument('-o', '--outputDir', help='directory where the output will be written',
                        default='.')
    parser.add_argument('-b', '--labels', help='text file with prediction file name (without path) and label '
                                               'desired in output organized in two columns')
    parser.add_argument('--prefix', help='prefix of the output files basename', default=None)
    parser.add_argument('--suffix', help='suffix of the output files basename', default=None)

    # log options
    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args() if arglist is None else parser.parse_args(arglist)
    return args


def strip_split(string: str, expected=None):
    if string:
        s = string.strip('\n').rsplit('\t')
        if expected and len(s) != expected:
            logging.error('expecting %i elements from split, got %i', expected, len(s))
            raise ValueError('expecting %i elements from split, got %i', expected, len(s))
        return s


def parse_prediction_file(predfile: str, reference_ids):
    reference_ids = set(reference_ids)

    with open(predfile) as fhandle:
        faiter = (x[1] for x in groupby(fhandle, lambda line: line[0] == ">"))
        for header in faiter:
            header = next(header).strip()[1:]
            body = next(faiter)
            if header in reference_ids:
                try:
                    positions, residues, scores, states = zip(*map(strip_split, body))
                except ValueError:
                    logging.error('error while parsing prediction %s', header)
                    raise ValueError('error while parsing prediction %s', header)

                pred_entry = PredictionEntry(positions, residues, scores, states)
                yield header, pred_entry


def iterate_prediction_files(pred_files: list):
    for pred_file in pred_files:
        pred_name = os.path.basename(os.path.splitext(pred_file)[0])
        yield pred_file, pred_name


def build_output_basename(reference: str, strict_negs: bool, outdir: str, prefix: str = None, suffix: str = None):
    basename = '{}_'.format(prefix) if prefix else ''
    basename += '{}-negs_'.format('str' if strict_negs is True else 'len')
    basename += os.path.splitext(os.path.basename(reference))[0].replace('_', '-')
    basename = os.path.join(os.path.abspath(outdir), basename)
    basename += '{}_'.format(suffix) if suffix else ''
    logging.info('output basename %s', basename)
    return basename


def evaluate(reference_dataset: ReferencePool, pred_file: str, lbl: str):
    evl, roc_evl = None, None
    prediction = Prediction()
    reference = Reference()

    logging.info('evaluating {}'.format(lbl))

    for accession, pred in parse_prediction_file(pred_file, reference_dataset.keys()):
        prediction.add_accession(accession)
        reference.add_accession(accession)

        ref_entry_states = reference_dataset[accession]['states']
        nan_indexes = np.isnan(ref_entry_states)
        ref_entry_states = ref_entry_states[~nan_indexes]

        reference.states.append(ref_entry_states)

        if 'scores' in pred:
            prediction.scores.append(pred['scores'][~nan_indexes])
        if 'states' in pred:
            prediction.states.append(pred['states'][~nan_indexes])

    prediction.set_coverage(reference.accessions_set)
    prediction.set_merged_states()
    reference.set_merged_states()

    if prediction.coverage != 0:
        evl = Evaluation()
        evl.calc_overall_scores(reference.mstates, prediction.mstates)
        evl.calc_curves(reference.mstates, prediction.mscores)
        evl.calc_average_instance_scores(reference.states, prediction.states,
                                         save_perinstance_scores=True,
                                         accessions=reference.accessions)

        youdensj_roc = evl.curves.roc.youdensj
        if youdensj_roc is not None and not np.isnan(youdensj_roc):
            roc_evl = Evaluation()
            prediction.apply_cutoff(youdensj_roc)
            roc_evl.calc_overall_scores(reference.mstates, prediction.mstates)
            roc_evl.calc_average_instance_scores(reference.states, prediction.states,
                                                 save_perinstance_scores=True,
                                                 accessions=reference.accessions)

    else:
        logging.warning('predictor %s had no predictions for entries in reference', lbl)

    return evl, roc_evl, prediction.coverage


def main(cm_args):
    # make reference path absolute
    cm_args.reference = os.path.abspath(cm_args.reference)
    # check if reference file exist: raise an error if it doesn't
    if not os.path.isfile(cm_args.reference):
        logging.critical('File not found %s', cm_args.reference)
        raise FileNotFoundError(cm_args.reference)
    # setup logger configurations
    set_logger(cm_args.log, cm_args.logLevel)
    # build output basename based on command line arguments
    output_basename = build_output_basename(cm_args.reference, cm_args.strictNegatives, cm_args.outputDir,
                                            cm_args.prefix, cm_args.suffix)
    # load labels to be used as classifier names in output
    labels = load_names(cm_args.labels) if cm_args.labels else dict()
    # log initialization data
    logging.info('reference: %s', cm_args.reference)
    logging.info('predicition(s): %s', ' '.join(cm_args.prediction))
    logging.info('negative definition: %s', 'strict' if cm_args.strictNegatives is True else 'lenient')

    ref_pool = ReferencePool(cm_args.reference, strict_negatives=cm_args.strictNegatives)

    metrics_table = pd.DataFrame()
    metrics_table_r = pd.DataFrame()
    perinstance_table = pd.DataFrame()
    roc_points = list()
    prc_points = list()

    for prediction_file, predname in iterate_prediction_files(cm_args.prediction):
        label = labels.get(os.path.basename(prediction_file), predname)
        evaluation, roc_evaluation, coverage = evaluate(ref_pool, prediction_file, label)

        if evaluation is not None:
            restabl = evaluation.get_scores_asdict()
            restabl.update({'Cov': coverage})
            metrics_table = pd.concat([metrics_table, pd.DataFrame(restabl, index=[label])],
                                      axis=0, sort=False)

            pis = pd.DataFrame(evaluation.get_perinstance_scores(insert=label))
            colnames = pis.iloc[0][2:]
            pis = pis.iloc[1:].rename(columns=colnames).set_index([0, 1]).astype(np.float)
            perinstance_table = perinstance_table.append(pis, sort=False)

        if roc_evaluation is not None:
            roc, prc = evaluation.get_curves_repr().split('\n')
            roc_points.append('{} {}'.format(label, roc))
            prc_points.append('{} {}'.format(label, prc))

            restablr = roc_evaluation.get_scores_asdict()
            restablr.update({'Cov': coverage})
            metrics_table_r = pd.concat([metrics_table_r, pd.DataFrame(restablr, index=[label])],
                                        axis=0, sort=False)

    metrics_table.to_csv('{}_scores.csv'.format(output_basename), float_format='%.3f')
    metrics_table_r.to_csv('{}_redefScores.csv'.format(output_basename), float_format='%.3f')
    perinstance_table.to_csv('{}_perInstanceScores.csv'.format(output_basename),
                             float_format='%.3f')

    if roc_points:
        with open('{}_rocPoints.txt'.format(output_basename), 'w') as f:
            f.write('\n'.join(roc_points))

    if prc_points:
        with open('{}_prcPoints.txt'.format(output_basename), 'w') as f:
            f.write('\n'.join(prc_points))


if __name__ == '__main__':
    # parse command line arguments
    main(parse_args())
