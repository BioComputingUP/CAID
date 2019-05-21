# module imports
import os
import sys
import numpy as np
import logging
import warnings
import pandas as pd
import argparse
# function imports
from itertools import groupby
from typing import Iterable, Generator, Tuple

# TODO: is this a solution? seems so
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

# relative imports
from bvaluation.assessment import set_logger, load_names
from bvaluation.assessment.dataset import ReferencePool, PredictionEntry
from bvaluation.assessment.evaluation import Evaluation
from bvaluation.assessment.dataset import Reference, Prediction


warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


def parse_args(arglist: list = None):
    """
    Parse command line arguments or arguments passed as a list

    :param arglist: list of string representing command line arguments
    :return: arguments namespace
    """
    parser = argparse.ArgumentParser(
        prog='evaluate.py', description="Binary evaluation package",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # positional arguments
    parser.add_argument('reference',
                        help='reference file to which predictions are to be compared')

    parser.add_argument('prediction', nargs='+',
                        help='file(s) containing binary predictions for reference targets.')

    # computation option
    parser.add_argument('-r', '--replace_undefined', choices=['0', '1'], default=None,
                        help='replace value for undefined positions (-) in reference. '
                             'By default not applied')
    # output options
    parser.add_argument('-o', '--outdir', default='.',
                        help='directory where the output will be written')
    parser.add_argument('-b', '--labels',
                        help='text file with prediction file name (without path) and label desired '
                             'in output organized in two columns')

    # output filename options
    parser.add_argument('--prefix', help='prefix of the output files basename', default=None)
    parser.add_argument('--suffix', help='suffix of the output files basename', default=None)
    parser.add_argument('-n', '--name_structure', nargs='*', default=['a', 'b'], choices=['a', 'b'],
                        help='strucutre of the output file name. a: undefined replace value; '
                             'b: reference basename. By default \'a, b\'.')

    # log options
    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--log_level", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args() if arglist is None else parser.parse_args(arglist)
    return args


def strip_split(string: str, expected: int = None) -> list:
    """
    remove newlines and split by tabs, raise an error for unexpected number of elements from split

    :param string: string to modify
    :param expected: number of elements expected splitting the string by tabs
    :return: splitted string
    """
    if string:
        splitted = string.strip('\n').rsplit('\t')
        if expected and len(splitted) != expected:
            logging.error('expecting %i elements from split, got %i', expected, len(splitted))
            raise ValueError('expecting %i elements from split, got %i', expected, len(splitted))
        return splitted


def parse_prediction_file(predfile: str, reference_ids: Iterable) -> Generator:
    """
    Lazy load predictions for one target at the time from file

    :param predfile: file containing predictions for targets
    :param reference_ids: filter targets by their presence in reference
    :return: predictions per target
    :rtype: Generator
    """
    reference_ids = set(reference_ids)

    with open(predfile) as fhandle:
        faiter = (x[1] for x in groupby(fhandle, lambda line: line[0] == ">"))
        for header in faiter:
            header = next(header).strip()[1:]
            body = next(faiter)
            if header in reference_ids:
                try:
                    positions, _, scores, states = zip(*map(strip_split, body))
                except ValueError:
                    logging.error('error while parsing prediction %s', header)
                    raise ValueError('error while parsing prediction %s', header)

                pred_entry = PredictionEntry(positions, scores, states)
                yield header, pred_entry


def iterate_prediction_files(pred_files: list) -> Generator:
    """
    Lazy load prediction files

    :param pred_files: list of prediction files
    :type pred_files:
    :return: prediction files
    :rtype: Generator
    """
    for pred_file in pred_files:
        pred_name = os.path.basename(os.path.splitext(pred_file)[0])
        yield pred_file, pred_name


def build_output_basename(reference: str, undef_repl: int, outdir: str,
                          name_structure: Iterable[str],
                          prefix: str = None, suffix: str = None) -> str:
    """
    Compose output file basename based on input and parameters

    :param reference: reference file
    :param undef_repl: value (if any) with which to replace undefined in reference
    :param outdir: path to output dir
    :param name_structure: structure of the output file basename
    :param prefix: prefix of the output basename
    :param suffix: suffix of the output basename
    :return: outout basename
    """
    basename = list()

    name_composer = {'a': 'urv-{}'.format(undef_repl),
                     'b': os.path.splitext(os.path.basename(reference))[0]}

    if undef_repl is None:
        del name_composer['a']

    if prefix is not None:
        basename.append(prefix)

    for label in name_structure:
        if label in name_composer:
            basename.append(name_composer[label])
        else:
            logging.warning('A relevant option is absent from the name composition')

    if suffix is not None:
        basename.append(suffix)

    basename = os.path.join(os.path.abspath(outdir), '_'.join(basename))
    logging.info('output basename %s', basename)
    return basename


def evaluate(reference_dataset: ReferencePool,
             pred_file: str, lbl: str) -> Tuple[Evaluation, Evaluation, float]:
    """
    Evaluate predictions coming from a single file against the reference

    :param reference_dataset: reference pool
    :param pred_file: path to prediction file
    :param lbl: evaluation label
    :return: evaluation, ROC-based evaluation, prediction coverage
    """
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


def bvaluation(reference: str, prediction: list, outdir: str, prefix: str = None,
               suffix: str = None, replace_undefined: int = None,
               name_structure: Iterable[str] = ('a', 'b'), labels: str = None, log: str = None, log_level="ERROR"):
    """
    Evaluate one or more prediction files against a given reference

    :param reference: reference file
    :param prediction: one or more prediction files
    :param outdir: path to output dir
    :param prefix: prefix of output basename
    :param suffix: suffix of output basename
    :param replace_undefined: value (if any) with which to replace undefined elements in reference
    :param name_structure: structure (presence and order) of elements in the output basename
    :param labels: file containing labels associated to prediction files
    :param log: log file
    :param log_level: log level
    """
    # make reference path absolute
    reference = os.path.abspath(reference)
    # check if reference file exist: raise an error if it doesn't
    if not os.path.isfile(reference):
        logging.critical('File not found %s', reference)
        raise FileNotFoundError(reference)
    # setup logger configurations
    set_logger(log, log_level)
    # build output basename based on command line arguments
    output_basename = build_output_basename(reference, replace_undefined,
                                            outdir, name_structure,
                                            prefix, suffix)
    # load labels to be used as classifier names in output
    labels = load_names(labels) if labels else dict()
    # log initialization data
    logging.info('reference: %s', reference)
    logging.info('predicition(s): %s', ' '.join(prediction))
    logging.info('undefined replaced with %s', replace_undefined)

    ref_pool = ReferencePool(reference, replace_undefined)

    metrics_table = pd.DataFrame()
    metrics_table_r = pd.DataFrame()
    perinstance_table = pd.DataFrame()
    roc_points = list()
    prc_points = list()

    for prediction_file, predname in iterate_prediction_files(prediction):
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
    bvaluation(**vars(parse_args()))
