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
from typing import Iterable, Generator, Tuple, List

# TODO: is this a solution? seems so
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

# relative imports
from bvaluation.assessment import set_logger, load_names
from bvaluation.assessment.dataset import PredictionEntry
from bvaluation.assessment.evaluation import Evaluation, ConsensusEvaluation
from bvaluation.assessment.dataset import Reference, Prediction, ReferencePool

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


def parse_reference(ref_file: str, pttrn: dict=None) -> (dict, set):
    """Load reference as dict from reference Fasta, return dict and set of accessions

    Reference dict is built to be later converted in pd.DataFrame

    :param ref_file: reference file
    :param pttrn: pattern dictionary. Each state of the reference will be interpreted looking for a
    corresponding value in `pptrn` dict
    :return:
    """

    pattern = {'0': 0.0, '1': 1.0, '-': np.nan} if pttrn is None else pttrn

    ref = {('ref', 'states'): {}, ('ref', 'seq'): {}}
    accs = set()

    with open(ref_file) as f:
        faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
        for header in faiter:
            header = next(header)
            if header[0] != '#':
                acc, *desc = header[1:].strip().split()
                seq, states = map(str.strip, next(faiter))
                states = np.array([pattern.get(s, np.nan) for s in states], dtype=np.float)

                accs.add(acc)
                for i, (st, aa), in enumerate(zip(states, seq)):
                    ref[('ref', 'states')][(acc, i)] = float(st)
                    ref[('ref', 'seq')][(acc, i)] = aa

    return ref, accs


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


def parse_prediction(predfile: str, reference_ids: Iterable) -> Generator:
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


def build_output_basename(reference: str, outdir: str,
                          name_structure: List[str],
                          prefix: str = None, suffix: str = None) -> str:
    """
    Compose output file basename based on input and parameters

    :param reference: reference file
    :param outdir: path to output dir
    :param name_structure: structure of the output file basename
    :param prefix: prefix of the output basename
    :param suffix: suffix of the output basename
    :return: outout basename
    """
    basename = list()

    name_composer = {'b': os.path.splitext(os.path.basename(reference))[0]}

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


def evaluate(reference_dataset, pred_file: str, lbl: str, opt_thr='roc'):
    """
    Evaluate predictions coming from a single file against the reference

    :param reference_dataset: reference pool
    :param pred_file: path to prediction file
    :param lbl: evaluation label
    :return: evaluation, ROC-based evaluation, prediction coverage
    """
    evl, roc_evl, fmax_evl = None, None, None
    roc_thr, fmax_thr = None, None
    predobj, roc_predobj, fmax_predobj = None, None, None
    prediction = Prediction()
    reference = Reference()
    reference.accessions = list(reference_dataset.keys())
    reference.accessions_set = set(reference.accessions)

    logging.info('evaluating {}'.format(lbl))

    for accession, pred in parse_prediction(pred_file, reference_dataset.keys()):

        prediction.add_accession(accession)
        # reference.add_accession(accession)

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

    if prediction.coverage[0] > 0:
        evl = Evaluation()
        evl.calc_overall_scores(reference.mstates, prediction.mstates)
        evl.calc_curves(reference.mstates, prediction.mscores)
        evl.calc_average_instance_scores(reference.states, prediction.states,
                                         save_perinstance_scores=True,
                                         accessions=prediction.accessions)

        predobj = prediction.copy()

        roc_thr = evl.curves.roc.youdensj
        fmax_thr = evl.curves.prc.fmax

        if roc_thr is not None and not np.isnan(roc_thr):
            roc_evl = Evaluation()
            prediction.apply_cutoff(roc_thr)
            roc_evl.calc_overall_scores(reference.mstates, prediction.mstates)
            roc_evl.calc_average_instance_scores(reference.states, prediction.states,
                                                 save_perinstance_scores=True,
                                                 accessions=prediction.accessions)

            roc_predobj = prediction.copy()

        if fmax_thr is not None and not np.isnan(fmax_thr):
            fmax_evl = Evaluation()
            prediction.apply_cutoff(fmax_thr)
            fmax_evl.calc_overall_scores(reference.mstates, prediction.mstates)
            fmax_evl.calc_average_instance_scores(reference.states, prediction.states,
                                                 save_perinstance_scores=True,
                                                 accessions=prediction.accessions)

            fmax_predobj = prediction.copy()

    else:
        logging.warning('predictor %s had no predictions for entries in reference', lbl)

    return (evl, None, predobj), (roc_evl, roc_thr, roc_predobj), (fmax_evl, fmax_thr, fmax_predobj)


def bvaluation(reference: str, prediction: list, outdir: str, prefix: str = None,
               suffix: str = None, replace_undefined: int = None,
               name_structure: Iterable[str] = ('a', 'b'),
               labels: str = None, log: str = None, log_level="ERROR"):
    """
    Evaluate one or more prediction files against a given reference

    :param ref_pool: reference file
    :param prediction: one or more prediction files
    :param outdir: path to output dir
    :param prefix: prefix of output basename
    :param suffix: suffix of output basename
    :param name_structure: structure (presence and order) of elements in the output basename
    :param labels: file containing labels associated to prediction files
    :param log: log file
    :param log_level: log level
    """

    # make reference path absolute
    ref = os.path.abspath(reference)

    # check if reference file exist: raise an error if it doesn't
    if not os.path.isfile(reference):
        logging.critical('File not found %s', reference)
        raise FileNotFoundError(reference)

    # setup logger configurations
    set_logger(log, log_level)

    # build output basename based on command line arguments
    output_basename = build_output_basename(reference, outdir, name_structure, prefix, suffix)

    # load labels to be used as classifier names in output
    labels = load_names(labels) if labels else dict()
    # log initialization data
    logging.info('reference: %s', reference)
    logging.info('predicition(s): %s', ' '.join(prediction))
    logging.info('undefined replaced with %s', replace_undefined)

    ref_pool = ReferencePool(reference, replace_undefined)

    metrics_table = pd.DataFrame()
    metrics_table_roc = pd.DataFrame()
    metrics_table_fmax = pd.DataFrame()
    perinstance_table = pd.DataFrame(columns=ref_pool.keys())
    perinstance_table_roc = pd.DataFrame(columns=ref_pool.keys())
    perinstance_table_fmax = pd.DataFrame(columns=ref_pool.keys())
    pred_stack = pd.DataFrame(columns=ref_pool.keys())
    roc_points = list()
    prc_points = list()

    # copy reference dict to seed the table of reference + predictors
    for prediction_file, predname in iterate_prediction_files(prediction):
        label = labels.get(prediction_file, predname)
        normal_evl, roc_evl, fmax_evl  = evaluate(ref_pool, prediction_file, label, opt_thr='prc')

        evaluation, _, predobj = normal_evl
        roc_evaluation, roc_thr, roc_predobj = roc_evl
        fmax_evaluation, fmax_thr, fmax_predobj = fmax_evl
    
        if evaluation is not None:
            restabl = evaluation.get_scores_asdict()
            restabl.update(dict(zip(('npred', 'nref'), predobj.coverage)))
            metrics_table = pd.concat([metrics_table, pd.DataFrame(restabl, index=[label])],
                                       axis=0, sort=False)
    
            pis = pd.DataFrame(evaluation.get_perinstance_scores(insert=label))
            colnames = pis.iloc[0][2:]
            pis = pis.iloc[1:].rename(columns=colnames)
            perinstance_table = perinstance_table.append(pis, sort=False)

    
        if roc_evaluation is not None:
            roc, prc = evaluation.get_curves_repr().split('\n')
            roc_points.append('{} {}'.format(label, roc))
            prc_points.append('{} {}'.format(label, prc))
    
            roc_evaluation.curves = evaluation.curves
            restabl_roc = roc_evaluation.get_scores_asdict()
            restabl_roc.update(dict(zip(('npred', 'nref'), predobj.coverage)))
            metrics_table_roc = pd.concat([metrics_table_roc, pd.DataFrame(restabl_roc, index=[label])],
                                         axis=0, sort=False)
    
            pis_roc = pd.DataFrame(roc_evaluation.get_perinstance_scores(insert=label))
            colnames_r = pis_roc.iloc[0][2:]
            pis_roc = pis_roc.iloc[1:].rename(columns=colnames_r)
            perinstance_table_roc = perinstance_table_roc.append(pis_roc, sort=False)

        if fmax_evaluation is not None:
            roc, prc = evaluation.get_curves_repr().split('\n')
            roc_points.append('{} {}'.format(label, roc))
            prc_points.append('{} {}'.format(label, prc))

            fmax_evaluation.curves = evaluation.curves
            restabl_fmax = fmax_evaluation.get_scores_asdict()
            restabl_fmax.update(dict(zip(('npred', 'nref'), predobj.coverage)))
            metrics_table_fmax = pd.concat(
                [metrics_table_fmax, pd.DataFrame(restabl_fmax, index=[label])],
                axis=0, sort=False)

            pis_fmax = pd.DataFrame(fmax_evaluation.get_perinstance_scores(insert=label))
            colnames_r = pis_fmax.iloc[0][2:]
            pis_fmax = pis_fmax.iloc[1:].rename(columns=colnames_r)
            perinstance_table_fmax = perinstance_table_fmax.append(pis_fmax, sort=False)
    
        df = pd.DataFrame(predobj.zip()).transpose()
        
        df = df.rename(columns=df.iloc[0]).iloc[1:].rename(index={1:predname})
        pred_stack = pred_stack.append(df, sort=False)

    # calc consensus evaluation
    cons_eval = ConsensusEvaluation()
    pstack, ptargt, merged_reference = cons_eval.zip_predstack_ref(ref_pool, pred_stack)
    cm = cons_eval.get_consensus_confmat(pstack, merged_reference)

    cons_eval.save_pred_stack(pstack, ptargt, '{}_predstack.csv'.format(output_basename))
    pd.DataFrame(cm).to_csv('{}_consensusCM.csv'.format(output_basename))

    perinstance_table = perinstance_table.set_index([0, 1]).astype(float)
    if not perinstance_table_roc.empty:
        perinstance_table_roc = perinstance_table_roc.set_index([0, 1]).astype(float)
    if not perinstance_table_fmax.empty:
        perinstance_table_fmax = perinstance_table_fmax.set_index([0, 1]).astype(float)
        
    metrics_table.to_csv('{}_scores.csv'.format(output_basename), float_format='%.3f')
    metrics_table_roc.to_csv('{}_redefScoresROC.csv'.format(output_basename), float_format='%.3f')
    metrics_table_fmax.to_csv('{}_redefScoresFmax.csv'.format(output_basename), float_format='%.3f')

    perinstance_table.to_csv('{}_perInstanceScores.csv'.format(output_basename),
                             float_format='%.3f')
        
    if roc_points:
        with open('{}_rocPoints.txt'.format(output_basename), 'w') as f:
            f.write('\n'.join(roc_points))
            perinstance_table_roc.to_csv('{}_perInstanceRedefScoresROC.csv'.format(output_basename),
                                         float_format='%.3f')
    if prc_points:
        with open('{}_prcPoints.txt'.format(output_basename), 'w') as f:
            f.write('\n'.join(prc_points))
            perinstance_table_fmax.to_csv('{}_perInstanceRedefScoresFmax.csv'.format(output_basename),
                                          float_format='%.3f')


if __name__ == '__main__':
    # parse command line arguments
    bvaluation(**vars(parse_args()))
