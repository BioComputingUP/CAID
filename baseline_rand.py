import os
import copy
import logging
import warnings
import argparse
import numpy as np
import pandas as pd
from itertools import chain
import multiprocessing as mp
# relative imports
from bvaluation.assessment.dataset import Prediction
from bvaluation.assessment import parse_config, set_logger
from bvaluation.assessment.evaluation import Evaluation
from bvaluation.assessment.dataset import ReferencePool
from bvaluation.evaluate import build_output_basename

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def parse_args(wd):
    parser = argparse.ArgumentParser(
        prog='caid-assess', description="CAID: Critical Assessment of Intrinsic Disorder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('reference',
                        help='reference file to which predictions are to be compared')

    parser.add_argument('-r', '--replaceUndefined', choices=['0', '1'], default=None,
                        help='replace value for undefined positions (-) in reference. '
                             'By default not applied')
    parser.add_argument('-c', '--conf', type=str,
                        default=os.path.join(wd, 'config.ini'),
                        help="path to an alternative configuration file.")
    parser.add_argument('-o', '--outdir', default='.')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args


def save_prediction(ref_file, acc, seq, scores, states):
    lseq = len(seq)

    scores = scores if scores is not None else [''] * lseq
    states = states if states is not None else [''] * lseq

    assert lseq == len(scores) == len(states), 'sequence, states and scores must have the same len'

    ref_file.write('>{}\n'.format(acc))
    for i, (aa, sc, st) in enumerate(zip(seq, scores, states), 1):
        ref_file.write('{}\t{}\t{:.3f}\t{}\n'.format(i, aa, sc, st))


def shuffle_reference(reference, output_basename, iterations=100):
    outname = '{}_chunk'.format(output_basename)
    logging.info('chunk-wise reference shuffling %i times: %s', iterations, outname)
    save_eval(randomize(reference, chunkwise_baseline, iterations=iterations), outname)
    outname = '{}_chain'.format(output_basename)
    logging.info('chain-wise reference shuffling %i times: %s', iterations, outname)
    save_eval(randomize(reference, chainwise_baseline, iterations=iterations), outname)


def shuffle_chunk(ref):
    predicted_states = list()
    for states in ref:
        tmp = copy.copy(states)
        np.random.shuffle(tmp)
        predicted_states.append(tmp)
    return predicted_states


def shuffle_chain(ref):
    chained_states = list()
    split_indexes = list()
    lenghts = list()

    for states in ref:
        chained_states.extend(states)
        lenghts.append(len(states))
        split_indexes.append(split_indexes[-1] + len(states) if split_indexes else len(states))

    np.random.shuffle(chained_states)
    splitted_states = np.split(np.array(chained_states), split_indexes[:-1])

    return splitted_states


def chunkwise_baseline(reference):
    prediction = Prediction()
    prediction.states = shuffle_chunk(reference.states)
    ps = np.fromiter(chain(*prediction.states), dtype=np.float)
    evaluation = Evaluation()
    evaluation.calc_overall_scores(reference.mstates, ps)
    evaluation.calc_average_instance_scores(reference.states, prediction.states,
                                            accessions=reference.accessions)
    return evaluation


def chainwise_baseline(reference):
    prediction = Prediction()
    prediction.states = shuffle_chain(reference.states)
    ps = np.fromiter(chain(*prediction.states), dtype=np.float)
    evaluation = Evaluation()
    evaluation.calc_overall_scores(reference.mstates, ps)
    evaluation.calc_average_instance_scores(reference.states, prediction.states,
                                            accessions=reference.accessions)
    return evaluation


def randomize(reference, baseline_func, iterations=100):
    with mp.Pool(4) as pool:
        return pool.map(baseline_func, (reference for _ in range(iterations)))


def save_eval(randomized_refs, outbasename):
    scores = pd.DataFrame()
    perinstance_table = pd.DataFrame()

    for i, evl in enumerate(randomized_refs):
        newscores = evl.get_scores_asdict()
        newscores.update(dict(zip(('npred', 'nref'), (len(randomized_refs), len(randomized_refs)))))

        scores = pd.concat([scores, pd.DataFrame(newscores, index=[i])], axis=0, sort=False)

        pis = pd.DataFrame(evl.get_perinstance_scores(insert=str(i)))
        colnames = pis.iloc[0][2:]
        pis = pis.iloc[1:].rename(columns=colnames).set_index([0, 1]).astype(np.float)
        perinstance_table = perinstance_table.append(pis, sort=False)

    perinstance_table.to_csv('{}_perInstanceScores.csv'.format(outbasename))
    scores.to_csv('{}_scores.csv'.format(outbasename))


if __name__ == '__main__':
    args = parse_args(SCRIPT_DIR)
    conf = parse_config(args.conf)
    set_logger(args.log, args.logLevel)

    ref = ReferencePool(
        os.path.abspath(args.reference),
        undefined_replace_value=args.replaceUndefined).make_pure_reference()

    # outdir = dict(conf.items('data_directories'))['baseline']
    # out_basename = 'len-negs_' if args.strictNegatives is False else 'str-negs_'
    # out_basename += os.path.splitext(os.path.basename(args.reference))[0].replace('_', '-')
    # out_basename = os.path.join(outdir, out_basename)

    out_basename = build_output_basename(args.reference, args.replaceUndefined,
                                         args.outdir, ['a', 'b'])

    shuffle_reference(ref, out_basename)