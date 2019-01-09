import os
import copy
import logging
import numpy as np
import pandas as pd
from itertools import chain
import multiprocessing as mp
# relative imports
from caid.dataset import Prediction
from caid import parse_config, parse_args, set_logger
from caid.evaluation import Evaluation
from caid.dataset import ReferencePool

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


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
        scores = pd.concat([scores, pd.DataFrame(evl.get_scores(dstr=dict), index=[i])],
                           axis=0, sort=False)

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

    ref_pool = ReferencePool(os.path.abspath(args.reference),
                             strict_negatives=args.strictNegatives).make_pure_reference()

    out_basename = os.path.splitext(os.path.basename(args.reference))[0].replace('_', '-')
    outdir = dict(conf.items('prj_directories'))['baseline']
    out_basename = os.path.join(outdir, out_basename)
    shuffle_reference(ref_pool, out_basename)
