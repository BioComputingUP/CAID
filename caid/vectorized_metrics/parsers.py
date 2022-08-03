import os
from pathlib import Path
from itertools import groupby
import logging
import numpy as np


def parse_reference(ref_file: str, pttrn: dict=None) -> (dict, set):
    """Load reference as dict from reference Fasta, return dict and set of accessions

    Reference dict is built to be later converted in pd.DataFrame

    :param ref_file: reference file
    :param pttrn: pattern dictionary. Each state of the reference will be interpreted looking for a
    corresponding value in `pattern` dict
    :return:
    """
    ref_file = Path(ref_file)
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
                states = np.array([pattern.get(s, np.nan) for s in states], dtype=np.float64)

                accs.add(acc)
                for i, (st, aa), in enumerate(zip(states, seq)):
                    ref[('ref', 'states')][(acc, i)] = float(st)
                    ref[('ref', 'seq')][(acc, i)] = aa

    logging.debug("loaded reference as <dict>; {}".format(ref_file.stem))
    return ref, accs


def strip_split(string: str) -> list:
    """
    remove newlines and split by tabs, raise an error for unexpected number of elements from split

    :param string: string to modify
    :return: splitted string
    """
    splitted = string.strip('\n').rsplit('\t')
    return splitted


def parse_prediction(predfile, reference_ids, label=None, decimals=3, threshold=0.5, normalize=False):
    """

    :param predfile:
    :type predfile:
    :param reference_ids:
    :type reference_ids:
    :param label:
    :type label:
    :param decimals:
    :type decimals:
    :param threshold:
    :type threshold:
    :return:
    :rtype:
    """
    predfile = Path(predfile)
    label = label if label is not None else predfile.stem
    pred = {(label, 'states'): {}, (label, 'scores'): {}}

    reference_ids = set(reference_ids)

    with open(predfile.resolve()) as fhandle:
        faiter = (x[1] for x in groupby(fhandle, lambda line: line[0] == ">"))
        for acc in faiter:
            acc = next(acc).strip()[1:]
            body = next(faiter)

            if acc in reference_ids:
                positions, _, scores, states = zip(*map(strip_split, body))
                # check truth-ish of a value instead of that of its container since zip(*map(strip_split, body))
                # returns list of empty strings when values are missing
                if states[0]:
                    states = np.array(states, dtype=np.float64)
                    # when scores are missing, use states as scores
                    scores = np.array(scores, dtype=np.float64) if scores[0] else states
                else:
                    if scores[0]:
                        scores = np.array(scores, dtype=np.float64)
                        # when states are missing, apply threshold to scores to generate states
                        # if threshold is not passed, default to 0.5
                        states = np.greater_equal(scores, threshold).astype(np.float64)
                    else:
                        continue

                if normalize is True:
                    # normalize in range [0, 1]
                    if np.min(scores) < 0 or np.max(scores) > 1:
                        scores = (scores - np.min(scores)) / np.ptp(scores)

                # round scores to the number of decimals passed (default 3)
                scores = scores.round(decimals)

                # reshape so that casting to pd.DataFrame is very quick
                for i, (st, sc), in enumerate(zip(states, scores)):
                    pred[(label, 'states')][(acc, i)] = st
                    pred[(label, 'scores')][(acc, i)] = sc

    logging.debug('loaded prediciton as <dict>; {}'.format(label))
    return pred


def parse_thresholds(thr_file):
    thresholds = None

    if thr_file.resolve(strict=True):
        try:
            thresholds = {}
            with open(thr_file) as f:
                for line in f:
                    pred, thr = line.strip().split()
                    thresholds[pred] = float(thr)
        except IndexError:
            logging.error('Threshold file was not formatted properly, default threshold will be estimated from scores')

    return thresholds
