import os
import logging
import numpy as np
import copy
import argparse
from pathlib import Path
# relative imports
from caid import set_logger
from vectorized_metrics import bvaluation
from baseline_random import get_reference


def tokenize(seq, n=7):
    """
    Tokenize an iterable in a series of adjacent windows

    Yield a series of blobs of size n * 2 + 1 with the
    actual residue in the middle.

    For ending residues mirrored windows are added.

    :param n: n * 2 + 1 is the size of the window
    :type n: int
    :return: an iterator of blobs
    :rtype: iterator
    """
    lstates = len(seq)
    n = n if n < lstates else lstates - 1

    new_seq = copy.copy(seq[1:n + 1][::-1])  # Reversed start
    new_seq += copy.copy(seq)
    new_seq += copy.copy(seq[-n - 1:-1][::-1])  # Reversed end
    seq = new_seq

    for pos in range(n, len(seq) - n):
        yield seq[pos - n:pos + n + 1]


def get_scores_movingwindow(sequence):
    tokens = list(tokenize(sequence))
    return np.array(tokens).mean(axis=1)


def save_prediction(ref_file, acc, seq, scores, states):
    lseq = len(seq)

    scores = scores if scores is not None else [''] * lseq
    states = states if states is not None else [''] * lseq

    assert lseq == len(scores) == len(states), 'sequence, states and scores must have the same len'

    ref_file.write('>{}\n'.format(acc))
    for i, (aa, sc, st) in enumerate(zip(seq, scores, states), 1):
        ref_file.write('{}\t{}\t{}\t{}\n'.format(i, aa, '{:.3f}'.format(sc) if isinstance(sc, float) else '', st))


def naive_prediction(loaded_ref, outbase, invert):
    logging.info('building naive baseline from reverse input reference')
    fname = '{}.txt'.format(outbase)

    with open(fname, 'w') as fhandle:
        for acc, data in loaded_ref.groupby(level=0):

            if invert is True:
                states = [1 if r == 0 else 0 for r in data['ref']['states']]
            else:
                states = [0 if r == 0 else 1 for r in data['ref']['states']]

            # scores = get_scores_movingwindow(states)
            scores = [''] * len(states)
            save_prediction(fhandle, acc, data['ref']['seq'], scores, states)

    return fname


def build_args(options, predfiles, output_basename):
    optns = [options.reference] + predfiles + ['-l', options.log,
                                               '-ll', options.logLevel,
                                               '-o', output_basename,
                                               '--suffix', 'cons']
    if options.replaceUndefined is not None:
        optns.extend(['-r', options.replaceUndefined])

    return optns


def parse_args():
    parser = argparse.ArgumentParser(
        prog='caid-assess', description="CAID: Critical Assessment of Intrinsic Disorder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('reference',
                        help='reference file to which predictions are to be compared')
    parser.add_argument("--pRef", help='reference from which to derive naif predictions', required=True)

    parser.add_argument('-i', '--invert', default=False, action='store_true')
    parser.add_argument('-o', '--outdir', default='.')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    set_logger(args.log, args.logLevel)

    ref, refname = get_reference(args.reference)
    pred, predname = get_reference(args.pRef)

    predname = predname if refname != predname else "self"

    pred_fname = naive_prediction(pred, Path(args.outdir) / predname, args.invert)

    bvaluation(reference=args.reference, predictions=[pred_fname], outpath=args.outdir,
               run_tag="naive-{}".format(predname), dataset=True, target=True, bootstrap=True)
