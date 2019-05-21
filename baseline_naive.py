import os
import logging
import numpy as np
import copy
import argparse
# relative imports
from caid import parse_config, set_logger
from bvaluation.evaluate import parse_args as parse_eval_args, bvaluation, build_output_basename
from bvaluation.assessment.dataset import ReferencePool

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


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


def naif_prediction(loaded_ref, outbase):
    logging.info('building naif baseline from reverse input reference')
    fname = '{}.txt'.format(outbase)

    with open(fname, 'w') as fhandle:
        for acc, data in loaded_ref.items():
            states = [0 if r == 0 else 1 for r in data['states']]
            scores = get_scores_movingwindow(states)

            save_prediction(fhandle, acc, data['seq'], scores, states)

    return fname


def build_args(options, predfiles, output_basename):
    optns = [options.reference] + predfiles + ['-l', options.log,
                                               '-ll', options.logLevel,
                                               '-o', output_basename,
                                               '--suffix', 'cons']
    if options.replaceUndefined is not None:
        optns.extend(['-r', options.replaceUndefined])

    return optns


def parse_args(wd):
    parser = argparse.ArgumentParser(
        prog='caid-assess', description="CAID: Critical Assessment of Intrinsic Disorder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('reference',
                        help='reference file to which predictions are to be compared')
    parser.add_argument('-p', "--pRef", default=None, help='reference from which to derive naif predictions, if not set, reference argument will be used')

    parser.add_argument('-r', '--replaceUndefined', choices=['0', '1'], default=None,
                        help='replace value for undefined positions (-) in reference. '
                             'By default not applied')
    parser.add_argument('-o', '--outdir', default='.')
    parser.add_argument('-c', '--conf', type=str,
                        default=os.path.join(wd, 'config.ini'),
                        help="path to an alternative configuration file.")

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args(SCRIPT_DIR)
    conf = parse_config(args.conf)
    set_logger(args.log, args.logLevel)
    # pssm_dir = dict(conf.items('data_directories'))['pssm']

    ref = ReferencePool(os.path.abspath(args.reference),
                        undefined_replace_value=args.replaceUndefined)

    if args.pRef is not None:
        suffix = 'naive-{}'.format(os.path.basename(os.path.splitext(args.pRef)[0]).split('_')[-1])
        pred = ReferencePool(os.path.abspath(args.pRef),
                             undefined_replace_value=args.replaceUndefined)
    else:
        suffix = "naive-self"
        pred = ref

    output_basename = build_output_basename(args.reference, args.outdir, ['b'], suffix=suffix)
    pred_fname = naif_prediction(pred, output_basename)

    bvaluation(reference=args.reference, prediction=[pred_fname], outdir=args.outdir, suffix=suffix,
               replace_undefined=args.replaceUndefined, log=args.log, log_level=args.logLevel)
