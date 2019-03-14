import os
import logging
import numpy as np
import argparse
# relative imports
from caid import parse_config, set_logger
from bvaluation.evaluate import parse_args as parse_eval_args, bvaluation, build_output_basename
from bvaluation.assessment.dataset import ReferencePool

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

PSEUDOCOUNT = .0000001

blosum62_freqs = {'S': 0.059, 'F': 0.044, 'V': 0.072, 'R': 0.051, 'K': 0.056, 'A': 0.078,
                  'G': 0.083, 'N': 0.041, 'M': 0.024, 'E': 0.059, 'D': 0.052, 'C': 0.024,
                  'Q': 0.034, 'H': 0.025, 'P': 0.043, 'Y': 0.034, 'W': 0.014, 'T': 0.055,
                  'I': 0.062, 'L': 0.092}


def save_prediction(ref_file, acc, seq, scores, states):
    lseq = len(seq)

    scores = scores if scores is not None else [''] * lseq
    states = states if states is not None else [''] * lseq

    assert lseq == len(scores) == len(states), 'sequence, states and scores must have the same len'

    ref_file.write('>{}\n'.format(acc))
    for i, (aa, sc, st) in enumerate(zip(seq, scores, states), 1):
        ref_file.write('{}\t{}\t{:.3f}\t{}\n'.format(i, aa, sc, st))


def js_divergence(freqs, bg_distr):
    """ Return the Jensen-Shannon Divergence for the column with the background
    distribution bg_distr. sim_matrix is ignored. JSD is the default method."""

    assert len(freqs) == len(bg_distr), 'len of freqs and background are different'
    freqs[freqs == 0] = PSEUDOCOUNT
    fc = freqs / freqs.sum()
    r = (fc * 0.5) + (bg_distr * 0.5)
    distance = (fc * np.log2(fc / r) + bg_distr * np.log2(bg_distr / r)).sum() / 2

    return distance


def conservation_based_prediction(loaded_ref, outbase, background_freqs):
    logging.info('building baseline distance from blosum62 frequencies')
    fname = '{}_cons.txt'.format(outbase)

    with open(fname, 'w') as fhandle:
        for acc, data in loaded_ref.items():
            scores = list()
            with open(os.path.join(pssm_dir, '{}.fasta'.format(acc))) as asn:
                prior_freqs = None
                for line in asn:
                    line = line.strip().split()
                    ll = len(line)

                    if ll == 40:
                        columns = line[:20]
                        prior_freqs = np.array([background_freqs[aa] for aa in columns])

                    if ll == 44:
                        assert prior_freqs is not None, 'Prior not defined, check pssm file format'

                        res_pos, res, *body, inf, wgh = line
                        freqs = np.fromiter(map(lambda e: float(e) / 100, body[20:]), dtype=float)
                        js_div = js_divergence(freqs, prior_freqs)
                        scores.append(1 - js_div)
                states = [1 if js >= 0.4 else 0 for js in scores]
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

    parser.add_argument('-r', '--replaceUndefined', choices=['0', '1'], default=None,
                        help='replace value for undefined positions (-) in reference. '
                             'By default not applied')
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
    pssm_dir = dict(conf.items('data_directories'))['pssm']
    # outdir = dict(conf.items('data_directories'))['baseline']
    outdir = ''
    ref = ReferencePool(os.path.abspath(args.reference),
                        undefined_replace_value=args.replaceUndefined)

    output_basename = build_output_basename(args.reference,
                                            args.replaceUndefined, outdir, ['a', 'b'])

    # negs = 'len-negs' if args.strictNegatives is False else 'str-negs'
    # output_basename = os.path.join(outdir, '{}_{}_{}'.format(
    #     negs,
    #     os.path.splitext(os.path.basename(args.reference))[0].replace('_', '-'),
    #     code
    # ))

    pred_fname = conservation_based_prediction(ref, output_basename, blosum62_freqs)
    bvaluation(parse_eval_args(build_args(args, [pred_fname], outdir)))
