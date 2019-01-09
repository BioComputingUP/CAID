import os
import logging
import numpy as np
import pandas as pd
# relative imports
from caid import parse_config, parse_args, set_logger, load_names
from caid.dataset import ReferencePool
from bvaluation.evaluate import evaluate

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
    fname = '{}.txt'.format(outbase)

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


if __name__ == '__main__':
    args = parse_args(SCRIPT_DIR)
    conf = parse_config(args.conf)
    set_logger(args.log, args.logLevel)
    pssm_dir = dict(conf.items('prj_directories'))['pssm']
    outdir = dict(conf.items('prj_directories'))['baseline']
    names = load_names(os.path.join(SCRIPT_DIR, 'caid_names.txt'))

    code = 'cons'
    ref = ReferencePool(os.path.abspath(args.reference),
                        strict_negatives=args.strictNegatives)

    output_basename = os.path.join(outdir,
                                   os.path.splitext(
                                       os.path.basename(args.reference))[0]
                                   .replace('_', '-')) + '_{}'.format(code)

    pred_fname = conservation_based_prediction(ref, output_basename, blosum62_freqs)
    evaluation, roc_evaluation, coverage = evaluate(ref, pred_fname, code)

    metrics_table = pd.DataFrame()
    metrics_table_r = pd.DataFrame()
    perinstance_table = pd.DataFrame()
    roc_points = list()
    prc_points = list()

    if evaluation is not None:
        d = evaluation.get_scores_asdict()
        d.update({'Cov': coverage})
        metrics_table = pd.concat([metrics_table, pd.DataFrame(d, index=[code])],
                                  axis=0, sort=False)

        pis = pd.DataFrame(evaluation.get_perinstance_scores(insert=code))
        colnames = pis.iloc[0][2:]
        pis = pis.iloc[1:].rename(columns=colnames).set_index([0, 1]).astype(np.float)
        perinstance_table = perinstance_table.append(pis, sort=False)

    if roc_evaluation is not None:
        roc, prc = evaluation.get_curves_repr().split('\n')
        roc_points.append('{} {}'.format(code, roc))
        prc_points.append('{} {}'.format(code, prc))

        d = evaluation.get_scores_asdict()
        d.update({'Cov': coverage})
        metrics_table_r = pd.concat([metrics_table_r, pd.DataFrame(d, index=[code])],
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
