import re
import json
import logging
import argparse
import numpy as np
from pathlib import Path
from itertools import groupby

# relative imports
from caid import set_logger


def parse_args():
    parser = argparse.ArgumentParser(
        prog='caid-makeref', description='CAID reference builder',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('disprotJson', help='json file of disprot annotation (obtained from api)')

    parser.add_argument('-d', '--disorder', help='fasta file of disorder regions (obtained from api)', default=None)
    parser.add_argument('-i', '--interaction', help='fasta file of interaction partners (obtained from api)', default=None)
    parser.add_argument('-e', '--exclude', help='fasta file of entries to exclude (obtained from api)', default=None)
    parser.add_argument('-s', '--structure', nargs='+', help='fasta file of the structured regions (obtained from script)', default=None)
    parser.add_argument('-o', '--outdir', help='output directory', default='.')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument('-ll', '--logLevel', default="ERROR",
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args


def parse_fasta(fasta):
    logging.debug('parse fasta: {}'.format(fasta))
    with open(fasta) as f:
        faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
        for header in faiter:
            yield next(header)[1:].strip(), "".join(s.strip() for s in next(faiter))


def parse_json(json_file):
    logging.debug('parse json: {}'.format(json_file))
    with open(json_file) as j:
        return json.load(j)['data']


def build_sequence(header, out_seq, pos_seq=None, neg_seq=None):
    out_seq = np.fromiter(out_seq, dtype=np.dtype('U1'))

    if neg_seq is not None:
        neg_seq = np.fromiter(neg_seq, dtype=np.dtype('U1'))
        if neg_seq.size == out_seq.size:
            out_seq = np.where(neg_seq != '-', neg_seq, out_seq)
        else:
            logging.warning('{}: pool and negative sequences have different lengths: {} {}'.format(header, out_seq.size,
                                                                                                   neg_seq.size))
    if pos_seq is not None:
        pos_seq = np.fromiter(pos_seq, dtype=np.dtype('U1'))
        if pos_seq.size == out_seq.size:
            out_seq = np.where(pos_seq != '-', pos_seq, out_seq)
        else:
            logging.warning('{}: pool and positive sequences have different lengths: {} {}'.format(header, out_seq.size,
                                                                                                   pos_seq.size))

    return out_seq


def build_reference(disprot_ann, pool, positives, excluded, label_map, outfile, negatives=None, allow_negative_only=True):
    logging.info('building reference: {} from positives: {}; negatives: {}'.format(outfile, positives, negatives))

    positives = {re.search('DP[0-9]{5}', h).group(): s for h, s in parse_fasta(positives)}
    negatives = {re.search('DP[0-9]{5}', h).group(): s for h, s in parse_fasta(negatives)} if negatives else None
    excluded = {re.search('DP[0-9]{5}', h).group() for h, _ in parse_fasta(excluded)}

    # build dict of amino acid sequences from disprot json
    sequences = {e['disprot_id']: e['sequence'] for e in disprot_ann}
    names = {e['disprot_id']: e['name'] for e in disprot_ann}

    with open(outfile, 'w') as o:
        # get header and label sequence from pool fasta
        for header, seq in parse_fasta(pool):
            header = header.split("|")[1]

            # exclude viral polyproteins
            if "polyprotein" in names[header].lower() :
                logging.debug('{} excluded because polyprotein'.format(header))
                continue

            # exclude entries obtained from excluded fasta
            if header in excluded:
                logging.debug('{} excluded because in --exclude fasta'.format(header))
                continue

            seq = (label_map['-'] for _ in seq)
            pos = (label_map[l] if l != '-' else l for l in positives[header]) if header in positives else None
            neg = (label_map[l] if l != '-' else l for l in negatives[header]) if negatives is not None and header in negatives else None
            seq = build_sequence(header, seq, pos, neg)

            # consume generator to build string
            seq = ''.join(seq)

            # check if '-' is the only label in seq
            if set(seq) - {'-'}:
                if allow_negative_only is False and '1' not in set(seq):
                    continue
                o.write('>{}\n{}\n{}\n'.format(header, sequences[header], seq))


if __name__ == "__main__":
    args = parse_args()
    set_logger(args.log, args.logLevel)

    disprot_ann = parse_json(args.disprotJson)
    outdir = Path(args.outdir)

    build_reference(disprot_ann, pool=args.disorder, positives=args.disorder, excluded=args.exclude,
                    label_map={'-': '0', 'D': '1', 'S': '0'},
                    outfile=outdir / 'disprot-disorder.txt')

    for negative in args.structure:
        build_reference(disprot_ann, pool=args.disorder, positives=args.disorder, excluded=args.exclude,
                        label_map={'-': '-', 'D': '1', 'S': '0'},
                        outfile=outdir / 'disprot-disorder-{}.txt'.format(Path(negative).stem), negatives=negative)

        build_reference(disprot_ann, pool=args.disorder, positives=negative, excluded=args.exclude,
                        label_map={'-': '1', 'S': '0'},
                        outfile=outdir / '{}-reverse.txt'.format(Path(negative).stem))

    build_reference(disprot_ann, pool=args.interaction, positives=args.interaction, excluded=args.exclude,
                    label_map={'-': '0', 'I': '1'},
                    outfile=outdir / 'disprot-binding.txt', allow_negative_only=False)

    build_reference(disprot_ann, pool=args.disorder, positives=args.interaction, excluded=args.exclude,
                    label_map={'-': '0', 'I': '1'},
                    outfile=outdir / 'disprot-binding-all.txt')

    build_reference(disprot_ann, pool=args.interaction, positives=args.interaction, excluded=args.exclude,
                    label_map={'-': '-', 'I': '1', 'D': '0', 'S': '-'},
                    outfile=outdir / 'disprot-binding-disorder.txt', negatives=args.disorder, allow_negative_only=False)
