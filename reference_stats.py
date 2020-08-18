import re
import requests
import json
import logging
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import groupby

# relative imports
from caid import set_logger

RESIDUE_MAP = {'1': 1, '0': 0, '-': np.nan}


def parse_args():
    parser = argparse.ArgumentParser(
        prog='caid-makeref', description='CAID reference builder',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('disprotJson', help='json file of disprot annotation (obtained from api)')

    parser.add_argument('ref_files', nargs='+', help='reference file(s) built with make_references.py')
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

def get_species(taxonomy):
    s = None
    if taxonomy:
        s = taxonomy[-1]
    return s


def get_taxon(taxonomy, acc):
    try:
        if taxonomy is None:
            txtup = requests.get('https://www.uniprot.org/uniprot/{}.txt'.format(acc)).text

            if not txtup:

                acc = (re.search('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}',
                                requests.get('https://www.uniprot.org/uniparc/{}.xml'.format(acc)).text).group(0))
                txtup = requests.get('https://www.uniprot.org/uniprot/{}.txt'.format(acc)).text

            taxonomy = next(line for line in
                            requests.get('https://www.uniprot.org/uniprot/{}.txt'.format(acc))
                            .text.split('\n') if line.startswith('OC')).split(';')[0][2:].lstrip()
        else:
            taxonomy = taxonomy[0]

    except:
        taxonomy = None
        if acc == 'A0A0P6V0V4':
            taxonomy = 'Bacteria'

    return taxonomy


def calc_ref_stats(disprot_ann, ref_file):
    logging.info('calculating stats from: {} '.format(ref_file.stem))

    target_table = {}
    region_table = {}
    # get header and label sequence from pool fasta
    for header, seq in parse_fasta(ref_file):
        _, states = seq[:len(seq)//2], seq[len(seq)//2:]
        regions = [(grouper, list(group)) for grouper, group in groupby(states)]
        taxon = get_taxon(disprot_ann[header]['taxonomy'], disprot_ann[header]['acc'])
        species = get_species(disprot_ann[header]['taxonomy'])
        print(disprot_ann[header]['acc'], disprot_ann[header]['taxonomy'])
        exit()
        # if taxon is None:
        #     continue

        target_table.update({
            (ref_file.stem, header): {
                ('residues', 'positive'): states.count('1'),
                ('residues', 'negative'): states.count('0'),
                ('residues', 'undefined'): states.count('-'),
                ('residues', 'total'): len(states),
                ('regions', 'positive'): sum(1 for grouper, _ in regions if grouper == '1'),
                ('regions', 'negative'): sum(1 for grouper, _ in regions if grouper == '0'),
                ('regions', 'undefined'): sum(1 for grouper, _ in regions if grouper == '-'),
                ('regions', 'total'): sum(1 for _ in regions),
                ('data', 'taxon'): taxon,
                ('data', 'species'): species,
                ('data', 'idc'): states.count('1') / (len(states) - states.count('-'))
            }
        })

        for i, (label, reg) in enumerate(regions, 1):
            region_table.update({
                (ref_file.stem, header, i): {
                    'type': label,
                    'length': len(list(reg)),
                    'taxon': taxon,
                }
            })

    return target_table, region_table


if __name__ == "__main__":
    args = parse_args()
    set_logger(args.log, args.logLevel)

    disprot_ann = {e['disprot_id']: e for e in parse_json(args.disprotJson)}
    outdir = Path(args.outdir)

    overall_target_table = {}
    overall_region_table = {}

    for ref_file in args.ref_files:
        ref_file = Path(ref_file)
        t, r = calc_ref_stats(disprot_ann, ref_file)

        overall_target_table.update(t)
        overall_region_table.update(r)

    pd.DataFrame.from_dict(overall_target_table, orient='index').to_csv(outdir / 'references-stats.target.csv')
    pd.DataFrame.from_dict(overall_region_table, orient='index').to_csv(outdir / 'references-stats.region.csv')

