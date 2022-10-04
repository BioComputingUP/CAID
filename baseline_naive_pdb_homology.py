import numpy as np
import pandas as pd
from make_references import parse_fasta
from vectorized_metrics import bvaluation
import argparse
from pathlib import Path


def get_missing_residues():
    print('loading missing residues')
    miss = {}
    for header, data in parse_fasta('../data/annotations/ss_dis.txt'):
        pdbid, chain, annotation = header.split(':')
        if annotation == 'disorder':
            miss['{}_{}'.format(pdbid.lower(), chain)] = data

    return miss


def get_caid_entries():
    print('loading caid disorder')
    entries = set()
    for header, _ in parse_fasta('../data/annotations/disprot-2018-11-disorder.fasta'):
        entries.add(header.split('|')[1])

    return entries


def get_sequences():
    print('loading sequences')
    seqs = {}
    for header, seq in parse_fasta('../data/annotations/disprot-2018-11-seq.fasta'):
        seqs[header.split()[0]] = seq

    return seqs


def get_pdbs_in_date_range(from_date, to_date):
    print('loading pdbs in depostion date range {} - {}'.format(from_date, to_date))
    pdb_entries = pd.read_csv('../data/annotations/pdb_entries.idx', skiprows=[0, 1], delimiter='\t', header=None)
    pdb_entries[2] = pd.to_datetime(pdb_entries[2][pdb_entries[2].str.match('\d{2}\/\d{2}\/\d{2}')])
    pdb_entries = pdb_entries.dropna(how='any', axis=0).set_index(2)
    return set(pdb_entries.loc[from_date:to_date][0])


def align_observed_residues(missing_residues, blast_output, caid_entries, sequences,  minid, maxid, pdbs_in_date_range):
    # slice missing residues to match alignmnent
    # insert subject gaps in missing residues
    # map observed residues on uniprot sequence
    # remove query gaps from observed residues

    print('aligning missing residues to uniprot sequence')
    caid_remote_pdbs = {}
    missing_from_missing_residues_file = set()
    missing_from_caid_entrires = set()

    with open(blast_output) as blast_out:
        for i, line in enumerate(blast_out):
            columns = line.strip().split('\t')
            if len(columns) == 14:
                qseqid, sseqid, pident, positive, length, mismatch, qstart, qend, sstart, send, evalue, bitscore, qseq, sseq = columns
                if qseqid in caid_entries and sseqid.split('_')[0].upper() in pdbs_in_date_range:

                    seqlen = len(sequences[qseqid])
                    pident, evalue, bitscore = [float(n) for n in [pident, evalue, bitscore]]
                    qstart, qend, sstart, send = [int(n) for n in [qstart, qend, sstart, send]]
                    qseq, sseq = [np.array(list(s)) for s in [qseq, sseq]]

                    if minid <= pident <= maxid :
                        qgaps = qseq == '-'
                        sgaps = sseq == '-'

                        sgaps_idx = np.where(sseq == '-')[0]
                        sgaps_insert_idx = sgaps_idx - np.arange(len(sgaps_idx))

                        if sseqid in missing_residues:
                            missing = np.array(list(missing_residues[sseqid]))[sstart - 1: send]
                            missing = np.insert(missing, sgaps_insert_idx, 'X')
                            try:
                                missing = missing[np.invert(qgaps)]
                            except:
                                print(qseqid, sseqid)
                                print('{} + {} = {} = {}'.format(len(qgaps), qend - qstart + 1, len(qgaps) + (qend - qstart + 1), len(qseq)))
                                print('{} + {} = {} = {}'.format(len(sgaps), send - sstart + 1, len(sgaps) + (send - sstart + 1), len(sseq)))
                                print(''.join(missing))
                                print(''.join(qgaps.astype(int).astype(str)))
                                print(''.join(sgaps.astype(int).astype(str)))

                            observed = np.invert(missing == 'X').astype(int).tolist()
                            observed_mapped_on_uniprot = [0] * seqlen
                            observed_mapped_on_uniprot[qstart - 1: qend] = observed
                            observed_mapped_on_uniprot = np.array(observed_mapped_on_uniprot)

                            if caid_remote_pdbs.get(qseqid, {}).get(sseqid, {}).get('pident', 0) <= pident:
                                caid_remote_pdbs.setdefault(qseqid, {})[sseqid] = {
                                    'pident': pident,
                                    'observed': observed_mapped_on_uniprot
                                }

                        else:
                            missing_from_missing_residues_file.add(sseqid)
                else:
                    missing_from_caid_entrires.add(qseqid)

    print('{} missing from caid entries'.format(len(missing_from_caid_entrires)))
    print('{} missing from ss_dis.txt'.format(len(missing_from_missing_residues_file)))
    return caid_remote_pdbs


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('reference',
                        help='reference file to which predictions are to be compared')
    parser.add_argument('blast',
                        help='tabular blast output')

    parser.add_argument('-o', '--outdir', default='.')
    parser.add_argument('-n', '--minid', default=20, type=int)
    parser.add_argument('-x', '--maxid', default=30, type=int)
    parser.add_argument('-f', '--fromDate', default=None, help='YYYY-MM-DD')
    parser.add_argument('-t', '--toDate', default=None, help='YYYY-MM-DD')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    caid_entries = get_caid_entries()
    sequences = get_sequences()
    included_pdbs = get_pdbs_in_date_range(args.fromDate, args.toDate)
    caid_observed_mapped = align_observed_residues(get_missing_residues(), args.blast, caid_entries,
                                                   sequences, args.minid, args.maxid, included_pdbs)

    date_range = 'from-{}-to-{}'.format(
        args.fromDate if args.fromDate is not None else 'any',
        args.toDate if args.toDate is not None else 'any')

    baseline_basename = 'pdb-remote{}-{}-{}'.format(args.minid, args.maxid, date_range)
    baseline_filepath = Path(args.outdir) / (baseline_basename + '.txt')
    
    with open(baseline_filepath, 'w') as baseline_file:
        for disprot_id in caid_observed_mapped:
            observed_mapped_merged = np.any(np.vstack([caid_observed_mapped[disprot_id][a]['observed'] for a in caid_observed_mapped[disprot_id]]), axis=0)
            reference = np.invert(observed_mapped_merged).astype(int).astype(str)
            
            baseline_file.write('>{}\n{}\n'.format(
                disprot_id,
                '\n'.join(['{}\t{}\t\t{}'.format(idx, aa, state) for idx, aa, state in 
                           zip(range(reference.size), sequences[disprot_id], reference.tolist())])
            ))

    bvaluation(reference=args.reference, predictions=[baseline_filepath], outpath=args.outdir ,
               run_tag="naive-{}".format(baseline_basename, date_range), dataset=True, target=True, bootstrap=True)
