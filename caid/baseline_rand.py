import argparse
import os
from pathlib import Path
import multiprocessing as mp


from caid import set_logger
from vectorized_metrics import baseline_random


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('reference',
                        help='reference dir that contains the files to which predictions are to be compared')

    parser.add_argument('-o', '--outdir', default='.')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    arguments = parser.parse_args()
    return arguments


def generate_references(reference):
    print('Generating references for {}'.format(reference), flush=True)
    ref, refname = baseline_random.get_reference(reference)
    baseline_random.baseline_random(ref, n=100, basename=refname, outpath=args.outdir, target=True)
    baseline_random.baseline_shuffle_dataset(ref, n=100, basename=refname, outpath=args.outdir, target=True)
    # baseline_random.baseline_shuffle_targets(ref, n=100, basename=refname, outpath=args.outdir, target=True)
    # baseline_random.baseline_fixed_positive_fraction(ref, 0.354, n=100, basename=refname, outpath=args.outdir,
    #                                                  target=True)  # id content in DisProt 7.0 dataset
    print('Done with {}'.format(reference), flush=True)


if __name__ == '__main__':
    args = parse_args()
    set_logger(args.log, args.logLevel)
    references = list(Path(args.reference).glob('*.fasta'))

    os.makedirs(args.outdir, exist_ok=True)

    with mp.Pool(processes=mp.cpu_count() - 1) as p:
        p.map(generate_references, references)
