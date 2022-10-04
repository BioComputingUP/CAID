import argparse

# relative imports
import baseline_random
from caid import set_logger


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('reference',
                        help='reference file to which predictions are to be compared')

    parser.add_argument('-o', '--outdir', default='.')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':
    args = parse_args()
    set_logger(args.log, args.logLevel)
    ref, refname = baseline_random.get_reference(args.reference)
    baseline_random.baseline_random(ref, n=100, basename=refname, outpath=args.outdir, target=True)
    baseline_random.baseline_shuffle_dataset(ref, n=100, basename=refname, outpath=args.outdir, target=True)
    baseline_random.baseline_shuffle_targets(ref, n=100, basename=refname, outpath=args.outdir, target=True)
    baseline_random.baseline_fixed_positive_fraction(ref, 0.354, n=100, basename=refname, outpath=args.outdir, target=True)  # id content in DisProt 7.0 dataset
