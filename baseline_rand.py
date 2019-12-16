import argparse

# relative imports
import baseline_random


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
    #TODO: use refname to build output name
    ref, refname = baseline_random.get_reference(args.reference)
    baseline_random.baseline_random(ref.values[:, 0], basename=refname, outpath=args.outdir)
    baseline_random.baseline_shuffle_dataset(ref.values[:, 0], basename=refname, outpath=args.outdir)
    baseline_random.baseline_shuffle_targets(ref, basename=refname, outpath=args.outdir)
