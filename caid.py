import argparse
import logging
import os
import sys
from pathlib import Path

from vectorized_metrics.vectorized_metrics import bvaluation


def parse_args():
    parser = argparse.ArgumentParser(
        prog='caid-assess', description="CAID: Critical Assessment of Intrinsic Disorder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('reference_file', help='Path to the reference file')

    parser.add_argument('predictions', help="directory containing prediction file(s)")

    parser.add_argument('-o', '--outputDir', default='.',
                        help='directory where the output will be written (default: cwd)')

    parser.add_argument('-b', '--labels', default=None, help='filename with labels')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="WARNING",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels >= choice will be displayed')

    args = parser.parse_args()
    return args


def set_logger(logfile, level):
    logging.basicConfig(level=level,
                        format='%(asctime)s | %(module)-13s | %(levelname)-8s | %(message)s',
                        stream=open(logfile) if logfile else sys.stderr)


if __name__ == '__main__':
    args = parse_args()
    set_logger(args.log, args.logLevel)
    pred_paths = list(Path(args.predictions).glob('*.caid'))
    reference_path = Path(args.reference_file)

    os.makedirs(args.outputDir, exist_ok=True)

    bvaluation(reference_path, pred_paths, outpath=args.outputDir, dataset=True, bootstrap=True, target=True,
               accs_to_read=None)

    print('Done')
