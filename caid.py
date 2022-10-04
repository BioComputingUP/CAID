# module imports
import sys
import logging
import argparse
from pathlib import Path
# relative imports
from vectorized_metrics import bvaluation

# Get path where this piece of code is
# SCRIPT_DIR = Path(getsourcefile(lambda : 0)).resolve()


def parse_args():
    parser = argparse.ArgumentParser(prog='caid-assess', description="CAID: Critical Assessment of Intrinsic Disorder")

    parser.add_argument('reference', help='reference file')

    parser.add_argument('predictions', help="directory containing prediction file(s)")

    parser.add_argument('-o', '--outputDir', default='.',
                        help='directory where the output will be written (default: cwd)')

    parser.add_argument('-b', '--labels', default=None, help='filename with labels')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args


# def parse_config(config_file):
#     config_parser = configparser.ConfigParser()
#     config_parser.optionxform = str
#     config_parser.read(config_file)
#
#     return config_parser


def set_logger(logfile, level):
    logging.basicConfig(level=level,
                        format='%(asctime)s | %(module)-13s | %(levelname)-8s | %(message)s',
                        stream=open(logfile) if logfile else sys.stderr)


if __name__ == '__main__':
    args = parse_args()
    set_logger(args.log, args.logLevel)
    bvaluation(args.reference, list(Path(args.predictions).glob('**/*')),
               outpath=args.outputDir, dataset=True, bootstrap=True, target=True)
