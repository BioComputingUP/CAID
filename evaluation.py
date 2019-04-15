# module imports
import os
import re
import logging
import warnings
import argparse
# relative imports
from caid import parse_config, set_logger
from bvaluation.evaluate import bvaluation

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Get path where this piece of code is
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def parse_args(wd: str):
    parser = argparse.ArgumentParser(
        prog='caid-assess', description="CAID: Critical Assessment of Intrinsic Disorder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('mode', choices=['disorder', 'binding', 'test'],
                        help='category for which to produce caid-uniformed prediction files')

    parser.add_argument('reference',
                        help='reference file to which predictions are to be compared')

    parser.add_argument('-r', '--replaceUndefined', choices=['0', '1'], default=None,
                        help='replace value for undefined positions (-) in reference. '
                             'By default not applied')

    parser.add_argument('-o', '--outputDir', help='directory where the output will be written',
                        default='.')
    parser.add_argument('-b', '--labels', default=None, help='filename with labels')
    # config file
    parser.add_argument('-c', '--conf', type=str,
                        default=os.path.join(wd, 'config.ini'),
                        help="path to an alternative configuration file.")
    parser.add_argument('-p', '--pattern', default=None)
    # log options
    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args


def build_args(options, predfiles):
    optns = [options.reference] + predfiles + ['-l', options.log,
                                               '-ll', options.logLevel,
                                               '-o', options.outputDir]
    if options.replaceUndefined is not None:
        optns.extend(['-r', options.replaceUndefined])
    if options.labels is not None:
        optns.extend(['-b', options.labels])
    return optns


def make_pfile_list(pred_dir, pattern, mode):
    if args.pattern is None:
        pattern = mode
        logging.warning('no --pattern passed, defaulting to mode: %s', mode)

    logging.debug('pattern: %s', pattern)
    r = re.compile(pattern)

    return [os.path.join(pred_dir, pred_fname) for pred_fname in filter(r.search,
                                                                        os.listdir(pred_dir))]


if __name__ == '__main__':
    args = parse_args(SCRIPT_DIR)
    conf = parse_config(args.conf)
    set_logger(args.log, args.logLevel)

    pred_dir = conf.get('data_directories', 'results')
    pred_files = make_pfile_list(pred_dir, args.pattern, args.mode)

    bvaluation(args.reference, pred_files, args.outputDir,
               replace_undefined=args.replaceUndefined, labels=args.labels,
               log=args.log, log_level=args.logLevel)

