# module imports
import os
import re
import logging
import warnings
# relative imports
from caid import parse_args, parse_config, set_logger
from bvaluation.evaluate import parse_args as parse_eval_args, main

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Get path where this piece of code is
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def build_args(options, predfiles):
    optns = [options.reference] + predfiles + ['-l', options.log,
                                               '-ll', options.logLevel,
                                               '-o', options.outputDir,
                                               '--prefix', options.mode]
    if options.strictNegatives is True:
        optns.append('-s')

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

    pred_dir = conf.get('prj_directories', 'results')
    pred_files = make_pfile_list(pred_dir, args.pattern, args.mode)

    main(parse_eval_args(build_args(args, pred_files)))
