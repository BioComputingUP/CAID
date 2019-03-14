# module imports
import os
import warnings
import argparse
import numpy as np
import matplotlib.pyplot as plt
import asyncio
# relative imports
from caid import parse_config, set_logger
from bvaluation.evaluate import parse_args as parse_eval_args, bvaluation

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Get path where this piece of code is
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def parse_args(wd):
    parser = argparse.ArgumentParser(
        prog='caid-assess', description="CAID: Critical Assessment of Intrinsic Disorder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--inputDir', required=True)
    parser.add_argument('-b', '--baselineDir', default=None)
    parser.add_argument('-c', '--conf', type=str,
                        default=os.path.join(wd, 'config.ini'),
                        help="path to an alternative configuration file.")
    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args


def coroutine(func):
    def start(*args,**kwargs):
        corout = func(*args,**kwargs)
        next(corout)
        return corout
    return start


def parse_curve_file(curve_file):
    with open(curve_file) as f:
        for line in f:
            code, _, auc, *points = line.split()
            auc = float(auc)
            x, y, thr = zip(*map(lambda s: map(float, s.split(',')), points))

            yield code, auc, np.array(x), np.array(y), np.array(thr)


def draw_curve(fname, basename, baselines, result_label):
        # initialize a figure and an ax
        fig, ax = plt.subplots()
        # plot diagonal
        ax.plot([0, 1], [0, 1], linestyle='--', color='k')
        # plot all curves
        for label, auc, x, y, thr in parse_curve_file(fname):
            ax.plot(x, y, label='{} AUC: {:.2f}'.format(label, auc))

        if baselines is not None:
            baseline_roc = '{}_{}_{}.txt'.format(basename, 'cons', result_label)

            if baseline_roc in baselines:
                baseline_roc = os.path.join(args.baselineDir, baseline_roc)

                assert os.path.isfile(baseline_roc) is True

                # plot baseline curve from its own file
                for curve_data in parse_curve_file(baseline_roc):
                    label, auc_roc, x, y, thr = curve_data
                    ax.plot(x, y, label='{} AUC: {:.2f}'.format(label, auc_roc))
        ax.legend()
        plt.show()


if __name__ == '__main__':
    args = parse_args(SCRIPT_DIR)
    conf = parse_config(args.conf)
    set_logger(args.log, args.logLevel)

    baselines_ids = ['chain', 'cons', 'chunk']

    for fname in os.listdir(args.inputDir):
        basename = '_'.join(fname.split('_')[:-1])
        fname = os.path.join(os.path.abspath(args.inputDir), fname)

        # condition allows everything to work even if baselines are in the same folder as predictors
        if basename.split('_')[-1] not in baselines_ids:
            if args.baselineDir is not None:
                baselines = [f for f in os.listdir(args.baselineDir) if basename in f
                             and any(e in f for e in baselines_ids)]
            else:
                baselines = None

            if 'rocPoints' in fname:
                draw_curve(fname, basename, baselines, 'rocPoints')
            if 'prcPoints' in fname:
                draw_curve(fname, basename, baselines, 'prcPoints')
