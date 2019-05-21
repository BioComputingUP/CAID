import os
import sys
import argparse

# TODO: is this a solution? seems so
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from bvaluation.plotting.curve import draw_curves
from bvaluation.plotting.scatter import draw_auc_scatterplot
from bvaluation.plotting.heatmap import draw_per_instance_heatmap


def check_args_consistency(plots, roc_files, prc_files, instance_scores, **kwargs):
    """
    Check if all necessary options have been set to to run properly

    :param plots: type of desired plot
    :param roc_files: files for roc-related plots
    :param prc_files: files for prc-related plots
    :param instance_scores: file for per instance heatmap
    :param kwargs: all remaining keyword arguments passed but not subject to checks
    :return:
    """
    if 'roc' in plots and roc_files is None:
        raise ValueError("no 'roc_files' provided with plot: 'roc'")

    if 'prc' in plots and prc_files is None:
        raise ValueError("no 'prc_files' provided with plot: 'prc'")

    if 'auc' in plots:
        if prc_files is None:
            raise ValueError("no 'prc_files' provided with plot: 'auc'")
        if roc_files is None:
            raise ValueError("no 'roc_files' provided with plot: 'auc'")

    if 'pih' in plots and instance_scores is None:
        raise ValueError("no 'instance_scores' provided with plot: 'pih'")


def parse_args(arglist=None):
    """

    :param arglist:
    :return:
    """
    parser = argparse.ArgumentParser(
        prog='evaluate.py', description="Binary evaluation package",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # plot type
    parser.add_argument('plots', choices=['roc', 'prc', 'auc', 'pih'], nargs='+',
                        help='type of plot to be drawn. Multiple choices allowed')

    # plot input files
    parser.add_argument('-r', '--roc_files', default=None, nargs='+',
                        help='file containing ROC points. Multiple files allowed')
    parser.add_argument('-p', '--prc_files', default=None, nargs='+',
                        help='file containing PRC points. Multiple files allowed')
    parser.add_argument('-i', '--instance_scores', default=None,
                        help='file containing scores per each reference target (instance)')

    # output options
    parser.add_argument('-o', '--output_dir', help='directory where the output will be written',
                        default='.')
    parser.add_argument('-b', '--labels', help='text file with prediction file name (without path) '
                                               'and label desired in output organized in two '
                                               'columns')

    # log options
    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--log_level", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    parsed_args = parser.parse_args() if arglist is None else parser.parse_args(arglist)
    # check if all necessary options are present
    check_args_consistency(**vars(parsed_args))

    return parsed_args


def draw(plots, roc_files, prc_files, instance_scores, output_dir, labels, log, log_level):
    """
    Wrapper of drawing functions.

    Main level function used to call one or more drawing functions.

    :param plots: type of desired plot
    :param roc_files: files containing roc points
    :param prc_files: files containing prc points
    :param instance_scores: files containing scores for each target
    :param output_dir: path to output dir
    :param labels: file containing labels
    :param log: log file
    :param log_level: log level
    """
    check_args_consistency(**locals())

    if 'roc' in plots:
        roc_basename = os.path.join(output_dir, '_'.join(roc_files[0].split('_')[:2])) + '_roc'
        draw_curves(roc_files, outbasename=roc_basename)

    if 'prc' in plots:
        prc_basename = os.path.join(output_dir, '_'.join(prc_files[0].split('_')[:2])) + '_prc'
        draw_curves(prc_files, outbasename=prc_basename)

    if 'auc' in plots:
        auc_basename = os.path.join(output_dir, '_'.join(roc_files[0].split('_')[:2])) + '_auc-points'
        draw_auc_scatterplot(roc_files, prc_files, cutoff=0.2, outbasename=auc_basename)

    if 'pih' in plots:
        pih_basename = os.path.join(output_dir, '_'.join(instance_scores.split('_')[:2])) + '_heatmap'
        draw_per_instance_heatmap(instance_scores, outbasename=pih_basename)


if __name__ == "__main__":
    args = parse_args()
    draw(**vars(args))
