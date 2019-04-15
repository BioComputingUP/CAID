import numpy as np
import matplotlib.pyplot as plt

from typing import Generator


def parse_curve_file(curve_file: str) -> Generator:
    """
    Lazy load curve points for a prediction method

    :param curve_file: file containing curve points
    :return: generator yielding (code, auc, x, y, threshold)
    """
    with open(curve_file) as f:
        for line in f:
            code, _, auc, *points = line.split()
            auc = float(auc)
            x, y, thr = zip(*map(lambda s: map(float, s.split(',')), points))

            yield code, auc, np.array(x), np.array(y), np.array(thr)


def draw_curves(curve_files: list, method_label: dict = None, outbasename: str = 'curve'):
    """
    Draw ROC or Precision-Recall curve

    :param curve_files: file(s) containing curve points
    :param method_label: label of the prediction method
    :param outbasename: basename of the output file
    """
    # transform method_label default to empty dict to exploit .get default value
    method_label = dict() if method_label is None else method_label

    # initialize a figure and an ax
    fig, ax = plt.subplots()

    # plot diagonal
    ax.plot([0, 1], [0, 1], linestyle='--', color='k')

    # plot all curves from all curve files
    for fname in curve_files:
        for label, auc, x, y, thr in parse_curve_file(fname):

            # replace label (first col of curve file) with label from label's file
            label = method_label.get(label, label)

            # draw curve
            ax.plot(x, y, label='{} AUC: {:.2f}'.format(label, auc))

    ax.legend()

    plt.savefig('{}.png'.format(outbasename), bbox_inches='tight', dpi=300)
