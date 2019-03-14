import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import auc as calc_auc
# relative imports
from plotting.curve import parse_curve_file

def get_auc_table(curve_files, labels=None, cutoff=None):
    # init table
    auc_tab = list()

    # get curves from all curve files
    for fname in curve_files:
        for label, auc, x, y, _ in parse_curve_file(fname):

            # replace label (first col of curve file) with label from label's file
            label = labels.get(label, label)

            if cutoff is not None:
                # apply cutoff to x and y
                y = y[x <= cutoff]
                x = x[x <= cutoff]

                if x.size != 0 and y.size != 0:
                    # recalculate auc for coordinates under a specified cutoff
                    auc = calc_auc(x, y)
                    # rescale auc based on cutoff
                    auc /= cutoff
                else:
                    logging.warning("%s: no points in either x or y above cutoff", label)
                    auc = np.nan

            # generate table row
            auc_tab.append([label, auc])

    return auc_tab


def draw_auc_scatterplot(roc_curve_files, prc_curve_files, labels=None, cutoff=0.2, outbasename: str = 'auc-points'):
    # transform method_label default to empty dict to exploit .get default value
    labels = dict() if labels is None else labels

    # get auc tables for prc (x) and roc (y) curves
    prc = get_auc_table(prc_curve_files, labels, cutoff)
    roc = get_auc_table(roc_curve_files, labels, cutoff)

    # convert tables into pandas Dataframe, use label as index
    prc = pd.DataFrame(prc, columns=['label', 'PRC AUC']).set_index('label')
    roc = pd.DataFrame(roc, columns=['label', 'ROC AUC']).set_index('label')

    # join tables per method label (x and y become columns 1 and 2)
    # if labels don't match fill missing values with nan
    scatter_table = pd.concat([prc, roc], axis=1, sort=False)

    # draw scatterplot
    ax = scatter_table.plot.scatter(x='PRC AUC', y='ROC AUC')

    # annotate each point with its label
    for k, v in scatter_table.iterrows():
        ax.annotate(k, v)

    # add diagonal dotted line (incidentally set x and y limits)
    ax.plot([0, 0.5], [0, 0.5], color='k', linestyle='--')

    plt.savefig('{}.png'.format(outbasename), bbox_inches='tight', dpi=300)
