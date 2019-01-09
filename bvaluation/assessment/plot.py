import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class RocPlot(object):
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.codes = list()
        self.points = list()

    def add_roc(self, fpr, tpr, code=None, predname=None, auc=None, coverage=None, style='-', lw=1):
        lbl = ''
        if predname != 'Conservation':
            lw = 3
        lbl += predname if predname else ''
        lbl += 'AUC {:.2f}' if auc else ''
        lbl += 'COV {:.2f}' if coverage else ''
        lbl = lbl if lbl else None

        self.ax.plot(fpr, tpr, lw=lw, label=lbl, linestyle=style)

        if code:
            self.codes.append(code)

        self.points.append('{} {} {} {} {}'.format(
            predname if predname else 'x',
            code if code else 'x', auc, coverage,
            ' '.join('{},{}'.format(f, t) for f, t in zip(fpr, tpr)))
        )

    def save_figure(self, basename='', hide_after_first=5):
        self.ax.plot([0, 1], [0, 1], color='lightgrey', lw=1, linestyle='--')
        handles, labels = self.ax.get_legend_handles_labels()

        if len(self.codes) == len(handles):
            labels, handles, codes = zip(*sorted(
                zip(labels, handles, self.codes),
                key=lambda t: float(t[0].split('AUC ')[1].split()[0]),
                reverse=True))

            newlables = list()
            for i, (label, code) in enumerate(zip(labels, codes)):
                if i >= hide_after_first and 'Conservation' not in label:
                    label = label.split()
                    prefix = '*' if label[0][0] == '*' else ''
                    label[0] = prefix + code
                    label = ' '.join(label)
                newlables.append(label)
            labels = newlables

        else:
            logging.warning('different number of codes and labels. Resorting to bare labels')

            labels, handles = zip(*sorted(zip(labels, handles),
                                          key=lambda t: float(t[0].split('AUC ')[1].split()[0]),
                                          reverse=True))

        # font = font_manager.FontProperties(family='Andale Mono')
        self.ax.legend(handles, labels, loc="lower right", prop={'size': 8})
        self.ax.set_xlabel('FPR')
        self.ax.set_ylabel('TPR')

        plt.tight_layout()
        fname = '{}_rocs.png'.format(basename)
        logging.info('saving ROC plot to: %s', fname)
        plt.savefig(fname, dpi=300, layout='tight')
        plt.close()

    def save_points(self, basename):
        fname = '{}_rocpoints.txt'.format(basename)
        with open(fname, 'w') as f:
            f.write('\n'.join(self.points))


class PrecisionRecallCurve(object):
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.codes = list()
        self.points = list()

    def add_roc(self, precision, recall, code=None, predname=None,
                auc=None, coverage=None, style='-', lw=1):
        lbl = ''
        if predname != 'Conservation':
            lw = 3
        lbl += predname if predname else ''
        lbl += 'AUC {:.2f}' if auc else ''
        lbl += 'COV {:.2f}' if coverage else ''
        lbl = lbl if lbl else None

        self.ax.plot(precision, recall, lw=lw, label=lbl, linestyle=style)

        if code:
            self.codes.append(code)

        self.points.append('{} {} {}'.format(
            predname if predname else 'x',
            code if code else 'x',
            ' '.join('{},{}'.format(f, t) for f, t in zip(precision, recall)))
        )

    def save_figure(self, basename='', hide_after_first=5):
        self.ax.plot([0, 1], [0, 1], color='lightgrey', lw=1, linestyle='--')
        handles, labels = self.ax.get_legend_handles_labels()

        if len(self.codes) == len(handles):
            labels, handles, codes = zip(*sorted(
                zip(labels, handles, self.codes),
                key=lambda t: float(t[0].split('AUC ')[1].split()[0]),
                reverse=True))

            newlables = list()
            for i, (label, code) in enumerate(zip(labels, codes)):
                if i >= hide_after_first and 'Conservation' not in label:
                    label = label.split()
                    prefix = '*' if label[0][0] == '*' else ''
                    label[0] = prefix + code
                    label = ' '.join(label)
                newlables.append(label)
            labels = newlables

        else:
            logging.warning('different number of codes and labels. Resorting to bare labels')

            labels, handles = zip(*sorted(zip(labels, handles),
                                          key=lambda t: float(t[0].split('AUC ')[1].split()[0]),
                                          reverse=True))

        # font = font_manager.FontProperties(family='Andale Mono')
        self.ax.legend(handles, labels, loc="lower right", prop={'size': 8})
        self.ax.set_xlabel('Precision')
        self.ax.set_ylabel('Recall')

        plt.tight_layout()
        fname = '{}_prcs.png'.format(basename)
        logging.info('saving precision-recall curve plot to: %s', fname)
        plt.savefig(fname, dpi=300, layout='tight')
        plt.close()

    def save_points(self, basename):
        fname = '{}_prpoints.txt'.format(basename)
        with open(fname, 'w') as f:
            f.write('\n'.join(self.points))


def stack_as_barplot(s, output_basename):
    averaged_stack = np.nanmean(s, axis=1)
    pd.DataFrame(data=averaged_stack, columns=['acc', 'mcc']).plot.hist(alpha=0.5, bins=20)
    plt.tight_layout()
    plt.savefig(output_basename + '_omd.png', dpi=300)
