import copy
import logging
import numpy as np
from math import sqrt
from sklearn import metrics


class Metric(object):
    name = None
    label = None

    def __init__(self, robust=False, name=None, label=None):
        self.name = name if name is not None else self.name
        self.label = label if label is not None else self.label
        self.robust = robust
        self.pseudocount = 0.000001
        self.amount = None
        self.average = None
        self.std = None

    def set_avg_and_std(self, amount_array):
        self.average = np.nanmean(amount_array)
        self.std = np.nanstd(amount_array)


class AreaUnderTheCurve(object):
    def __init__(self, name=None, label=None):
        self.amount = None
        self.name = name
        self.label = label
        self.average = None
        self.std = None

    def from_curve_points(self, x, y):
        if x.size != 0 and y.size != 0:
            x = copy.copy(x)
            y = copy.copy(y)
            self.amount = metrics.auc(sorted(x), sorted(y))
        else:
            self.amount = np.nan


class BalancedAccuracy(Metric):
    name = 'Balanced accuracy'
    label = 'BAc'

    def __init__(self, robust=False, name=None, label=None):
        super(BalancedAccuracy, self).__init__(robust, name, label)

    def from_states(self, ytrue, ypred):
        self.amount = metrics.balanced_accuracy_score(ytrue, ypred)

    def from_confusion_matrix(self, conf_mat):
        tn, fp, fn, tp = conf_mat.ravel()
        p = tp + fn
        n = tn + fp

        if self.robust is True:
            p = p if p != 0 else self.pseudocount
            n = n if n != 0 else self.pseudocount

        try:
            self.amount = (tp / p + tn / n) / 2
        except ZeroDivisionError:
            self.amount = np.nan


class MCC(Metric):
    name = 'Matthew\'s Correlation Coefficient'
    label = 'MCC'

    def __init__(self, robust=False, name=None, label=None):
        super(MCC, self).__init__(robust, name, label)

    def from_states(self, ytrue, ypred):
        self.amount = metrics.matthews_corrcoef(ytrue, ypred)

    def from_confusion_matrix(self, conf_mat):
        tn, fp, fn, tp = conf_mat.ravel()
        denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)

        if self.robust is True:
            denom = denom if denom > 0 else self.pseudocount

        try:
            self.amount = (tp * tn - fp * fn) / sqrt(denom)
        except ZeroDivisionError:
            self.amount = np.nan


class Precision(Metric):
    name = 'Precision'
    alternative_names = ['Positive Predictive Value']
    label = 'Pre'

    def __init__(self, robust=False, name=None, label=None):
        super(Precision, self).__init__(robust, name, label)

    def from_states(self, ytrue, ypred):
        self.amount = metrics.precision_score(ytrue, ypred)

    def from_confusion_matrix(self, conf_mat):
        tn, fp, fn, tp = conf_mat.ravel()
        pos_assignments = tp + fp

        if self.robust is True:
            pos_assignments = pos_assignments if pos_assignments != 0 else self.pseudocount

        try:
            self.amount = tp / pos_assignments
        except ZeroDivisionError:
            self.amount = np.nan


class Recall(Metric):
    name = 'Recall'
    alternative_names = ['Sensitivity', 'Hit rate', 'True Positive Rate']
    label = 'Rec'

    def __init__(self, robust=False, name=None, label=None):
        super(Recall, self).__init__(robust, name, label)

    def from_states(self, ytrue, ypred):
        self.amount = metrics.recall_score(ytrue, ypred)

    def from_confusion_matrix(self, conf_mat):
        tn, fp, fn, tp = conf_mat.ravel()
        p = tp + fn

        if self.robust is True:
            p = p if p != 0 else self.pseudocount

        try:
            self.amount = tp / p
        except ZeroDivisionError:
            self.amount = np.nan


class FScore(Metric):
    name = 'F1 score'
    label = 'F1s'

    def __init__(self, robust=False, name=None, label=None):
        super(FScore, self).__init__(robust, name, label)

    def from_states(self, ytrue, ypred):
        self.amount = metrics.f1_score(ytrue, ypred)

    def from_confusion_matrix(self, conf_mat):
        tn, fp, fn, tp = conf_mat.ravel()
        denom = 2 * tp + fp + fn

        if self.robust is True:
            denom = denom if denom != 0 else self.pseudocount

        try:
            self.amount = (2 * tp) / denom
        except ZeroDivisionError:
            self.amount = np.nan

    def from_precision_recall(self, precision, recall):
        self.amount = 2 * (precision * recall) / (precision + recall)


class FalsePositiveRate(Metric):
    name = 'False Positive Rate'
    alternative_names = ['Fall-out']
    label = 'FPR'

    def __init__(self, robust=False, name=None, label=None):
        super(FalsePositiveRate, self).__init__(robust, name, label)

    def from_confusion_matrix(self, conf_mat):
        tn, fp, fn, tp = conf_mat.ravel()
        n = tn + fp

        if self.robust is True:
            n = n if n != 0 else self.pseudocount

        try:
            self.amount = fp / n
        except ZeroDivisionError:
            self.amount = np.nan

    def from_states(self, ytrue, ypred):
        self.from_confusion_matrix(metrics.confusion_matrix(ytrue, ypred))



