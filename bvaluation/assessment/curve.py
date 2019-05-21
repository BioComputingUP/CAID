import numpy as np
from sklearn import metrics
# relative imports
from bvaluation.assessment.metric import AreaUnderTheCurve


class Curve(object):
    """
    Evaluation curve. Aggregates scores from different evaluation metrics.

        Evaluation metrics are stored in x and y vectors, which are functions of the thresholds vector. An optimal
        threshold can be calculated as the highest Youden's J index. The area under the curve can be calculated from x
        and y vectors.
    """
    def __init__(self, name: str = None, label: str = None):
        """
        Initialize Curve object. Set kwargs as attributes.

        :param name: extended name of the curve
        :param label: label of the curve, used in repr method
        """
        self.x = None
        self.y = None
        self.thresholds = None
        self.youdensj = None
        self.name = name
        self.label = label
        self.auc = AreaUnderTheCurve(name='AUC_' + self.name, label='AUC_' + self.label)

    def cutoff_youdens_j(self):
        """
        Calculate and set threshold corresponding to highest Youden's J index.

        For each (x, y) coordinates coming from a threshold, calculate Youden's J. Set as optimal threshold the one
        associated to the highest Youden's J.
        """
        if self.x.size != 0 and self.x.size != 0 and self.thresholds.size != 0:
            j_scores = self.y - self.x
            j_ordered = sorted(zip(j_scores, self.thresholds))
            self.youdensj = j_ordered[-1][1]

    def auc_from_points(self):
        """
        Calculate and set the area under the curve defined by x and y vectors.

        Set AUC calculated from x, y vectors. Calls `sklearn.metrics.auc` function. Find documentation at
        https://scikit-learn.org/stable/modules/generated/sklearn.metrics.auc.html
        """
        self.auc.from_curve_points(self.x, self.y)


class ROC(Curve):
    """Receiving Operator Characteristic for binary classification.
    <br>x axis = FPR
    <br>y axis = TPR
    <br>More about ROC at https://en.wikipedia.org/wiki/Receiver_operating_characteristic"""
    def __init__(self, name=None, label=None):
        super(ROC, self).__init__(name, label)

    def calc_points(self, ytrue, yscore):
        if yscore.size != 0:
            self.x, self.y, self.thresholds = metrics.roc_curve(ytrue, yscore, pos_label=1)
            self.cutoff_youdens_j()
            self.auc_from_points()
        else:
            self.x = np.array([])
            self.y = np.array([])
            self.thresholds = np.array([])
            self.auc.amount = np.nan


class PrecisionRecallCurve(Curve):
    """
    Precision-Recall curve for binary classification.
    <br>x axis = Recall
    <br>y axis = Precision
    <br>More about ROC at https://en.wikipedia.org/wiki/Receiver_operating_characteristic
    """
    def __init__(self, name=None, label=None):
        super(PrecisionRecallCurve, self).__init__(name, label)
        self.fmax = None

    def calc_points(self, ytrue, yscore):
        if yscore.size != 0:
            self.y, self.x, self.thresholds = metrics.precision_recall_curve(ytrue,
                                                                             yscore,
                                                                             pos_label=1)
            self.cutoff_youdens_j()
            self.cutoff_fmax()
            self.auc_from_points()
        else:
            self.x = np.array([])
            self.y = np.array([])
            self.thresholds = np.array([])
            self.auc.amount = np.nan

    def cutoff_fmax(self):
        if self.x.size != 0 and self.x.size != 0 and self.thresholds.size != 0:
            f_scores = 2 * ((self.y * self.x)/(self.y + self.x))
            f_ordered = sorted(zip(f_scores, self.thresholds))
            self.fmax = f_ordered[-1][1]
