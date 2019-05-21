import inspect
import numpy as np
from sklearn import metrics
# relative imports
from bvaluation.assessment.metric import Precision, Recall, FScore, MCC, FalsePositiveRate, BalancedAccuracy, TruePositive, TrueNegative, FalseNegative, FalsePositive
from bvaluation.assessment.curve import ROC, PrecisionRecallCurve


class Curves(object):
    def __init__(self):
        self.roc = ROC(name='Receiving Operator Characteristic', label='ROC')
        self.prc = PrecisionRecallCurve(name='Precision Recall Curve', label='PRC')
        self._is_calculated = False

    def calc_curves(self, ytrue, yscore):
        self.roc.calc_points(ytrue, yscore)
        self.prc.calc_points(ytrue, yscore)
        self._is_calculated = True

    def get_members(self):
        variables = inspect.getmembers(self, lambda a: not inspect.ismethod(a))
        return [
            (a[1].auc.name, a[1].auc.label, a[1].auc.amount, a[1].name,
             a[1].label, a[1].x, a[1].y, a[1].thresholds)
            for a in variables if a[0][0] != '_']

    def __repr__(self):
        return '\n'.join('{} {} {}\n{}\n{}'.format(*m[:3], ','.join(np.around(m[3], 3).astype(str)),
                                                   ','.join(np.around(m[4], 3).astype(str))) for m
                         in self.get_members())

    def __str__(self):
        return '\n'.join('{:<37} {:<7} {:.2f}'.format(*m[:3]) for m in self.get_members())


class Scores(object):
    def __init__(self):
        self.mcc = MCC()
        self.fpr = FalsePositiveRate()
        self.recall = Recall() # sensitivity
        self.fscore = FScore()
        self.bal_acc = BalancedAccuracy()
        self.recall_n = Recall(label=Recall.label + '_n') # specificity
        self.fscore_n = FScore(label=FScore.label + '_n')
        self.precision = Precision()
        self.precision_n = Precision(label=Precision.label + '_n')
        self.tp = TruePositive()
        self.fn = FalseNegative()
        self.tn = TrueNegative()
        self.fp = FalsePositive()
        self._is_calculated = False

    def states_based_metrics(self, confusion_matrix=None, ytrue=None, ypred=None):
        if confusion_matrix is not None:

            self.recall.from_confusion_matrix(confusion_matrix)
            self.bal_acc.from_confusion_matrix(confusion_matrix)
            self.precision.from_confusion_matrix(confusion_matrix)
            self.mcc.from_confusion_matrix(confusion_matrix)
            self.fpr.from_confusion_matrix(confusion_matrix)
            self.fscore.from_confusion_matrix(confusion_matrix)

            self.tn.from_confusion_matrix(confusion_matrix)
            self.fp.from_confusion_matrix(confusion_matrix)
            self.fn.from_confusion_matrix(confusion_matrix)
            self.tp.from_confusion_matrix(confusion_matrix)

            confusion_matrix_n = np.flip(confusion_matrix)
            self.precision_n.from_confusion_matrix(confusion_matrix_n)
            self.fscore_n.from_confusion_matrix(confusion_matrix_n)
            self.recall_n.from_confusion_matrix(confusion_matrix_n)

        elif ytrue is not None and ypred is not None:

            # self.recall.from_states(ytrue, ypred)
            self.bal_acc.from_states(ytrue, ypred)
            # self.precision.from_states(ytrue, ypred)
            self.mcc.from_states(ytrue, ypred)
            self.fpr.from_states(ytrue, ypred)
            self.tn.from_states(ytrue, ypred)
            self.fp.from_states(ytrue, ypred)
            self.fn.from_states(ytrue, ypred)
            self.tp.from_states(ytrue, ypred)
            # self.fscore.from_states(ytrue, ypred)
            report = metrics.classification_report(ytrue, ypred, output_dict=True)
            preport = report.get('1.0', {})
            nreport = report.get('0.0', {})
            self.recall.amount = preport.get('recall', np.nan)
            self.precision.amount = preport.get('precision', np.nan)
            self.fscore.amount = preport.get('f1-score', np.nan)

            self.recall_n.amount = nreport.get('recall', np.nan)
            self.precision_n.amount = nreport.get('precision', np.nan)
            self.fscore_n.amount = nreport.get('f1-score', np.nan)



        self._is_calculated = True

    def get_members(self):
        variables = inspect.getmembers(self, lambda a: not inspect.ismethod(a))
        return sorted([(a[1].name, a[1].label, a[1].amount) for a in variables if a[0][0] != '_'], key=lambda t: len(t[1]))

    def __repr__(self):
        return '\n'.join('{} {} {}'.format(*m) for m in self.get_members())

    def __str__(self):
        return '\n'.join('{:<37} {:<7} {:.2f}'.format(*m) for m in self.get_members())


class AverageInstanceScores(Scores):
    def __init__(self):
        super(AverageInstanceScores, self).__init__()

    def from_iterable(self, ytrue_list, ypred_list):
        bal_acc_array = list()
        precision_array = list()
        recall_array = list()
        mcc_array = list()
        fpr_array = list()
        fscore_array = list()
        precision_n_array = list()
        recall_n_array = list()
        fscore_n_array = list()
        tn_array = list()
        fp_array = list()
        fn_array = list()
        tp_array = list()


        ypred_list = ypred_list if ypred_list is not None else [np.array([])] * len(ytrue_list)


        for ytrue, ypred in zip(ytrue_list, ypred_list):
            confusion_matrix = metrics.confusion_matrix(ytrue, ypred,
                                                        labels=(0, 1)).astype(np.float)
            self.states_based_metrics(
                # confusion_matrix=confusion_matrix,
                ytrue=ytrue,
                ypred=ypred
            )
            precision_array.append(self.precision.amount)
            bal_acc_array.append(self.bal_acc.amount)
            recall_array.append(self.recall.amount)
            mcc_array.append(self.mcc.amount)
            fpr_array.append(self.fpr.amount)
            fscore_array.append(self.fscore.amount)
            precision_n_array.append(self.precision_n.amount)
            recall_n_array.append(self.recall_n.amount)
            fscore_n_array.append(self.fscore_n.amount)

            tn, fp, fn, tp = confusion_matrix.ravel()
            tn_array.append(tn)
            fp_array.append(fp)
            fn_array.append(fn)
            tp_array.append(tp)

        precision_array = np.array(precision_array, dtype=np.float)
        bal_acc_array = np.array(bal_acc_array, dtype=np.float)
        recall_array = np.array(recall_array, dtype=np.float)
        mcc_array = np.array(mcc_array, dtype=np.float)
        fpr_array = np.array(fpr_array, dtype=np.float)
        fscore_array = np.array(fscore_array, dtype=np.float)
        precision_n_array = np.array(precision_n_array, dtype=np.float)
        recall_n_array = np.array(recall_n_array, dtype=np.float)
        fscore_n_array = np.array(fscore_n_array, dtype=np.float)
        tn_array = np.array(tn_array, dtype=int)
        fp_array = np.array(fp_array, dtype=int)
        fn_array = np.array(fn_array, dtype=int)
        tp_array = np.array(tp_array, dtype=int)

        self.precision.set_avg_and_std(precision_array)
        self.bal_acc.set_avg_and_std(bal_acc_array)
        self.recall.set_avg_and_std(recall_array)
        self.mcc.set_avg_and_std(mcc_array)
        self.fpr.set_avg_and_std(fpr_array)
        self.fscore.set_avg_and_std(fscore_array)
        self.precision_n.set_avg_and_std(precision_n_array)
        self.recall_n.set_avg_and_std(recall_n_array)
        self.fscore_n.set_avg_and_std(fscore_n_array)
        self._is_calculated = True

        scores_array = ScoreArray()
        scores_array.set_arrays(mcc=mcc_array, fpr=fpr_array, recall=recall_array,
                                fscore=fscore_array, bal_acc=bal_acc_array,
                                recall_n=recall_n_array, fscore_n=fscore_n_array,
                                precision=precision_array, precision_n=precision_n_array,
                                tn=tn_array, fp=fp_array, fn=fn_array, tp=tp_array)
        return scores_array

    def get_members(self):
        variables = inspect.getmembers(self, lambda a: not inspect.ismethod(a))
        return [(a[1].name, a[1].label, a[1].average, a[1].std)
                for a in variables if a[0][0] != '_']

    def __repr__(self):
        return '\n'.join('{} {} {} {}'.format(*m) for m in self.get_members())

    def __str__(self):
        return '\n'.join('{:<37} {:<7} {:.2f} {:.2f}'.format(*m) for m in self.get_members())


class ScoreArray(object):
    def __init__(self):
        self.accessions = None
        self.mcc = None
        self.fpr = None
        self.recall = None
        self.fscore = None
        self.bal_acc = None
        self.recall_n = None
        self.fscore_n = None
        self.precision = None
        self.precision_n = None
        self.tp = None
        self.fn = None
        self.tn = None
        self.fp = None
        self._is_calculated = False

    def set_arrays(self, **kwargs):
        for key, val in kwargs.items():
            if hasattr(self, key):
                self.__setattr__(key, val)

    def get_members(self):
        variables = inspect.getmembers(self, lambda a: not inspect.ismethod(a))
        return [a for a in variables if a[0][0] != '_']

    def get_scorearray_asdict(self, key='acc'):
        if key == 'acc' and self.accessions is not None:
            d = dict()
            members = self.get_members()

            accessions = dict(members)['accessions']
            colnames, columns = zip(*list(members))
            for acc, colname, col in zip(accessions, colnames, columns):
                if colname != 'accessions':
                    d[acc] = {colname: list(np.around(col, 3))}

            return d

    def __repr__(self):
        print(type(self))
        return '\n'.join('{} {}'.format(
            attr[0],
            '{}'.format(
                ' '.join('{:.3f}'.format(e) if isinstance(e, float) else e
                         for e in attr[1]))) for attr in self.get_members())

    def __str__(self):
        return '\n'.join('{:<11} {}'.format(
            name,
            '{:<4.2f} {:>4.2f} {:>4.2f} ... {}'.format(
                *elems[:3],
                len(elems)) if isinstance(elems[0], float) else '{} {} {} ... {}'.format(
                *elems[:3],
                len(elems))) for name, elems in self.get_members())
