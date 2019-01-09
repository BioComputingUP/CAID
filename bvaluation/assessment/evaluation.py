import logging
import numpy as np
from sklearn import metrics
from collections import OrderedDict
# relative imports
from bvaluation.assessment.scores import Scores, Curves, AverageInstanceScores, ScoreArray


class EvaluationArray(object):
    def __init__(self, eval_list=None):
        self.body = list(eval_list) if eval_list is not None else list()
        self.avg = None
        self.std = None

    def __add__(self, eval_obj):
        if isinstance(eval_obj, Evaluation):
            return self.body.append(eval_obj)
        else:
            raise TypeError('unsupported operand type(s) for +: \'EvaluationArray\' and \'%s\'',
                            type(eval_obj))


class Evaluation(object):
    """
    Evaluation of a prediction against a reference
    """
    def __init__(self):
        self.overall_scores = Scores()
        self.curves = Curves()
        self.avg_instance_scores = AverageInstanceScores()
        self.per_instance_scores = ScoreArray()

    def calc_overall_scores(self, reference, prediction):
        if len(reference) == len(prediction):
            confusion_matrix = metrics.confusion_matrix(reference, prediction,
                                                    labels=(0, 1)).astype(np.float)
            self.overall_scores.states_based_metrics(confusion_matrix=confusion_matrix)

    def calc_curves(self, reference, scores):
        self.curves.calc_curves(reference, scores)

    def calc_average_instance_scores(self, reference, prediction, save_perinstance_scores=True,
                                     accessions=None):
        score_array = self.avg_instance_scores.from_iterable(reference, prediction)

        if save_perinstance_scores is True:
            if accessions:
                score_array.accessions = accessions
            self.per_instance_scores = score_array
            self.per_instance_scores._is_calculated = True

    def get_scores_aslist(self):
        vals = list()
        lbls = list()

        if self.overall_scores._is_calculated is True:
            _, labels, scores = zip(*self.overall_scores.get_members())
            lbls.extend(labels)
            vals.extend(scores)

        if self.avg_instance_scores._is_calculated is True:
            _, labels, *scores = zip(*self.avg_instance_scores.get_members())
            lbls.extend(sum((('{}_avg'.format(l), '{}_std'.format(l)) for l in labels), ()))
            vals.extend(sum(scores, ()))

        if self.curves._is_calculated is True:
            _, labels, scores, *_ = zip(*self.curves.get_members())
            lbls.extend(labels)
            vals.extend(scores)

        return lbls, vals

    def get_scores_asdict(self):
        return OrderedDict(self.get_scores_aszip())

    def get_scores_aszip(self, lazy=True):
        z = zip(*self.get_scores_aslist())
        return list(z) if lazy is False else z

    def __repr__(self):
        rpr = '\n'.join('{} {}'.format(label, score)
                        for label, score in zip(*self.get_scores_aslist()))
        rpr += '\n'
        rpr += self.get_curves_repr()

        return rpr

    def __str__(self):
        return '\n'.join('{:<10} {:.3f}'.format(l, s) for l, s in zip(*self.get_scores_aslist()))

    def get_scores(self, dstr=dict):
        if dstr is dict:
            return self.get_scores_asdict()
        elif dstr is zip:
            return self.get_scores_aszip()
        elif dstr is list:
            return self.get_scores_aslist()
        else:
            logging.error('Data structure not implemented')

    def get_curves_aslist(self):
        crvs = list()
        sorting = ['ROC', 'PRC']

        if self.curves._is_calculated is True:
            for curve in sorted(self.curves.get_members(), key=lambda c: sorting.index(c[4])):
                _, _, auc, _, label, x, y, thr = curve
                crvs.append((label, auc, list(zip(np.around(x, 3).tolist(),
                                                  np.around(y, 3).tolist(),
                                                  np.around(thr, 3).tolist()))))
        return crvs

    def get_curves_asdict(self):
        return {label: (auc, points) for label, auc, points in self.get_curves_aslist()}

    def get_curves(self, dstr=list):
        if dstr is list:
            return self.get_scores_aslist()
        if dstr is dict:
            return self.get_curves_asdict()
        pass

    def get_curves_repr(self):
        return '\n'.join('{} {:.3f} {}'.format(
            label, auc, ' '.join('{},{},{}'.format(x, y, t) for (x, y, t) in points))
                  for label, auc, points in self.get_curves_aslist())

    def get_perinstance_scores(self, insert=None):
        if insert is None:
            scores = [[e[0]] + list(e[1]) for e in self.per_instance_scores.get_members()]
        else:
            scores = [[insert, e[0]] + list(e[1]) for e in self.per_instance_scores.get_members()]

        return scores
