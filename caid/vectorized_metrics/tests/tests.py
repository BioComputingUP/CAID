import pytest
import numpy as np
import pandas as pd
from sklearn import metrics

# local imports
import helpfuncs
import test_results
import vectorized_metrics
import parsers


@pytest.mark.parametrize("rfile", ['ref.test.txt'])
def test_parse_reference(rfile):
    assert parsers.parse_reference(rfile) == test_results.parsed_reftets


def test_parse_predictions():
    assert parsers.parse_prediction("p1.test.txt", {'001', '002', '003'}) == test_results.parsed_p1
    assert parsers.parse_prediction("p2.test.txt", {'001', '002', '003'}) == test_results.parsed_p2


@pytest.mark.parametrize('ytrue,expected', [
    (np.array([0, 0, 0, 0]), np.array([[1, 2, 3, 4], [0, 0, 0, 0], [.6, .3, .2, .1]])),
    (np.array([0, 0, 0, 1]), np.array([[0, 1, 2, 3], [1, 1, 1, 1], [.6, .3, .2, .1]])),
    (np.array([0, 0, 1, 1]), np.array([[0, 1, 2, 2], [1, 1, 1, 2], [.6, .3, .2, .1]])),
    (np.array([0, 1, 1, 1]), np.array([[0, 1, 1, 1], [1, 1, 2, 3], [.6, .3, .2, .1]])),
    (np.array([1, 1, 1, 1]), np.array([[0, 0, 0, 0], [1, 2, 3, 4], [.6, .3, .2, .1]]))])
@pytest.mark.parametrize('yscore', [np.array([.3, .2, .1, .6])])
def test_binary_clf_curve(ytrue, yscore, expected):
    assert np.array_equal(np.array(vectorized_metrics.binary_clf_curve(ytrue, yscore)), expected)
    assert np.array_equal(np.array(vectorized_metrics.binary_clf_curve(ytrue, yscore)),
                          np.array(helpfuncs.binary_clf_curve(ytrue, yscore, 1)))


@pytest.mark.parametrize("ytrue", [np.array([0, 0, 0, 0]),
                                   np.array([0, 0, 0, 1]),
                                   np.array([0, 0, 1, 1]),
                                   np.array([0, 1, 1, 1]),
                                   np.array([1, 1, 1, 1])])
@pytest.mark.parametrize("yscore", [np.array([.3, .2, .1, .6])])
def test_roc(ytrue, yscore):
    assert np.allclose(np.array(vectorized_metrics.roc(*vectorized_metrics.binary_clf_curve(ytrue, yscore))),
                       np.array(metrics.roc_curve(ytrue, yscore, pos_label=1, drop_intermediate=False)),
                       equal_nan=True)


@pytest.mark.parametrize("ytrue", [np.array([0, 0, 0, 0]),
                                   np.array([0, 0, 0, 1]),
                                   np.array([0, 0, 1, 1]),
                                   np.array([0, 1, 1, 1]),
                                   np.array([1, 1, 1, 1])])
@pytest.mark.parametrize("yscore", [np.array([.3, .2, .1, .6])])
def test_pr(ytrue, yscore):
    assert(all(np.allclose(x, y, equal_nan=True) for x, y in zip(
        vectorized_metrics.pr(*vectorized_metrics.binary_clf_curve(ytrue, yscore)),
        metrics.precision_recall_curve(ytrue, yscore, pos_label=1)
    )))


@pytest.mark.parametrize("ytrue", [np.array([0, 0, 0, 0]),
                                   np.array([0, 0, 0, 1]),
                                   np.array([0, 0, 1, 1]),
                                   np.array([0, 1, 1, 1]),
                                   np.array([1, 1, 1, 1])])
@pytest.mark.parametrize("yscore", [np.array([.3, .2, .1, .6])])
def test_confmat(ytrue, yscore):
    cms = []
    for thr in np.sort(np.unique(yscore))[::-1].tolist():
        ypred = np.greater_equal(yscore, thr).astype(np.int)
        cms.append(metrics.confusion_matrix(ytrue, ypred, labels=(0, 1)).ravel())
    assert np.array_equal(np.array(cms).T,
                          vectorized_metrics.confmat(*vectorized_metrics.binary_clf_curve(ytrue, yscore)[:-1]))

def test_get_default_threshold():
    a = np.array([[0, 1, 1, 0], [.3, .5, .6, .4]]).T
    assert vectorized_metrics.calculate_default_threshold(a) == .5
