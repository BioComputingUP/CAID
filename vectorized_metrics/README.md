## Index
* Introduction
* Dependencies
* Installation
* Assumptions, approximations and workings

## Introduction
`vectorized_cls_metrics` is a library the vectorizes the calculation of binary classification metrics. It revolves 
around `numpy` and builds on the `binary_cls_curve` function from `scikit-learn.metrics`. 

While the latter focuses on
providing a robust and consistent API to obtain classification metrics, `vectorized_cls_metrics` focuses on speed and 
outputs classification metrics for all thresholds in a fraction of time that would normally take with the `sklearn` api 
called in a for-loop.

![performance](results.png)
Time required by `vectorized_cls_metrics` and `sklearn` to perform the same calculation. Reproduce this result with the 
`stats.ipynb` jupyter notebook.  

## Dependencies
* `numpy`
* `pandas`
* `scipy`

## Installation
Git clone this repository

## Assumptions, approximations and workings
See the *vectorizedClsMetrics.pdf* for info on the general approach to metrics calculation.

* while parsing reference states, if a character is found that is not in {0, 1, -} a nan is returned.

* while parsing prediction, if scores are missing, states will be used as scores. if states are missing, they are
generated applying a threshold to scores. The threshold can be passed to parsing function; if it is not passed and 
states are missing, the threshold defaults to 0.5. If both scores and states are missing, the target is excluded from
the analysis.

* prediction scores are rounded to the third decimal figure. This sets the number of possible thresholds to 1000. Rounding
errors are then expected 1/10,000 times: one error each 10,000 labels.

* roc curve uses one additional threshold calculated as max(score) + 1, to ensure the plot starts at point (0, 0)

* precision and recall returned by pr function have one
additional last value (respectively 1, 0) to ensure the plot starts at point (1, 0). An additional threshold of
max(score) + 1 is added to returned thresholds. This is done for two reasons: to ensure consistency with thresholds
returned by roc function; to ensure that returned arrays have same shapes. However it is very important to notice that
this is a mock threshold and it's only there as a placeholder. Differently from the roc function, this threshold is not
used to calculate metrics and metrics associated with this threshold could not be obtained by actually applying the
threshold.

* default threshold is an estimate based on prediction scores. It is calculated as the min score of the positive labels.
Since in metrics calculation thresholds are always applied using a greater-equal strategy, the estimate is functionally
identical to the real value of the threshold.

* targets whose predictions have a different length compared to their references are excluded from the analysis.

* balanced accuracy is the arithmetic mean between tpr and tnr. In case either positive or negative labels are missing
from a target, balanced accuracy is equal to the true rate of the other label.

* fbeta is 0 if denom is 0

* mcc is 0 if denom is 0

* in per-target table missing metrics score are forward filled: nans are replaced with the value from the immediate lower
threshold. Since threshold are applied using a greater-equal strategy, Given two consecutive scores assigned by a
predictor: x and y, where x<y; any value calculated with a threshold between x and y is equal to that calculated using x
as threshold.

* more than 100 bootstrapping samples are hardly beneficial for error estimation. ref:
https://www.misq.org/skin/frontend/default/misq/pdf/appendices/2012/V36I3_Appendices/GoodhueResearchNoteAppendix.pdf

* bootstrapping samples whole dataset with replacements 100 times. Re-sampling is done at the label level (for simplicity)

* confidence intervals are calculated on t distribution at 0.05 alpha.

* shuffled baselines don't have prediction scores since they are just a reshuffling of the reference. Their shuffled
prediction states are however treated as scores for the sake of computation. Default threshold is therefore set to 1.0
since threhsholding is always applied following a greater-equal strategy.