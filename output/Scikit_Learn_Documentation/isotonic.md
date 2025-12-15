# 1.15.Isotonic regression#

Source URL: https://scikit-learn.org/stable/modules/isotonic.html
Date Scraped: 2025-12-11

---

* [User Guide](../user_guide.html)
* [1. Supervised learning](../supervised_learning.html)
* 1.15. Isotonic regression

The class [`IsotonicRegression`](generated/sklearn.isotonic.IsotonicRegression.html#sklearn.isotonic.IsotonicRegression "sklearn.isotonic.IsotonicRegression") fits a non-decreasing real function to
1-dimensional data. It solves the following problem:

\[\min \sum\_i w\_i (y\_i - \hat{y}\_i)^2\]

subject to \(\hat{y}\_i \le \hat{y}\_j\) whenever \(X\_i \le X\_j\),
where the weights \(w\_i\) are strictly positive, and both `X` and `y` are
arbitrary real quantities.

The `increasing` parameter changes the constraint to
\(\hat{y}\_i \ge \hat{y}\_j\) whenever \(X\_i \le X\_j\). Setting it to
‘auto’ will automatically choose the constraint based on [Spearman’s rank
correlation coefficient](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient).

[`IsotonicRegression`](generated/sklearn.isotonic.IsotonicRegression.html#sklearn.isotonic.IsotonicRegression "sklearn.isotonic.IsotonicRegression") produces a series of predictions
\(\hat{y}\_i\) for the training data which are the closest to the targets
\(y\) in terms of mean squared error. These predictions are interpolated
for predicting to unseen data. The predictions of [`IsotonicRegression`](generated/sklearn.isotonic.IsotonicRegression.html#sklearn.isotonic.IsotonicRegression "sklearn.isotonic.IsotonicRegression")
thus form a function that is piecewise linear:

[![../_images/sphx_glr_plot_isotonic_regression_001.png](../_images/sphx_glr_plot_isotonic_regression_001.png)](../auto_examples/miscellaneous/plot_isotonic_regression.html)

Examples

* [Isotonic Regression](../auto_examples/miscellaneous/plot_isotonic_regression.html#sphx-glr-auto-examples-miscellaneous-plot-isotonic-regression-py)

[previous

1.14. Semi-supervised learning](semi_supervised.html "previous page")
[next

1.16. Probability calibration](calibration.html "next page")

### This Page

* [Show Source](../_sources/modules/isotonic.rst.txt)