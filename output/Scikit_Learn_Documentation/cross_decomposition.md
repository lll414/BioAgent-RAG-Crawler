# 1.8.Cross decomposition#

Source URL: https://scikit-learn.org/stable/modules/cross_decomposition.html
Date Scraped: 2025-12-11

---

* [User Guide](../user_guide.html)
* [1. Supervised learning](../supervised_learning.html)
* 1.8. Cross decomposition

The cross decomposition module contains **supervised** estimators for
dimensionality reduction and regression, belonging to the “Partial Least
Squares” family.

[![../_images/sphx_glr_plot_compare_cross_decomposition_001.png](../_images/sphx_glr_plot_compare_cross_decomposition_001.png)](../auto_examples/cross_decomposition/plot_compare_cross_decomposition.html)

Cross decomposition algorithms find the fundamental relations between two
matrices (X and Y). They are latent variable approaches to modeling the
covariance structures in these two spaces. They will try to find the
multidimensional direction in the X space that explains the maximum
multidimensional variance direction in the Y space. In other words, PLS
projects both `X` and `Y` into a lower-dimensional subspace such that the
covariance between `transformed(X)` and `transformed(Y)` is maximal.

PLS draws similarities with [Principal Component Regression](https://en.wikipedia.org/wiki/Principal_component_regression) (PCR), where
the samples are first projected into a lower-dimensional subspace, and the
targets `y` are predicted using `transformed(X)`. One issue with PCR is that
the dimensionality reduction is unsupervised, and may lose some important
variables: PCR would keep the features with the most variance, but it’s
possible that features with small variances are relevant for predicting
the target. In a way, PLS allows for the same kind of dimensionality
reduction, but by taking into account the targets `y`. An illustration of
this fact is given in the following example:
\* [Principal Component Regression vs Partial Least Squares Regression](../auto_examples/cross_decomposition/plot_pcr_vs_pls.html#sphx-glr-auto-examples-cross-decomposition-plot-pcr-vs-pls-py).

Apart from CCA, the PLS estimators are particularly suited when the matrix of
predictors has more variables than observations, and when there is
multicollinearity among the features. By contrast, standard linear regression
would fail in these cases unless it is regularized.

Classes included in this module are [`PLSRegression`](generated/sklearn.cross_decomposition.PLSRegression.html#sklearn.cross_decomposition.PLSRegression "sklearn.cross_decomposition.PLSRegression"),
[`PLSCanonical`](generated/sklearn.cross_decomposition.PLSCanonical.html#sklearn.cross_decomposition.PLSCanonical "sklearn.cross_decomposition.PLSCanonical"), [`CCA`](generated/sklearn.cross_decomposition.CCA.html#sklearn.cross_decomposition.CCA "sklearn.cross_decomposition.CCA") and [`PLSSVD`](generated/sklearn.cross_decomposition.PLSSVD.html#sklearn.cross_decomposition.PLSSVD "sklearn.cross_decomposition.PLSSVD")

## 1.8.1. PLSCanonical[#](#plscanonical "Link to this heading")

We here describe the algorithm used in [`PLSCanonical`](generated/sklearn.cross_decomposition.PLSCanonical.html#sklearn.cross_decomposition.PLSCanonical "sklearn.cross_decomposition.PLSCanonical"). The other
estimators use variants of this algorithm, and are detailed below.
We recommend section [[1]](#id6) for more details and comparisons between these
algorithms. In [[1]](#id6), [`PLSCanonical`](generated/sklearn.cross_decomposition.PLSCanonical.html#sklearn.cross_decomposition.PLSCanonical "sklearn.cross_decomposition.PLSCanonical") corresponds to “PLSW2A”.

Given two centered matrices \(X \in \mathbb{R}^{n \times d}\) and
\(Y \in \mathbb{R}^{n \times t}\), and a number of components \(K\),
[`PLSCanonical`](generated/sklearn.cross_decomposition.PLSCanonical.html#sklearn.cross_decomposition.PLSCanonical "sklearn.cross_decomposition.PLSCanonical") proceeds as follows:

Set \(X\_1\) to \(X\) and \(Y\_1\) to \(Y\). Then, for each
\(k \in [1, K]\):

* a) compute \(u\_k \in \mathbb{R}^d\) and \(v\_k \in \mathbb{R}^t\),
  the first left and right singular vectors of the cross-covariance matrix
  \(C = X\_k^T Y\_k\).
  \(u\_k\) and \(v\_k\) are called the *weights*.
  By definition, \(u\_k\) and \(v\_k\) are
  chosen so that they maximize the covariance between the projected
  \(X\_k\) and the projected target, that is \(\text{Cov}(X\_k u\_k,
  Y\_k v\_k)\).
* b) Project \(X\_k\) and \(Y\_k\) on the singular vectors to obtain
  *scores*: \(\xi\_k = X\_k u\_k\) and \(\omega\_k = Y\_k v\_k\)
* c) Regress \(X\_k\) on \(\xi\_k\), i.e. find a vector \(\gamma\_k
  \in \mathbb{R}^d\) such that the rank-1 matrix \(\xi\_k \gamma\_k^T\)
  is as close as possible to \(X\_k\). Do the same on \(Y\_k\) with
  \(\omega\_k\) to obtain \(\delta\_k\). The vectors
  \(\gamma\_k\) and \(\delta\_k\) are called the *loadings*.
* d) *deflate* \(X\_k\) and \(Y\_k\), i.e. subtract the rank-1
  approximations: \(X\_{k+1} = X\_k - \xi\_k \gamma\_k^T\), and
  \(Y\_{k + 1} = Y\_k - \omega\_k \delta\_k^T\).

At the end, we have approximated \(X\) as a sum of rank-1 matrices:
\(X = \Xi \Gamma^T\) where \(\Xi \in \mathbb{R}^{n \times K}\)
contains the scores in its columns, and \(\Gamma^T \in \mathbb{R}^{K
\times d}\) contains the loadings in its rows. Similarly for \(Y\), we
have \(Y = \Omega \Delta^T\).

Note that the scores matrices \(\Xi\) and \(\Omega\) correspond to
the projections of the training data \(X\) and \(Y\), respectively.

Step *a)* may be performed in two ways: either by computing the whole SVD of
\(C\) and only retaining the singular vectors with the biggest singular
values, or by directly computing the singular vectors using the power method (cf section 11.3 in [[1]](#id6)),
which corresponds to the `'nipals'` option of the `algorithm` parameter.

Transforming data[#](#transforming-data "Link to this dropdown")

To transform \(X\) into \(\bar{X}\), we need to find a projection
matrix \(P\) such that \(\bar{X} = XP\). We know that for the
training data, \(\Xi = XP\), and \(X = \Xi \Gamma^T\). Setting
\(P = U(\Gamma^T U)^{-1}\) where \(U\) is the matrix with the
\(u\_k\) in the columns, we have \(XP = X U(\Gamma^T U)^{-1} = \Xi
(\Gamma^T U) (\Gamma^T U)^{-1} = \Xi\) as desired. The rotation matrix
\(P\) can be accessed from the `x_rotations_` attribute.

Similarly, \(Y\) can be transformed using the rotation matrix
\(V(\Delta^T V)^{-1}\), accessed via the `y_rotations_` attribute.

Predicting the targets `Y`[#](#predicting-the-targets-y "Link to this dropdown")

To predict the targets of some data \(X\), we are looking for a
coefficient matrix \(\beta \in R^{d \times t}\) such that \(Y =
X\beta\).

The idea is to try to predict the transformed targets \(\Omega\) as a
function of the transformed samples \(\Xi\), by computing \(\alpha
\in \mathbb{R}\) such that \(\Omega = \alpha \Xi\).

Then, we have \(Y = \Omega \Delta^T = \alpha \Xi \Delta^T\), and since
\(\Xi\) is the transformed training data we have that \(Y = X \alpha
P \Delta^T\), and as a result the coefficient matrix \(\beta = \alpha P
\Delta^T\).

\(\beta\) can be accessed through the `coef_` attribute.

## 1.8.2. PLSSVD[#](#plssvd "Link to this heading")

[`PLSSVD`](generated/sklearn.cross_decomposition.PLSSVD.html#sklearn.cross_decomposition.PLSSVD "sklearn.cross_decomposition.PLSSVD") is a simplified version of [`PLSCanonical`](generated/sklearn.cross_decomposition.PLSCanonical.html#sklearn.cross_decomposition.PLSCanonical "sklearn.cross_decomposition.PLSCanonical")
described earlier: instead of iteratively deflating the matrices \(X\_k\)
and \(Y\_k\), [`PLSSVD`](generated/sklearn.cross_decomposition.PLSSVD.html#sklearn.cross_decomposition.PLSSVD "sklearn.cross_decomposition.PLSSVD") computes the SVD of \(C = X^TY\)
only *once*, and stores the `n_components` singular vectors corresponding to
the biggest singular values in the matrices `U` and `V`, corresponding to the
`x_weights_` and `y_weights_` attributes. Here, the transformed data is
simply `transformed(X) = XU` and `transformed(Y) = YV`.

If `n_components == 1`, [`PLSSVD`](generated/sklearn.cross_decomposition.PLSSVD.html#sklearn.cross_decomposition.PLSSVD "sklearn.cross_decomposition.PLSSVD") and [`PLSCanonical`](generated/sklearn.cross_decomposition.PLSCanonical.html#sklearn.cross_decomposition.PLSCanonical "sklearn.cross_decomposition.PLSCanonical") are
strictly equivalent.

## 1.8.3. PLSRegression[#](#plsregression "Link to this heading")

The [`PLSRegression`](generated/sklearn.cross_decomposition.PLSRegression.html#sklearn.cross_decomposition.PLSRegression "sklearn.cross_decomposition.PLSRegression") estimator is similar to
[`PLSCanonical`](generated/sklearn.cross_decomposition.PLSCanonical.html#sklearn.cross_decomposition.PLSCanonical "sklearn.cross_decomposition.PLSCanonical") with `algorithm='nipals'`, with 2 significant
differences:

* at step a) in the power method to compute \(u\_k\) and \(v\_k\),
  \(v\_k\) is never normalized.
* at step c), the targets \(Y\_k\) are approximated using the projection
  of \(X\_k\) (i.e. \(\xi\_k\)) instead of the projection of
  \(Y\_k\) (i.e. \(\omega\_k\)). In other words, the loadings
  computation is different. As a result, the deflation in step d) will also
  be affected.

These two modifications affect the output of `predict` and `transform`,
which are not the same as for [`PLSCanonical`](generated/sklearn.cross_decomposition.PLSCanonical.html#sklearn.cross_decomposition.PLSCanonical "sklearn.cross_decomposition.PLSCanonical"). Also, while the number
of components is limited by `min(n_samples, n_features, n_targets)` in
[`PLSCanonical`](generated/sklearn.cross_decomposition.PLSCanonical.html#sklearn.cross_decomposition.PLSCanonical "sklearn.cross_decomposition.PLSCanonical"), here the limit is the rank of \(X^TX\), i.e.
`min(n_samples, n_features)`.

[`PLSRegression`](generated/sklearn.cross_decomposition.PLSRegression.html#sklearn.cross_decomposition.PLSRegression "sklearn.cross_decomposition.PLSRegression") is also known as PLS1 (single targets) and PLS2
(multiple targets). Much like [`Lasso`](generated/sklearn.linear_model.Lasso.html#sklearn.linear_model.Lasso "sklearn.linear_model.Lasso"),
[`PLSRegression`](generated/sklearn.cross_decomposition.PLSRegression.html#sklearn.cross_decomposition.PLSRegression "sklearn.cross_decomposition.PLSRegression") is a form of regularized linear regression where the
number of components controls the strength of the regularization.

## 1.8.4. Canonical Correlation Analysis[#](#canonical-correlation-analysis "Link to this heading")

Canonical Correlation Analysis was developed prior and independently to PLS.
But it turns out that [`CCA`](generated/sklearn.cross_decomposition.CCA.html#sklearn.cross_decomposition.CCA "sklearn.cross_decomposition.CCA") is a special case of PLS, and corresponds
to PLS in “Mode B” in the literature.

[`CCA`](generated/sklearn.cross_decomposition.CCA.html#sklearn.cross_decomposition.CCA "sklearn.cross_decomposition.CCA") differs from [`PLSCanonical`](generated/sklearn.cross_decomposition.PLSCanonical.html#sklearn.cross_decomposition.PLSCanonical "sklearn.cross_decomposition.PLSCanonical") in the way the weights
\(u\_k\) and \(v\_k\) are computed in the power method of step a).
Details can be found in section 10 of [[1]](#id6).

Since [`CCA`](generated/sklearn.cross_decomposition.CCA.html#sklearn.cross_decomposition.CCA "sklearn.cross_decomposition.CCA") involves the inversion of \(X\_k^TX\_k\) and
\(Y\_k^TY\_k\), this estimator can be unstable if the number of features or
targets is greater than the number of samples.

References

[1]
([1](#id2),[2](#id3),[3](#id4),[4](#id5))

[A survey of Partial Least Squares (PLS) methods, with emphasis on the two-block
case](https://stat.uw.edu/sites/default/files/files/reports/2000/tr371.pdf),
JA Wegelin

Examples

* [Compare cross decomposition methods](../auto_examples/cross_decomposition/plot_compare_cross_decomposition.html#sphx-glr-auto-examples-cross-decomposition-plot-compare-cross-decomposition-py)
* [Principal Component Regression vs Partial Least Squares Regression](../auto_examples/cross_decomposition/plot_pcr_vs_pls.html#sphx-glr-auto-examples-cross-decomposition-plot-pcr-vs-pls-py)

[previous

1.7. Gaussian Processes](gaussian_process.html "previous page")
[next

1.9. Naive Bayes](naive_bayes.html "next page")

On this page

* [1.8.1. PLSCanonical](#plscanonical)
* [1.8.2. PLSSVD](#plssvd)
* [1.8.3. PLSRegression](#plsregression)
* [1.8.4. Canonical Correlation Analysis](#canonical-correlation-analysis)

### This Page

* [Show Source](../_sources/modules/cross_decomposition.rst.txt)