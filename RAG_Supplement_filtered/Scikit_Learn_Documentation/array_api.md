# 12.1.Array API support (experimental)#

Source URL: https://scikit-learn.org/stable/modules/array_api.html
Date Scraped: 2025-12-11

---

* [User Guide](../user_guide.html)
* [12. Dispatching](../dispatching.html)
* 12.1. Array API support (experimental)

The [Array API](https://data-apis.org/array-api/latest/) specification defines
a standard API for all array manipulation libraries with a NumPy-like API.
Scikit-learn vendors pinned copies of
[array-api-compat](https://github.com/data-apis/array-api-compat)
and [array-api-extra](https://github.com/data-apis/array-api-extra).

Scikit-learn’s support for the array API standard requires the environment variable
`SCIPY_ARRAY_API` to be set to `1` before importing `scipy` and `scikit-learn`:

```
export SCIPY_ARRAY_API=1
```

Please note that this environment variable is intended for temporary use.
For more details, refer to SciPy’s [Array API documentation](https://docs.scipy.org/doc/scipy/dev/api-dev/array_api.html#using-array-api-standard-support).

Some scikit-learn estimators that primarily rely on NumPy (as opposed to using
Cython) to implement the algorithmic logic of their `fit`, `predict` or
`transform` methods can be configured to accept any Array API compatible input
data structures and automatically dispatch operations to the underlying namespace
instead of relying on NumPy.

At this stage, this support is **considered experimental** and must be enabled
explicitly by the `array_api_dispatch` configuration. See below for details.

Note

Currently, only `array-api-strict`, `cupy`, and `PyTorch` are known to work
with scikit-learn’s estimators.

The following video provides an overview of the standard’s design principles
and how it facilitates interoperability between array libraries:

* [Scikit-learn on GPUs with Array API](https://www.youtube.com/watch?v=c_s8tr1AizA)
  by [Thomas Fan](https://github.com/thomasjpfan) at PyData NYC 2023.

## 12.1.1. Enabling array API support[#](#enabling-array-api-support "Link to this heading")

The configuration `array_api_dispatch=True` needs to be set to `True` to enable array
API support. We recommend setting this configuration globally to ensure consistent
behaviour and prevent accidental mixing of array namespaces.
Note that in the examples below, we use a context manager ([`config_context`](generated/sklearn.config_context.html#sklearn.config_context "sklearn.config_context"))
to avoid having to reset it to `False` at the end of every code snippet, so as to
not affect the rest of the documentation.

Scikit-learn accepts [array-like](../glossary.html#term-array-like) inputs for all `metrics`
and some estimators. When `array_api_dispatch=False`, these inputs are converted
into NumPy arrays using [`numpy.asarray`](https://numpy.org/doc/stable/reference/generated/numpy.asarray.html#numpy.asarray "(in NumPy v2.3)") (or [`numpy.array`](https://numpy.org/doc/stable/reference/generated/numpy.array.html#numpy.array "(in NumPy v2.3)")).
While this will successfully convert some array API inputs (e.g., JAX array),
we generally recommend setting `array_api_dispatch=True` when using array API inputs.
This is because NumPy conversion can often fail, e.g., torch tensor allocated on GPU.

## 12.1.2. Example usage[#](#example-usage "Link to this heading")

The example code snippet below demonstrates how to use [CuPy](https://cupy.dev/) to run
[`LinearDiscriminantAnalysis`](generated/sklearn.discriminant_analysis.LinearDiscriminantAnalysis.html#sklearn.discriminant_analysis.LinearDiscriminantAnalysis "sklearn.discriminant_analysis.LinearDiscriminantAnalysis") on a GPU:

```
>>> from sklearn.datasets import make_classification
>>> from sklearn import config_context
>>> from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
>>> import cupy

>>> X_np, y_np = make_classification(random_state=0)
>>> X_cu = cupy.asarray(X_np)
>>> y_cu = cupy.asarray(y_np)
>>> X_cu.device
<CUDA Device 0>

>>> with config_context(array_api_dispatch=True):
...     lda = LinearDiscriminantAnalysis()
...     X_trans = lda.fit_transform(X_cu, y_cu)
>>> X_trans.device
<CUDA Device 0>
```

After the model is trained, fitted attributes that are arrays will also be
from the same Array API namespace as the training data. For example, if CuPy’s
Array API namespace was used for training, then fitted attributes will be on the
GPU. We provide an experimental `_estimator_with_converted_arrays` utility that
transfers an estimator attributes from Array API to an ndarray:

```
>>> from sklearn.utils._array_api import _estimator_with_converted_arrays
>>> cupy_to_ndarray = lambda array : array.get()
>>> lda_np = _estimator_with_converted_arrays(lda, cupy_to_ndarray)
>>> X_trans = lda_np.transform(X_np)
>>> type(X_trans)
<class 'numpy.ndarray'>
```

### 12.1.2.1. PyTorch Support[#](#pytorch-support "Link to this heading")

PyTorch Tensors can also be passed directly:

```
>>> import torch
>>> X_torch = torch.asarray(X_np, device="cuda", dtype=torch.float32)
>>> y_torch = torch.asarray(y_np, device="cuda", dtype=torch.float32)

>>> with config_context(array_api_dispatch=True):
...     lda = LinearDiscriminantAnalysis()
...     X_trans = lda.fit_transform(X_torch, y_torch)
>>> type(X_trans)
<class 'torch.Tensor'>
>>> X_trans.device.type
'cuda'
```

## 12.1.3. Support for `Array API`-compatible inputs[#](#support-for-array-api-compatible-inputs "Link to this heading")

Estimators and other tools in scikit-learn that support Array API compatible inputs.

### 12.1.3.1. Estimators[#](#estimators "Link to this heading")

* [`decomposition.PCA`](generated/sklearn.decomposition.PCA.html#sklearn.decomposition.PCA "sklearn.decomposition.PCA") (with `svd_solver="full"`, `svd_solver="covariance_eigh"`, or
  `svd_solver="randomized"` (`svd_solver="randomized"` only if `power_iteration_normalizer="QR"`))
* [`linear_model.Ridge`](generated/sklearn.linear_model.Ridge.html#sklearn.linear_model.Ridge "sklearn.linear_model.Ridge") (with `solver="svd"`)
* [`linear_model.RidgeCV`](generated/sklearn.linear_model.RidgeCV.html#sklearn.linear_model.RidgeCV "sklearn.linear_model.RidgeCV") (with `solver="svd"`, see [Note on device support for float64](#device-support-for-float64))
* [`linear_model.RidgeClassifier`](generated/sklearn.linear_model.RidgeClassifier.html#sklearn.linear_model.RidgeClassifier "sklearn.linear_model.RidgeClassifier") (with `solver="svd"`)
* [`linear_model.RidgeClassifierCV`](generated/sklearn.linear_model.RidgeClassifierCV.html#sklearn.linear_model.RidgeClassifierCV "sklearn.linear_model.RidgeClassifierCV") (with `solver="svd"`, see [Note on device support for float64](#device-support-for-float64))
* [`discriminant_analysis.LinearDiscriminantAnalysis`](generated/sklearn.discriminant_analysis.LinearDiscriminantAnalysis.html#sklearn.discriminant_analysis.LinearDiscriminantAnalysis "sklearn.discriminant_analysis.LinearDiscriminantAnalysis") (with `solver="svd"`)
* [`naive_bayes.GaussianNB`](generated/sklearn.naive_bayes.GaussianNB.html#sklearn.naive_bayes.GaussianNB "sklearn.naive_bayes.GaussianNB")
* [`preprocessing.Binarizer`](generated/sklearn.preprocessing.Binarizer.html#sklearn.preprocessing.Binarizer "sklearn.preprocessing.Binarizer")
* [`preprocessing.KernelCenterer`](generated/sklearn.preprocessing.KernelCenterer.html#sklearn.preprocessing.KernelCenterer "sklearn.preprocessing.KernelCenterer")
* [`preprocessing.LabelBinarizer`](generated/sklearn.preprocessing.LabelBinarizer.html#sklearn.preprocessing.LabelBinarizer "sklearn.preprocessing.LabelBinarizer") (with `sparse_output=False`)
* [`preprocessing.LabelEncoder`](generated/sklearn.preprocessing.LabelEncoder.html#sklearn.preprocessing.LabelEncoder "sklearn.preprocessing.LabelEncoder")
* [`preprocessing.MaxAbsScaler`](generated/sklearn.preprocessing.MaxAbsScaler.html#sklearn.preprocessing.MaxAbsScaler "sklearn.preprocessing.MaxAbsScaler")
* [`preprocessing.MinMaxScaler`](generated/sklearn.preprocessing.MinMaxScaler.html#sklearn.preprocessing.MinMaxScaler "sklearn.preprocessing.MinMaxScaler")
* [`preprocessing.Normalizer`](generated/sklearn.preprocessing.Normalizer.html#sklearn.preprocessing.Normalizer "sklearn.preprocessing.Normalizer")
* [`preprocessing.PolynomialFeatures`](generated/sklearn.preprocessing.PolynomialFeatures.html#sklearn.preprocessing.PolynomialFeatures "sklearn.preprocessing.PolynomialFeatures")
* [`preprocessing.StandardScaler`](generated/sklearn.preprocessing.StandardScaler.html#sklearn.preprocessing.StandardScaler "sklearn.preprocessing.StandardScaler") (see [Note on device support for float64](#device-support-for-float64))
* [`mixture.GaussianMixture`](generated/sklearn.mixture.GaussianMixture.html#sklearn.mixture.GaussianMixture "sklearn.mixture.GaussianMixture") (with `init_params="random"` or
  `init_params="random_from_data"` and `warm_start=False`)

### 12.1.3.2. Meta-estimators[#](#meta-estimators "Link to this heading")

Meta-estimators that accept Array API inputs conditioned on the fact that the
base estimator also does:

* [`calibration.CalibratedClassifierCV`](generated/sklearn.calibration.CalibratedClassifierCV.html#sklearn.calibration.CalibratedClassifierCV "sklearn.calibration.CalibratedClassifierCV") (with `method="temperature"`)
* [`model_selection.GridSearchCV`](generated/sklearn.model_selection.GridSearchCV.html#sklearn.model_selection.GridSearchCV "sklearn.model_selection.GridSearchCV")
* [`model_selection.RandomizedSearchCV`](generated/sklearn.model_selection.RandomizedSearchCV.html#sklearn.model_selection.RandomizedSearchCV "sklearn.model_selection.RandomizedSearchCV")
* [`model_selection.HalvingGridSearchCV`](generated/sklearn.model_selection.HalvingGridSearchCV.html#sklearn.model_selection.HalvingGridSearchCV "sklearn.model_selection.HalvingGridSearchCV")
* [`model_selection.HalvingRandomSearchCV`](generated/sklearn.model_selection.HalvingRandomSearchCV.html#sklearn.model_selection.HalvingRandomSearchCV "sklearn.model_selection.HalvingRandomSearchCV")

### 12.1.3.3. Metrics[#](#metrics "Link to this heading")

* [`sklearn.metrics.accuracy_score`](generated/sklearn.metrics.accuracy_score.html#sklearn.metrics.accuracy_score "sklearn.metrics.accuracy_score")
* [`sklearn.metrics.balanced_accuracy_score`](generated/sklearn.metrics.balanced_accuracy_score.html#sklearn.metrics.balanced_accuracy_score "sklearn.metrics.balanced_accuracy_score")
* [`sklearn.metrics.brier_score_loss`](generated/sklearn.metrics.brier_score_loss.html#sklearn.metrics.brier_score_loss "sklearn.metrics.brier_score_loss")
* `sklearn.metrics.cluster.calinski_harabasz_score`
* [`sklearn.metrics.cohen_kappa_score`](generated/sklearn.metrics.cohen_kappa_score.html#sklearn.metrics.cohen_kappa_score "sklearn.metrics.cohen_kappa_score")
* [`sklearn.metrics.confusion_matrix`](generated/sklearn.metrics.confusion_matrix.html#sklearn.metrics.confusion_matrix "sklearn.metrics.confusion_matrix")
* [`sklearn.metrics.d2_brier_score`](generated/sklearn.metrics.d2_brier_score.html#sklearn.metrics.d2_brier_score "sklearn.metrics.d2_brier_score")
* [`sklearn.metrics.d2_log_loss_score`](generated/sklearn.metrics.d2_log_loss_score.html#sklearn.metrics.d2_log_loss_score "sklearn.metrics.d2_log_loss_score")
* [`sklearn.metrics.d2_tweedie_score`](generated/sklearn.metrics.d2_tweedie_score.html#sklearn.metrics.d2_tweedie_score "sklearn.metrics.d2_tweedie_score")
* [`sklearn.metrics.det_curve`](generated/sklearn.metrics.det_curve.html#sklearn.metrics.det_curve "sklearn.metrics.det_curve")
* [`sklearn.metrics.explained_variance_score`](generated/sklearn.metrics.explained_variance_score.html#sklearn.metrics.explained_variance_score "sklearn.metrics.explained_variance_score")
* [`sklearn.metrics.f1_score`](generated/sklearn.metrics.f1_score.html#sklearn.metrics.f1_score "sklearn.metrics.f1_score")
* [`sklearn.metrics.fbeta_score`](generated/sklearn.metrics.fbeta_score.html#sklearn.metrics.fbeta_score "sklearn.metrics.fbeta_score")
* [`sklearn.metrics.hamming_loss`](generated/sklearn.metrics.hamming_loss.html#sklearn.metrics.hamming_loss "sklearn.metrics.hamming_loss")
* [`sklearn.metrics.jaccard_score`](generated/sklearn.metrics.jaccard_score.html#sklearn.metrics.jaccard_score "sklearn.metrics.jaccard_score")
* [`sklearn.metrics.log_loss`](generated/sklearn.metrics.log_loss.html#sklearn.metrics.log_loss "sklearn.metrics.log_loss")
* [`sklearn.metrics.max_error`](generated/sklearn.metrics.max_error.html#sklearn.metrics.max_error "sklearn.metrics.max_error")
* [`sklearn.metrics.mean_absolute_error`](generated/sklearn.metrics.mean_absolute_error.html#sklearn.metrics.mean_absolute_error "sklearn.metrics.mean_absolute_error")
* [`sklearn.metrics.mean_absolute_percentage_error`](generated/sklearn.metrics.mean_absolute_percentage_error.html#sklearn.metrics.mean_absolute_percentage_error "sklearn.metrics.mean_absolute_percentage_error")
* [`sklearn.metrics.mean_gamma_deviance`](generated/sklearn.metrics.mean_gamma_deviance.html#sklearn.metrics.mean_gamma_deviance "sklearn.metrics.mean_gamma_deviance")
* [`sklearn.metrics.mean_pinball_loss`](generated/sklearn.metrics.mean_pinball_loss.html#sklearn.metrics.mean_pinball_loss "sklearn.metrics.mean_pinball_loss")
* [`sklearn.metrics.mean_poisson_deviance`](generated/sklearn.metrics.mean_poisson_deviance.html#sklearn.metrics.mean_poisson_deviance "sklearn.metrics.mean_poisson_deviance") (requires [enabling array API support for SciPy](https://docs.scipy.org/doc/scipy/dev/api-dev/array_api.html#using-array-api-standard-support))
* [`sklearn.metrics.mean_squared_error`](generated/sklearn.metrics.mean_squared_error.html#sklearn.metrics.mean_squared_error "sklearn.metrics.mean_squared_error")
* [`sklearn.metrics.mean_squared_log_error`](generated/sklearn.metrics.mean_squared_log_error.html#sklearn.metrics.mean_squared_log_error "sklearn.metrics.mean_squared_log_error")
* [`sklearn.metrics.mean_tweedie_deviance`](generated/sklearn.metrics.mean_tweedie_deviance.html#sklearn.metrics.mean_tweedie_deviance "sklearn.metrics.mean_tweedie_deviance")
* [`sklearn.metrics.median_absolute_error`](generated/sklearn.metrics.median_absolute_error.html#sklearn.metrics.median_absolute_error "sklearn.metrics.median_absolute_error")
* [`sklearn.metrics.multilabel_confusion_matrix`](generated/sklearn.metrics.multilabel_confusion_matrix.html#sklearn.metrics.multilabel_confusion_matrix "sklearn.metrics.multilabel_confusion_matrix")
* [`sklearn.metrics.pairwise.additive_chi2_kernel`](generated/sklearn.metrics.pairwise.additive_chi2_kernel.html#sklearn.metrics.pairwise.additive_chi2_kernel "sklearn.metrics.pairwise.additive_chi2_kernel")
* [`sklearn.metrics.pairwise.chi2_kernel`](generated/sklearn.metrics.pairwise.chi2_kernel.html#sklearn.metrics.pairwise.chi2_kernel "sklearn.metrics.pairwise.chi2_kernel")
* [`sklearn.metrics.pairwise.cosine_similarity`](generated/sklearn.metrics.pairwise.cosine_similarity.html#sklearn.metrics.pairwise.cosine_similarity "sklearn.metrics.pairwise.cosine_similarity")
* [`sklearn.metrics.pairwise.cosine_distances`](generated/sklearn.metrics.pairwise.cosine_distances.html#sklearn.metrics.pairwise.cosine_distances "sklearn.metrics.pairwise.cosine_distances")
* `sklearn.metrics.pairwise.pairwise_distances` (only supports “cosine”, “euclidean”, “manhattan” and “l2” metrics)
* [`sklearn.metrics.pairwise.euclidean_distances`](generated/sklearn.metrics.pairwise.euclidean_distances.html#sklearn.metrics.pairwise.euclidean_distances "sklearn.metrics.pairwise.euclidean_distances") (see [Note on device support for float64](#device-support-for-float64))
* [`sklearn.metrics.pairwise.laplacian_kernel`](generated/sklearn.metrics.pairwise.laplacian_kernel.html#sklearn.metrics.pairwise.laplacian_kernel "sklearn.metrics.pairwise.laplacian_kernel")
* [`sklearn.metrics.pairwise.linear_kernel`](generated/sklearn.metrics.pairwise.linear_kernel.html#sklearn.metrics.pairwise.linear_kernel "sklearn.metrics.pairwise.linear_kernel")
* [`sklearn.metrics.pairwise.manhattan_distances`](generated/sklearn.metrics.pairwise.manhattan_distances.html#sklearn.metrics.pairwise.manhattan_distances "sklearn.metrics.pairwise.manhattan_distances")
* [`sklearn.metrics.pairwise.paired_cosine_distances`](generated/sklearn.metrics.pairwise.paired_cosine_distances.html#sklearn.metrics.pairwise.paired_cosine_distances "sklearn.metrics.pairwise.paired_cosine_distances")
* [`sklearn.metrics.pairwise.paired_euclidean_distances`](generated/sklearn.metrics.pairwise.paired_euclidean_distances.html#sklearn.metrics.pairwise.paired_euclidean_distances "sklearn.metrics.pairwise.paired_euclidean_distances")
* [`sklearn.metrics.pairwise.pairwise_kernels`](generated/sklearn.metrics.pairwise.pairwise_kernels.html#sklearn.metrics.pairwise.pairwise_kernels "sklearn.metrics.pairwise.pairwise_kernels")
* [`sklearn.metrics.pairwise.polynomial_kernel`](generated/sklearn.metrics.pairwise.polynomial_kernel.html#sklearn.metrics.pairwise.polynomial_kernel "sklearn.metrics.pairwise.polynomial_kernel")
* [`sklearn.metrics.pairwise.rbf_kernel`](generated/sklearn.metrics.pairwise.rbf_kernel.html#sklearn.metrics.pairwise.rbf_kernel "sklearn.metrics.pairwise.rbf_kernel") (see [Note on device support for float64](#device-support-for-float64))
* [`sklearn.metrics.pairwise.sigmoid_kernel`](generated/sklearn.metrics.pairwise.sigmoid_kernel.html#sklearn.metrics.pairwise.sigmoid_kernel "sklearn.metrics.pairwise.sigmoid_kernel")
* [`sklearn.metrics.precision_score`](generated/sklearn.metrics.precision_score.html#sklearn.metrics.precision_score "sklearn.metrics.precision_score")
* [`sklearn.metrics.precision_recall_curve`](generated/sklearn.metrics.precision_recall_curve.html#sklearn.metrics.precision_recall_curve "sklearn.metrics.precision_recall_curve")
* [`sklearn.metrics.precision_recall_fscore_support`](generated/sklearn.metrics.precision_recall_fscore_support.html#sklearn.metrics.precision_recall_fscore_support "sklearn.metrics.precision_recall_fscore_support")
* [`sklearn.metrics.r2_score`](generated/sklearn.metrics.r2_score.html#sklearn.metrics.r2_score "sklearn.metrics.r2_score")
* [`sklearn.metrics.recall_score`](generated/sklearn.metrics.recall_score.html#sklearn.metrics.recall_score "sklearn.metrics.recall_score")
* [`sklearn.metrics.roc_curve`](generated/sklearn.metrics.roc_curve.html#sklearn.metrics.roc_curve "sklearn.metrics.roc_curve")
* [`sklearn.metrics.root_mean_squared_error`](generated/sklearn.metrics.root_mean_squared_error.html#sklearn.metrics.root_mean_squared_error "sklearn.metrics.root_mean_squared_error")
* [`sklearn.metrics.root_mean_squared_log_error`](generated/sklearn.metrics.root_mean_squared_log_error.html#sklearn.metrics.root_mean_squared_log_error "sklearn.metrics.root_mean_squared_log_error")
* [`sklearn.metrics.zero_one_loss`](generated/sklearn.metrics.zero_one_loss.html#sklearn.metrics.zero_one_loss "sklearn.metrics.zero_one_loss")

### 12.1.3.4. Tools[#](#tools "Link to this heading")

* [`preprocessing.label_binarize`](generated/sklearn.preprocessing.label_binarize.html#sklearn.preprocessing.label_binarize "sklearn.preprocessing.label_binarize") (with `sparse_output=False`)
* [`model_selection.cross_val_predict`](generated/sklearn.model_selection.cross_val_predict.html#sklearn.model_selection.cross_val_predict "sklearn.model_selection.cross_val_predict")
* [`model_selection.train_test_split`](generated/sklearn.model_selection.train_test_split.html#sklearn.model_selection.train_test_split "sklearn.model_selection.train_test_split")
* [`utils.check_consistent_length`](generated/sklearn.utils.check_consistent_length.html#sklearn.utils.check_consistent_length "sklearn.utils.check_consistent_length")

Coverage is expected to grow over time. Please follow the dedicated [meta-issue on GitHub](https://github.com/scikit-learn/scikit-learn/issues/22352) to track progress.

## 12.1.4. Input and output array type handling[#](#input-and-output-array-type-handling "Link to this heading")

Estimators and scoring functions are able to accept input arrays
from different array libraries and/or devices. When a mixed set of input arrays is
passed, scikit-learn converts arrays as needed to make them all consistent.

For estimators, the rule is **“everything follows** `X` **“** - mixed array inputs are
converted so that they all match the array library and device of `X`.
For scoring functions the rule is **“everything follows** `y_pred` **“** - mixed array
inputs are converted so that they all match the array library and device of `y_pred`.

When a function or method has been called with array API compatible inputs, the
convention is to return arrays from the same array library and on the same
device as the input data.

### 12.1.4.1. Estimators[#](#id2 "Link to this heading")

When an estimator is fitted with an array API compatible `X`, all other
array inputs, including constructor arguments, (e.g., `y`, `sample_weight`)
will be converted to match the array library and device of `X`, if they do not already.
This behaviour enables switching from processing on the CPU to processing
on the GPU at any point within a pipeline.

This allows estimators to accept mixed input types, enabling `X` to be moved
to a different device within a pipeline, without explicitly moving `y`.
Note that scikit-learn pipelines do not allow transformation of `y` (to avoid
[leakage](../common_pitfalls.html#data-leakage)).

Take for example a pipeline where `X` and `y` both start on CPU, and go through
the following three steps:

* [`TargetEncoder`](generated/sklearn.preprocessing.TargetEncoder.html#sklearn.preprocessing.TargetEncoder "sklearn.preprocessing.TargetEncoder"), which will transform categorial
  `X` but also requires `y`, meaning both `X` and `y` need to be on CPU.
* [`FunctionTransformer(func=partial(torch.asarray, device="cuda"))`](generated/sklearn.preprocessing.FunctionTransformer.html#sklearn.preprocessing.FunctionTransformer "sklearn.preprocessing.FunctionTransformer"),
  which moves `X` to GPU, to improve performance in the next step.
* [`Ridge`](generated/sklearn.linear_model.Ridge.html#sklearn.linear_model.Ridge "sklearn.linear_model.Ridge"), whose performance can be improved when
  passed arrays on a GPU, as they can handle large matrix operations very efficiently.

`X` initially contains categorical string data (thus needs to be on CPU), which is
target encoded to numerical values in [`TargetEncoder`](generated/sklearn.preprocessing.TargetEncoder.html#sklearn.preprocessing.TargetEncoder "sklearn.preprocessing.TargetEncoder").
`X` is then explicitly moved to GPU to improve the performance of
[`Ridge`](generated/sklearn.linear_model.Ridge.html#sklearn.linear_model.Ridge "sklearn.linear_model.Ridge"). `y` cannot be transformed by the pipeline
(recall scikit-learn pipelines do not allow transformation of `y`) but as
[`Ridge`](generated/sklearn.linear_model.Ridge.html#sklearn.linear_model.Ridge "sklearn.linear_model.Ridge") is able to accept mixed input types,
this is not a problem and the pipeline is able to be run.

The fitted attributes of an estimator fitted with an array API compatible `X`, will
be arrays from the same library as the input and stored on the same device.
The `predict` and `transform` method subsequently expect
inputs from the same array library and device as the data passed to the `fit`
method.

### 12.1.4.2. Scoring functions[#](#scoring-functions "Link to this heading")

When an array API compatible `y_pred` is passed to a scoring function,
all other array inputs (e.g., `y_true`, `sample_weight`) will be converted
to match the array library and device of `y_pred`, if they do not already.
This allows scoring functions to accept mixed input types, enabling them to be
used within a [meta-estimator](../glossary.html#term-meta-estimator) (or function that accepts estimators), with a
pipeline that moves input arrays between devices (e.g., CPU to GPU).

For example, to be able to use the pipeline described above within e.g.,
[`cross_validate`](generated/sklearn.model_selection.cross_validate.html#sklearn.model_selection.cross_validate "sklearn.model_selection.cross_validate") or
[`GridSearchCV`](generated/sklearn.model_selection.GridSearchCV.html#sklearn.model_selection.GridSearchCV "sklearn.model_selection.GridSearchCV"), the scoring function internally
called needs to be able to accept mixed input types.

The output type of scoring functions depends on the number of output values.
When a scoring function returns a scalar value, it will return a Python
scalar (typically a `float` instance) instead of an array scalar value.
For scoring functions that support [multiclass](../glossary.html#term-multiclass) or [multioutput](../glossary.html#term-multioutput),
an array from the same array library and device as `y_pred` will be returned when
multiple values need to be output.

## 12.1.5. Common estimator checks[#](#common-estimator-checks "Link to this heading")

Add the `array_api_support` tag to an estimator’s set of tags to indicate that
it supports the array API. This will enable dedicated checks as part of the
common tests to verify that the estimators’ results are the same when using
vanilla NumPy and array API inputs.

To run these checks you need to install
[array-api-strict](https://data-apis.org/array-api-strict/) in your
test environment. This allows you to run checks without having a
GPU. To run the full set of checks you also need to install
[PyTorch](https://pytorch.org/), [CuPy](https://cupy.dev/) and have
a GPU. Checks that can not be executed or have missing dependencies will be
automatically skipped. Therefore it’s important to run the tests with the
`-v` flag to see which checks are skipped:

```
pip install array-api-strict  # and other libraries as needed
pytest -k "array_api" -v
```

Running the scikit-learn tests against `array-api-strict` should help reveal
most code problems related to handling multiple device inputs via the use of
simulated non-CPU devices. This allows for fast iterative development and debugging of
array API related code.

However, to ensure full handling of PyTorch or CuPy inputs allocated on actual GPU
devices, it is necessary to run the tests against those libraries and hardware.
This can either be achieved by using
[Google Colab](https://gist.github.com/EdAbati/ff3bdc06bafeb92452b3740686cc8d7c)
or leveraging our CI infrastructure on pull requests (manually triggered by maintainers
for cost reasons).

### 12.1.5.1. Note on MPS device support[#](#note-on-mps-device-support "Link to this heading")

On macOS, PyTorch can use the Metal Performance Shaders (MPS) to access
hardware accelerators (e.g. the internal GPU component of the M1 or M2 chips).
However, the MPS device support for PyTorch is incomplete at the time of
writing. See the following github issue for more details:

* [pytorch/pytorch#77764](https://github.com/pytorch/pytorch/issues/77764)

To enable the MPS support in PyTorch, set the environment variable
`PYTORCH_ENABLE_MPS_FALLBACK=1` before running the tests:

```
PYTORCH_ENABLE_MPS_FALLBACK=1 pytest -k "array_api" -v
```

At the time of writing all scikit-learn tests should pass, however, the
computational speed is not necessarily better than with the CPU device.

### 12.1.5.2. Note on device support for `float64`[#](#note-on-device-support-for-float64 "Link to this heading")

Certain operations within scikit-learn will automatically perform operations
on floating-point values with `float64` precision to prevent overflows and ensure
correctness (e.g., [`metrics.pairwise.euclidean_distances`](generated/sklearn.metrics.pairwise.euclidean_distances.html#sklearn.metrics.pairwise.euclidean_distances "sklearn.metrics.pairwise.euclidean_distances"),
[`preprocessing.StandardScaler`](generated/sklearn.preprocessing.StandardScaler.html#sklearn.preprocessing.StandardScaler "sklearn.preprocessing.StandardScaler")). However,
certain combinations of array namespaces and devices, such as `PyTorch on MPS`
(see [Note on MPS device support](#mps-support)) do not support the `float64` data type. In these cases,
scikit-learn will revert to using the `float32` data type instead. This can result in
different behavior (typically numerically unstable results) compared to not using array
API dispatching or using a device with `float64` support.

[previous

12. Dispatching](../dispatching.html "previous page")
[next

13. Choosing the right estimator](../machine_learning_map.html "next page")

On this page

* [12.1.1. Enabling array API support](#enabling-array-api-support)
* [12.1.2. Example usage](#example-usage)
  + [12.1.2.1. PyTorch Support](#pytorch-support)
* [12.1.3. Support for `Array API`-compatible inputs](#support-for-array-api-compatible-inputs)
  + [12.1.3.1. Estimators](#estimators)
  + [12.1.3.2. Meta-estimators](#meta-estimators)
  + [12.1.3.3. Metrics](#metrics)
  + [12.1.3.4. Tools](#tools)
* [12.1.4. Input and output array type handling](#input-and-output-array-type-handling)
  + [12.1.4.1. Estimators](#id2)
  + [12.1.4.2. Scoring functions](#scoring-functions)
* [12.1.5. Common estimator checks](#common-estimator-checks)
  + [12.1.5.1. Note on MPS device support](#note-on-mps-device-support)
  + [12.1.5.2. Note on device support for `float64`](#note-on-device-support-for-float64)

### This Page

* [Show Source](../_sources/modules/array_api.rst.txt)