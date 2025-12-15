# 13.Choosing the right estimator#

Source URL: https://scikit-learn.org/stable/machine_learning_map.html
Date Scraped: 2025-12-11

---

* [User Guide](user_guide.html)
* 13. Choosing the right estimator

Often the hardest part of solving a machine learning problem can be finding the right
estimator for the job. Different estimators are better suited for different types of
data and different problems.

The flowchart below is designed to give users a bit of a rough guide on how to approach
problems with regard to which estimators to try on your data. Click on any estimator in
the chart below to see its documentation. The **Try next** orange arrows are to be read as
“if this estimator does not achieve the desired outcome, then follow the arrow and try
the next one”. Use scroll wheel to zoom in and out, and click and drag to pan around.
You can also download the chart: [`ml_map.svg`](_downloads/b82bf6cd7438a351f19fac60fbc0d927/ml_map.svg).




xml version="1.0" encoding="UTF-8"?




**START**

START

>50

samples

>50...

get

more

data

get...

NO

NO

predicting a

category

predicting...

YES

YES

do you have

labeled

data

do you hav...

YES

YES

predicting a

quantity

predicting...

NO

NO

just

looking

just...

NO

NO

predicting

structure

predicting...

NO

NO

tough

luck

tough...

<100K

samples

<100K...

YES

YESSGD

Classifier

SGD...

NO

NOLinear

SVC

Linear...

YES

YES

text

data

text...Kernel

Approximation

Kernel...KNeighbors

Classifier

KNeighbors...

NO

NOSVC

SVCEnsemble

Classifiers

Ensemble...Naive

Bayes

Naive...

YES

YES

classification

classification

number of

categories

known

number of...

NO

NO

<10K

samples

<10K...

<10K

samples

<10K...

NO

NO

NO

NO

YES

YESMeanShift

MeanShiftVBGMM

VBGMM

YES

YESMiniBatch

KMeans

MiniBatch...

NO

NO

clustering

clusteringKMeans

KMeans

YES

YESSpectral

Clustering

Spectral...GMM

GMM

<100K

samples

<100K...

YES

YES

few features

should be

important

few features...

YES

YESSGD

Regressor

SGD...

NO

NOLasso

LassoElasticNet

ElasticNet

YES

YESRidgeRegression

RidgeRegressionSVR(kernel="linear")

SVR(kernel="linea...

NO

NOSVR(kernel="rbf")

SVR(kernel="rbf...Ensemble

Regressors

Ensemble...

regression

regressionRandomized

PCA

Randomized...

YES

YES

<10K

samples

<10K...Kernel

Approximation

Kernel...

NO

NOIsoMap

IsoMapSpectral

Embedding

Spectral...

YES

YESLLE

LLE

dimensionality

reduction

dimensionality...

scikit-learn

algorithm cheat sheet

scikit-learn...

TRY

NEXT

TRY...

TRY

NEXT

TRY...

TRY

NEXT

TRY...

TRY

NEXT

TRY...

TRY

NEXT

TRY...

TRY

NEXT

TRY...

TRY

NEXT

TRY...Text is not SVG - cannot display

[previous

12.1. Array API support (experimental)](modules/array_api.html "previous page")
[next

14. External Resources, Videos and Talks](presentations.html "next page")