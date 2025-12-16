Source URL: https://www.huber.embl.de/msmb/12-chap.html
Date Scraped: 2025-12-12

---

# 12  Supervised Learning

[![](imgs/BuildWall.png)](https://www.huber.embl.de/msmb/imgs/BuildWall.png)

[![](imgs/EWall.png)](https://www.huber.embl.de/msmb/imgs/EWall.png)

In a supervised learning setting, we have a yardstick or plumbline to judge how well we are doing: the response itself.

A frequent question in biological and biomedical applications is whether a property of interest (say, disease type, cell type, the prognosis of a patient) can be “predicted”, given one or more other properties, called the **predictors**. Often we are motivated by a situation in which the property to be predicted is unknown (it lies in the future, or is hard to measure), while the predictors are known. The crucial point is that we *learn* the prediction rule from a set of *training data* in which the property of interest is also known. Once we have the rule, we can either apply it to new data, and make actual predictions of unknown outcomes; or we can dissect the rule with the aim of better understanding the underlying biology.

Compared to unsupervised learning and what we have seen in Chapters [5](https://www.huber.embl.de/msmb/05-chap.html), [7](https://www.huber.embl.de/msmb/07-chap.html) and [9](https://www.huber.embl.de/msmb/09-chap.html), where we do not know what we are looking for or how to decide whether our result is “right”, we are on much more solid ground with supervised learning: the objective is clearly stated, and there are straightforward criteria to measure how well we are doing.

The central issues in **supervised learning**[1](https://www.huber.embl.de/msmb/12-chap.html#fn1) are **overfitting** and **generalizability**: did we just learn the training data “by heart” by constructing a rule that has 100% accuracy on the training data, but would perform poorly on any new data? Or did our rule indeed pick up some of the pertinent patterns in the system being studied, which will also apply to yet unseen new data? ([Figure 12.1](https://www.huber.embl.de/msmb/12-chap.html#fig-overfitting-1))

1 Sometimes the term **statistical learning** is used, more or less exchangeably.

[![](12-chap_files/figure-html/fig-overfitting-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-overfitting-1-1.png "Figure 12.1: An example for overfitting: two regression lines are fit to data in the (x, y)-plane (black points). We can think of such a line as a rule that predicts the y-value, given an x-value. Both lines are smooth, but the fits differ in what is called their bandwidth, which intuitively can be interpreted their stiffness. The blue line seems overly keen to follow minor wiggles in the data, while the orange line captures the general trend but is less detailed. The effective number of parameters needed to describe the blue line is much higher than for the orange line. Also, if we were to obtain additional data, it is likely that the blue line would do a worse job than the orange line in modeling the new data. We’ll formalize these concepts –training error and test set error– later in this chapter. Although exemplified here with line fitting, the concept applies more generally to prediction models.")

Figure 12.1: An example for **overfitting**: two regression lines are fit to data in the \((x, y)\)-plane (black points). We can think of such a line as a rule that predicts the \(y\)-value, given an \(x\)-value. Both lines are smooth, but the fits differ in what is called their **bandwidth**, which intuitively can be interpreted their stiffness. The blue line seems overly keen to follow minor wiggles in the data, while the orange line captures the general trend but is less detailed. The effective number of parameters needed to describe the blue line is much higher than for the orange line. Also, if we were to obtain additional data, it is likely that the blue line would do a **worse** job than the orange line in modeling the new data. We’ll formalize these concepts –training error and test set error– later in this chapter. Although exemplified here with line fitting, the concept applies more generally to prediction models.

## 12.1 Goals for this chapter

In this chapter we will:

* See exemplary applications that motivate the use of supervised learning methods.
* Learn what discriminant analysis does.
* Define measures of performance.
* Encounter the curse of dimensionality and see what overfitting is.
* Find out about regularization – in particular, penalization – and understand the concepts of generalizability and model complexity.
* See how to use cross-validation to tune parameters of algorithms.
* Discuss method hacking.

## 12.2 What are the data?

The basic data structure for both supervised and unsupervised learning is (at least conceptually) a dataframe, where each row corresponds to an object and the columns are different features (usually numerical values) of the objects[2](https://www.huber.embl.de/msmb/12-chap.html#fn2). While in unsupervised learning we aim to find (dis)similarity relationships between the objects based on their feature values (e.g., by clustering or ordination), in supervised learning we aim to find a mathematical function (or a computational algorithm) that predicts the value of one of the features from the other features. Many implementations require that there are no missing values, whereas other methods can be made to work with some amount of missing data.

2 This is a simplified description. Machine learning is a huge field, and lots of generalizations of this simple conceptual picture have been made. Already the construction of relevant features is an art by itself — we have seen examples with images of cells in [Chapter 11](https://www.huber.embl.de/msmb/11-chap.html), and more generally there are lots of possibilities to extract features from images, sounds, movies, free text, \(...\) Moreover, there is a variant of machine learning methods called **kernel methods** that do not need features at all; instead, kernel methods use distances or measures of similarity between objects. It may be easier, for instance, to define a measure of similarity between two natural language text objects than to find relevant numerical features to represent them. Kernel methods are beyond the scope of this book.

The feature that we select over all the others with the aim of predicting is called the **objective** or the **response**. Sometimes the choice is natural, but sometimes it is also instructive to reverse the roles, especially if we are interested in dissecting the prediction function for the purpose of biological understanding, or in disentangling correlations from causation.

The framework for supervised learning covers both continuous and categorical response variables. In the continuous case we also call it **regression**, in the categorical case, **classification**. It turns out that this distinction is not a detail, as it has quite far-reaching consequences for the choice of loss function ([Section 12.5](https://www.huber.embl.de/msmb/12-chap.html#sec-supervised-objective)) and thus the choice of algorithm ([Friedman 1997](https://www.huber.embl.de/msmb/16-chap.html#ref-friedmanbiasvariance01)).

The first question to consider in any supervised learning task is how the number of objects compares to the number of predictors. The more objects, the better, and much of the hard work in supervised learning has to do with overcoming the limitations of having a finite (and typically, too small) training set.

[![](imgs/fourquad.png)](https://www.huber.embl.de/msmb/imgs/fourquad.png "Figure 12.2: In supervised learning, we assign two different roles to our variables. We have labeled the explanatory variables X and the response variable(s) Y. There are also two different sets of observations: the training set X_\ell and Y_\ell and the test set X_v and Y_v. (The subscripts refer to alternative names for the two sets: “learning” and “validation”.)")

Figure 12.2: In supervised learning, we assign two different roles to our variables. We have labeled the explanatory variables \(X\) and the response variable(s) \(Y\). There are also two different sets of observations: the training set \(X\_\ell\) and \(Y\_\ell\) and the test set \(X\_v\) and \(Y\_v\). (The subscripts refer to alternative names for the two sets: “learning” and “validation”.)

Task

Give examples where we have encountered instances of supervised learning with a categorical response in this book.

### 12.2.1 Motivating examples

#### Predicting diabetes type

The `diabetes` dataset ([Reaven and Miller 1979](https://www.huber.embl.de/msmb/16-chap.html#ref-diabetes)) presents three different groups of diabetes patients and five clinical variables measured on them.

```
data("diabetes", package = "rrcov")
head(diabetes)
```

```
    rw fpg glucose insulin sspg  group
1 0.81  80     356     124   55 normal
2 0.95  97     289     117   76 normal
3 0.94 105     319     143  105 normal
4 1.04  90     356     199  108 normal
5 1.00  90     323     240  143 normal
6 0.76  86     381     157  165 normal
```

The univariate distributions (more precisely, some density estimates of them) are shown in [Figure 12.3](https://www.huber.embl.de/msmb/12-chap.html#fig-ldagroups-1).

```
library("reshape2")
ggplot(melt(diabetes, id.vars = "group"), aes(x = value, col = group)) +
 geom_density() + facet_wrap( ~variable, ncol = 1, scales = "free") +
 theme(legend.position = "bottom")
```

[![](12-chap_files/figure-html/fig-ldagroups-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-ldagroups-1-1.png "Figure 12.3: We see already from the one-dimensional distributions that some of the individual variables could potentially predict which group a patient is more likely to belong to. Our goal is to combine variables to improve over such one-dimensional prediction models.")

Figure 12.3: We see already from the one-dimensional distributions that some of the individual variables could potentially predict which group a patient is more likely to belong to. Our goal is to combine variables to improve over such one-dimensional prediction models.

The variables are explained in the manual page of the dataset, and in the paper ([Reaven and Miller 1979](https://www.huber.embl.de/msmb/16-chap.html#ref-diabetes)):

* rw: relative weight
* fpg: fasting plasma glucose
* glucose: area under the plasma glucose curve for the three hour oral glucose tolerance test (OGTT)
* insulin: area under the plasma insulin curve for the OGTT
* sspg: steady state plasma glucose response
* group: normal, chemical diabetes and overt diabetes

#### Predicting cellular phenotypes

Neumann et al. ([2010](https://www.huber.embl.de/msmb/16-chap.html#ref-Neumann:2010)) observed human cancer cells using live-cell imaging. The cells were genetically engineered so that their histones were tagged with a green fluorescent protein (GFP). A genome-wide RNAi library was applied to the cells, and for each siRNA perturbation, movies of a few hundred cells were recorded for about two days, to see what effect the depletion of each gene had on cell cycle, nuclear morphology and cell proliferation. Their paper reports the use of an automated image classification algorithm that quantified the visual appearance of each cell’s nucleus and enabled the prediction of normal mitosis states or aberrant nuclei. The algorithm was trained on the data from around 3000 cells that were annotated by a human expert. It was then applied to almost 2 billions images of nuclei ([Figure 12.4](https://www.huber.embl.de/msmb/12-chap.html#fig-cellshape)). Using automated image classification provided scalability (annotating 2 billion images manually would take a long time) and objectivity.

[![](imgs/Neumann2010Fig1b_web.jpg)](https://www.huber.embl.de/msmb/imgs/Neumann2010Fig1b_web.jpg "Figure 12.4: The data were images of 2\times10^9 nuclei from movies. The images were segmented to identify the nuclei, and numeric features were computed for each nucleus, corresponding to size, shape, brightness and lots of other more or less abstract quantitative summaries of the joint distribution of pixel intensities. From the features, the cells were classified into 16 different nuclei morphology classes, represented by the rows of the barplot. Representative images for each class are shown in black and white in the center column. The class frequencies, which are very unbalanced, are shown by the lengths of the bars.")

Figure 12.4: The data were images of \(2\times10^9\) nuclei from movies. The images were segmented to identify the nuclei, and numeric features were computed for each nucleus, corresponding to size, shape, brightness and lots of other more or less abstract quantitative summaries of the joint distribution of pixel intensities. From the features, the cells were classified into 16 different nuclei morphology classes, represented by the rows of the barplot. Representative images for each class are shown in black and white in the center column. The class frequencies, which are very unbalanced, are shown by the lengths of the bars.

#### Predicting embryonic cell states

We will revisit the mouse embryo data ([Ohnishi et al. 2014](https://www.huber.embl.de/msmb/16-chap.html#ref-Ohnishi2014)), which we have already seen in Chapters [3](https://www.huber.embl.de/msmb/03-chap.html), [5](https://www.huber.embl.de/msmb/05-chap.html) and [7](https://www.huber.embl.de/msmb/07-chap.html). We’ll try to predict cell state and genotype from the gene expression measurements in Sections [12.3.2](https://www.huber.embl.de/msmb/12-chap.html#sec-supervised-predictembryonicday) and [12.6.3](https://www.huber.embl.de/msmb/12-chap.html#sec-supervised-classifymousecells).

## 12.3 Linear discrimination

We start with one of the simplest possible discrimination problems[3](https://www.huber.embl.de/msmb/12-chap.html#fn3): we have objects described by two continuous features (so the objects can be thought of as points in the 2D plane) and falling into three groups. Our aim is to define class boundaries, which are lines in the 2D space.

3 Arguably the simplest possible problem is a single continuous feature, two classes, and the task of finding a single threshold to discriminate between the two groups – as in [Figure 6.2](https://www.huber.embl.de/msmb/06-chap.html#fig-testing-FDRvspstatic1).

### 12.3.1 Diabetes data

Let’s see whether we can predict the `group` from the `sspg` and `glucose` variables in the `diabetes` data. It’s always a good idea to first visualise the data ([Figure 12.5](https://www.huber.embl.de/msmb/12-chap.html#fig-scatterdiabetes-1)).

```
ggdb = ggplot(mapping = aes(x = sspg, y = glucose)) +
  geom_point(aes(colour = group), data = diabetes)
ggdb
```

[![](12-chap_files/figure-html/fig-scatterdiabetes-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-scatterdiabetes-1-1.png "Figure 12.5: Scatterplot of two of the variables in the diabetes data. Each point is a sample, and the color indicates the diabetes type as encoded in the group variable.")

Figure 12.5: Scatterplot of two of the variables in the `diabetes` data. Each point is a sample, and the color indicates the diabetes type as encoded in the `group` variable.

We’ll start with a method called **linear discriminant analysis** (**LDA**). This method is a foundation stone of classification, many of the more complicated (and sometimes more powerful) algorithms are really just generalizations of LDA.

```
library("MASS")
diabetes_lda = lda(group ~ sspg + glucose, data = diabetes)
diabetes_lda
```

```
Call:
lda(group ~ sspg + glucose, data = diabetes)

Prior probabilities of groups:
   normal  chemical     overt 
0.5241379 0.2482759 0.2275862 

Group means:
             sspg   glucose
normal   114.0000  349.9737
chemical 208.9722  493.9444
overt    318.8788 1043.7576

Coefficients of linear discriminants:
                LD1         LD2
sspg    0.005036943 -0.01539281
glucose 0.005461400  0.00449050

Proportion of trace:
   LD1    LD2 
0.9683 0.0317
```

```
ghat = predict(diabetes_lda)$class
table(ghat, diabetes$group)
```

```
ghat       normal chemical overt
  normal       69       12     1
  chemical      7       24     6
  overt         0        0    26
```

```
mean(ghat != diabetes$group)
```

```
[1] 0.1793103
```

Question 12.1

What do the different parts of the above output mean?

Now, let’s visualise the LDA result. We are going to plot the prediction regions for each of the three groups. We do this by creating a grid of points and using our prediction rule on each of them. We’ll then also dig a bit deeper into the mechanics of LDA and plot the class centers (`diabetes_lda$means`) and ellipses that correspond to the fitted covariance matrix (`diabetes_lda$scaling`). Assembling this visualization requires us to write a bit of code.

```
make1Dgrid = function(x) {
  rg = grDevices::extendrange(x)
  seq(from = rg[1], to = rg[2], length.out = 100)
}
```

Set up the points for prediction, a \(100 \times 100\) grid that covers the data range.

```
diabetes_grid = with(diabetes,
  expand.grid(sspg = make1Dgrid(sspg),
              glucose = make1Dgrid(glucose)))
```

Do the predictions.

```
diabetes_grid$ghat =
  predict(diabetes_lda, newdata = diabetes_grid)$class
```

The group centers.

```
centers = diabetes_lda$means
```

Compute the ellipse. We start from a unit circle (approximated by a polygon with 360 sides) and apply the corresponding affine transformation from the LDA output.

```
unitcircle = exp(1i * seq(0, 2*pi, length.out = 360)) |>
          (\(z) cbind(Re(z), Im(z)))() 
ellipse = unitcircle %*% solve(diabetes_lda$scaling) |> as_tibble()
```

All three ellipses, one for each group center.

```
library("dplyr")
ellipses = lapply(rownames(centers), function(gr) {
  mutate(ellipse,
     sspg    = sspg    + centers[gr, "sspg"],
     glucose = glucose + centers[gr, "glucose"],
     group   = gr)
}) |> bind_rows()
```

Now we are ready to plot ([Figure 12.6](https://www.huber.embl.de/msmb/12-chap.html#fig-modeldiabetes-1)).

```
ggdb + geom_raster(aes(fill = ghat),
            data = diabetes_grid, alpha = 0.25, interpolate = TRUE) +
    geom_point(data = as_tibble(centers), pch = "+", size = 8) +
    geom_path(aes(colour = group), data = ellipses) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
```

[![](12-chap_files/figure-html/fig-modeldiabetes-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-modeldiabetes-1-1.png "Figure 12.6: As Figure fig-scatterdiabetes-1, with the classification regions from the LDA model shown. The three ellipses represent the class centers and the covariance matrix of the LDA model; note that there is only one covariance matrix, which is the same for all three classes. Therefore also the sizes and orientations of the ellipses are the same for the three classes, only their centers differ. They represent contours of equal class membership probability.")

Figure 12.6: As [Figure 12.5](https://www.huber.embl.de/msmb/12-chap.html#fig-scatterdiabetes-1), with the classification regions from the LDA model shown. The three ellipses represent the class centers and the covariance matrix of the LDA model; note that there is only one covariance matrix, which is the same for all three classes. Therefore also the sizes and orientations of the ellipses are the same for the three classes, only their centers differ. They represent contours of equal class membership probability.

Question 12.2

Why is the boundary between the prediction regions for *chemical* and *overt* not perpendicular to the line between the group centers?

Solution

The boundaries would be perpendicular if the ellipses were circles. In general, a boundary is tangential to the contours of equal class probabilities, and due the elliptic shape of the contours, a boundary is in general not perpendicular to the line between centers.

Question 12.3

How confident would you be about the predictions in those areas of the 2D plane that are far from all of the cluster centers?

Solution

Predictions that are far from any cluster center should be assessed critically, as this amounts to an extrapolation into regions where the LDA model may not be very good and/or there may be no training data nearby to support the prediction. We could use the distance to the nearest center as a measure of confidence in the prediction for any particular point; although we will see that resampling and cross-validation based methods offer more generic and usually more reliable measures.

Question 12.4

Why is the boundary between the prediction regions for *normal* and *chemical* not half-way between the centers, but shifted in favor of *normal*? Hint: have a look at the `prior` argument of `lda`. Try again with uniform prior.

Solution

The result of the following code chunk is shown in [Figure 12.7](https://www.huber.embl.de/msmb/12-chap.html#fig-diabetes-lda-uniform-prior-1). The suffix `_up` is short for “uniform prior”.

```
diabetes_up = lda(group ~ sspg + glucose, data = diabetes,
  prior = (\(n) rep(1/n, n)) (nlevels(diabetes$group)))

diabetes_grid$ghat_up =
  predict(diabetes_up, newdata = diabetes_grid)$class

stopifnot(all.equal(diabetes_up$means, diabetes_lda$means))

ellipse_up  = unitcircle %*% solve(diabetes_up$scaling) |> as_tibble()
ellipses_up = lapply(rownames(centers), function(gr) {
  mutate(ellipse_up,
     sspg    = sspg    + centers[gr, "sspg"],
     glucose = glucose + centers[gr, "glucose"],
     group   = gr)
}) |> bind_rows()

ggdb + geom_raster(aes(fill = ghat_up),
            data = diabetes_grid, alpha = 0.4, interpolate = TRUE) +
    geom_point(data = data.frame(centers), pch = "+", size = 8) +
    geom_path(aes(colour = group), data = ellipses_up) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
```

[![](12-chap_files/figure-html/fig-diabetes-lda-uniform-prior-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-diabetes-lda-uniform-prior-1-1.png "Figure 12.7: As Figure fig-modeldiabetes-1, but with uniform class priors.")

Figure 12.7: As [Figure 12.6](https://www.huber.embl.de/msmb/12-chap.html#fig-modeldiabetes-1), but with uniform class priors.

The `stopifnot` line confirms that the class centers are the same, as they are independent of the prior. The joint covariance is not.

Question 12.5

Figures [12.6](https://www.huber.embl.de/msmb/12-chap.html#fig-modeldiabetes-1) and [12.7](https://www.huber.embl.de/msmb/12-chap.html#fig-diabetes-lda-uniform-prior-1) show both the fitted LDA model, through the ellipses, and the prediction regions, through the area coloring. What part of this visualization is generic for all sorts of classification methods, what part is method-specific?

Solution

The prediction regions can be shown for any classification method, including a “black box” method. The cluster centers and ellipses in Figures [12.6](https://www.huber.embl.de/msmb/12-chap.html#fig-modeldiabetes-1) and [12.7](https://www.huber.embl.de/msmb/12-chap.html#fig-diabetes-lda-uniform-prior-1) are method-specific.

Question 12.6

What is the difference in the prediction accuracy if we use all 5 variables instead of just `glucose` and `sspg`?

Solution

```
diabetes_lda5 = lda(group ~ rw + fpg + glucose + sspg + insulin, data = diabetes)
diabetes_lda5
```

```
Call:
lda(group ~ rw + fpg + glucose + sspg + insulin, data = diabetes)

Prior probabilities of groups:
   normal  chemical     overt 
0.5241379 0.2482759 0.2275862 

Group means:
                rw       fpg   glucose     sspg  insulin
normal   0.9372368  91.18421  349.9737 114.0000 172.6447
chemical 1.0558333  99.30556  493.9444 208.9722 288.0000
overt    0.9839394 217.66667 1043.7576 318.8788 106.0000

Coefficients of linear discriminants:
                  LD1          LD2
rw       1.3624356881 -3.784142444
fpg     -0.0336487883  0.036633317
glucose  0.0125763942 -0.007092017
sspg     0.0042431866  0.001134070
insulin -0.0001022245 -0.006173424

Proportion of trace:
   LD1    LD2 
0.8812 0.1188
```

```
ghat5 = predict(diabetes_lda5)$class
table(ghat5, diabetes$group)
```

```
ghat5      normal chemical overt
  normal       73        5     1
  chemical      3       31     5
  overt         0        0    27
```

```
mean(ghat5 != diabetes$group)
```

```
[1] 0.09655172
```

Question 12.7

Instead of approximating the prediction regions by classification from a grid of points, compute the separating lines explicitly from the linear determinant coefficients.

Solution

See Section 4.3, Equation (4.10) in ([Hastie, Tibshirani, and Friedman 2008](https://www.huber.embl.de/msmb/16-chap.html#ref-HastieTibshiraniFriedman)).

### 12.3.2 Predicting embryonic cell state from gene expression

Assume that we already know that the four genes *FN1*, *TIMD2*, *GATA4* and *SOX7* are relevant to the classification task[4](https://www.huber.embl.de/msmb/12-chap.html#fn4). We want to build a classifier that predicts the developmental time (embryonic days: E3.25, E3.5, E4.5). We load the data and select four corresponding probes.

4 Later in this chapter we will see methods that can drop this assumption and screen all available features.

```
library("Hiiragi2013")
```

```
In chunk 'loadHiiragi2': Warning: replacing previous import 'boot::logit' by 'gtools::logit' whenloading 'Hiiragi2013'
```

```
In chunk 'loadHiiragi2': Warning: replacing previous import 'boot::inv.logit' by 'gtools::inv.logit'when loading 'Hiiragi2013'
```

```
data("x")
probes = c("1426642_at", "1418765_at", "1418864_at", "1416564_at")
embryoCells = t(Biobase::exprs(x)[probes, ]) |> as_tibble() |>
  mutate(Embryonic.day = x$Embryonic.day) |>
  dplyr::filter(x$genotype == "WT")
```

We can use the Bioconductor annotation package associated with the microarray to verify that the probes correspond to the intended genes.

```
annotation(x)
```

```
[1] "mouse4302"
```

```
library("mouse4302.db")
anno = AnnotationDbi::select(mouse4302.db, keys = probes,
                             columns = c("SYMBOL", "GENENAME"))
anno
```

```
     PROBEID SYMBOL                                            GENENAME
1 1426642_at    Fn1                                       fibronectin 1
2 1418765_at  Timd2 T cell immunoglobulin and mucin domain containing 2
3 1418864_at  Gata4                              GATA binding protein 4
4 1416564_at   Sox7                SRY (sex determining region Y)-box 7
```

```
mt = match(anno$PROBEID, colnames(embryoCells))
colnames(embryoCells)[mt] = anno$SYMBOL
```

Now we are ready to visualize the data in a pairs plot ([Figure 12.8](https://www.huber.embl.de/msmb/12-chap.html#fig-HiiragiFourGenesPairs-1)).

```
library("GGally")
ggpairs(embryoCells, mapping = aes(col = Embryonic.day),
  columns = anno$SYMBOL, upper = list(continuous = "points"))
```

[![](12-chap_files/figure-html/fig-HiiragiFourGenesPairs-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-HiiragiFourGenesPairs-1-1.png "Figure 12.8: Expression values of the discriminating genes, with the prediction target Embryonic.day shown by color.")

Figure 12.8: Expression values of the discriminating genes, with the prediction target Embryonic.day shown by color.

We can now call `lda` on these data. The linear combinations `LD1` and `LD2` that serve as discriminating variables are given in the slot `ed_lda$scaling` of the output from `lda`.

```
ec_lda = lda(Embryonic.day ~ Fn1 + Timd2 + Gata4 + Sox7,
             data = embryoCells)
round(ec_lda$scaling, 1)
```

```
       LD1  LD2
Fn1   -0.2  0.4
Timd2  0.5  0.0
Gata4 -0.1  0.6
Sox7  -0.7 -0.5
```

For the visualization of the learned model in [Figure 12.9](https://www.huber.embl.de/msmb/12-chap.html#fig-edcontour-1), we need to build the prediction regions and their boundaries by expanding the grid in the space of the two new coordinates LD1 and LD2.

```
ec_rot = predict(ec_lda)$x |> as_tibble() |>
           mutate(ed = embryoCells$Embryonic.day)
ec_lda2 = lda(ec_rot[, 1:2], predict(ec_lda)$class)
ec_grid = with(ec_rot, expand.grid(
  LD1 = make1Dgrid(LD1),
  LD2 = make1Dgrid(LD2)))
ec_grid$edhat = predict(ec_lda2, newdata = ec_grid)$class
ggplot() +
  geom_point(aes(x = LD1, y = LD2, colour = ed), data = ec_rot) +
  geom_raster(aes(x = LD1, y = LD2, fill = edhat),
            data = ec_grid, alpha = 0.4, interpolate = TRUE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed()
```

[![](12-chap_files/figure-html/fig-edcontour-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-edcontour-1-1.png "Figure 12.9: LDA classification regions for Embryonic.day.")

Figure 12.9: LDA classification regions for Embryonic.day.

Question 12.8

Repeat these analyses using quadratic discriminant analysis (`qda`). What difference do you see in the shape of the boundaries?

Solution

See code below and [Figure 12.10](https://www.huber.embl.de/msmb/12-chap.html#fig-qdamouse-1).

```
library("gridExtra")

ec_qda = qda(Embryonic.day ~ Fn1 + Timd2 + Gata4 + Sox7,
             data = embryoCells)

variables = colnames(ec_qda$means)
pairs = combn(variables, 2)
lapply(seq_len(ncol(pairs)), function(i) {
  grid = with(embryoCells,
    expand.grid(x = make1Dgrid(get(pairs[1, i])),
                y = make1Dgrid(get(pairs[2, i])))) |>
    `colnames<-`(pairs[, i])

  for (v in setdiff(variables, pairs[, i]))
    grid[[v]] = median(embryoCells[[v]])

  grid$edhat = predict(ec_qda, newdata = grid)$class

  x <- pairs[1,i]
  y <- pairs[2,i]
  ggplot() + 
    geom_point(
      data = embryoCells,
      aes(x = .data[[x]], y = .data[[y]], colour = Embryonic.day)
    ) +
    geom_raster(
      aes(x = .data[[x]], y = .data[[y]], fill = edhat),
      data = grid, alpha = 0.4, interpolate = TRUE
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    if (i != ncol(pairs)) theme(legend.position = "none")
}) |> (\(g) grid.arrange(grobs = g, ncol = 2))()
```

[![](12-chap_files/figure-html/fig-qdamouse-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-qdamouse-1-1.png "Figure 12.10: QDA for the mouse cell data. Shown are all pairwise plots of the four features. In each plot, the other two features are set to the median.")

Figure 12.10: QDA for the mouse cell data. Shown are all pairwise plots of the four features. In each plot, the other two features are set to the median.

Question 12.9

What happens if you call `lda` or `qda` with a lot more genes, say the first 1000, in the Hiiragi dataset?

Solution

```
lda(t(Biobase::exprs(x))[, 1:1000], x$Embryonic.day)
```

```
warnings()
qda(t(Biobase::exprs(x))[, 1:1000], x$Embryonic.day)
```

```
Error in qda.default(x, grouping, ...): some group is too small for 'qda'
```

The `lda` function manages to fit a model, but complains (with the warning) about the fact that there are more variables than replicates, which means that the variables are not linearly independent, and thus are redundant of each other. The `qda` function aborts with an error, since the QDA model with so many parameters cannot be fitted from the available data (at least, without making further assumptions, such as some sort of regularization, which it is not equipped for).

## 12.4 Machine learning vs rote learning

Computers are really good at memorizing facts. In the worst case, a machine learning algorithm is a roundabout way of doing this[5](https://www.huber.embl.de/msmb/12-chap.html#fn5). The central goal in statistical learning, however, is *generalizability*. We want an algorithm that is able to generalize, i.e., interpolate and extrapolate from given data to make good predictions about future data.

5 The not-so roundabout way is database technologies.

Let’s look at the following example. We generate random data (`rnorm`) for `n` objects, with different numbers of features (given by `p`). We train a LDA on these data and compute the **misclassification rate**, i.e., the fraction of times the prediction is wrong (`pred != resp`).

```
p = 2:21
n = 20

mcl = lapply(p, function(pp) {
  replicate(100, {
    xmat = matrix(rnorm(n * pp), nrow = n)
    resp = sample(c("apple", "orange"), n, replace = TRUE)
    fit  = lda(xmat, resp)
    pred = predict(fit)$class
    mean(pred != resp)
  }) |> mean() |> (\(x) tibble(mcl = x, p = pp))()
}) |> bind_rows()
```

```
ggplot(mcl, aes(x = p, y = mcl)) + 
  geom_line() + geom_point() +
  ylab("Misclassification rate")
```

[![](12-chap_files/figure-html/fig-learnbyheart-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-learnbyheart-1-1.png "Figure 12.11: Misclassification rate of LDA applied to random data. While the number of observations n is held constant (at 20), we are increasing the number of features p starting from 2 up to 21. The misclassification rate becomes almost zero as p approaches 20. The LDA model becomes so elaborate and over-parameterized that it manages to learn the random labels “by heart”. (As p becomes even larger, the “performance” degrades again somewhat, apparently due to numerical properties of the lda implementation used here.)")

Figure 12.11: Misclassification rate of LDA applied to random data. While the number of observations `n` is held constant (at 20), we are increasing the number of features `p` starting from 2 up to 21. The misclassification rate becomes almost zero as `p` approaches 20. The LDA model becomes so elaborate and over-parameterized that it manages to learn the random labels “by heart”. (As `p` becomes even larger, the “performance” degrades again somewhat, apparently due to numerical properties of the `lda` implementation used here.)

Question 12.10

What is the purpose of the `replicate` loop in the above code? What happens if you omit it (or replace the 100 by 1)?

Solution

For each single replicate, the curve is a noisier version of [Figure 12.11](https://www.huber.embl.de/msmb/12-chap.html#fig-learnbyheart-1). Averaging the measured misclassifications rate over 100 replicates makes the estimate more stable. We can do this since we are working with simulated data.

[Figure 12.11](https://www.huber.embl.de/msmb/12-chap.html#fig-learnbyheart-1) seems to imply that we can perfectly predict random labels from random data, if we only fit a complex enough model, i.e., one with many parameters. How can we overcome such an absurd conclusion? The problem with the above code is that the model performance is evaluated on the same data on which it was trained. This generally leads to positive bias, as you see in this crass example. How can we overcome this problem? The key idea is to assess model performance on different data than those on which the model was trained.

### 12.4.1 Cross-validation

A naive approach might be to split the data in two halves, and use the first half for learning (“training”) and the second half for assessment (“testing”). It turns out that this is needlessly variable and needlessly inefficient. It is needlessly variable, since by splitting the data only once, our results can be quite affected by how the split happens to fall. It seems better to do the splitting many times, and average. This will give us more stable results. It is needlessly inefficient, since the performance of machine learning algorithms depends on the number of observations, and the performance measured on half the data is likely[6](https://www.huber.embl.de/msmb/12-chap.html#fn6) to be worse than what it is with all the data. For this reason, it is better to use unequal sizes of training and test data. In the extreme case, we’ll use as much as \(n-1\) observations for training, and the remaining one for testing. After we’ve done this likewise for all observations, we can average our performance metric. This is called **leave-one-out cross-validation**.

6 Unless we have such an excess of data that it doesn’t matter.

[![](imgs/book_icon.png)](https://www.huber.embl.de/msmb/imgs/book_icon.png)

See Chapter *Model Assessment and Selection* in the book by Hastie, Tibshirani, and Friedman ([2008](https://www.huber.embl.de/msmb/16-chap.html#ref-HastieTibshiraniFriedman)) for further discussion on these trade-offs.

An alternative is **\(k\)-fold cross-validation**, where the observations are repeatedly split into a training set of size of around \(n(k-1)/k\) and a test set of size of around \(n/k\). Both alternatives have pros and contras, and there is not a universally best choice. An advantage of leave-one-out is that the amount of data used for training is close to the maximally available data; this is especially important if the sample size is limiting and “every little matters” for the algorithm. A drawback of leave-one-out is that the training sets are all very similar, so they may not model sufficiently well the kind of sampling changes to be expected if a new dataset came along. For large \(n\), leave-one-out cross-validation can be needlessly time-consuming.

```
estimate_mcl_loocv = function(x, resp) {
  vapply(seq_len(nrow(x)), function(i) {
    fit  = lda(x[-i, ], resp[-i])
    ptrn = predict(fit, newdata = x[-i,, drop = FALSE])$class
    ptst = predict(fit, newdata = x[ i,, drop = FALSE])$class
    c(train = mean(ptrn != resp[-i]), test = (ptst != resp[i]))
  }, FUN.VALUE = numeric(2)) |> rowMeans() |> t() |> as_tibble()
}

xmat = matrix(rnorm(n * last(p)), nrow = n)
resp = sample(c("apple", "orange"), n, replace = TRUE)

mcl = lapply(p, function(k) {
  estimate_mcl_loocv(xmat[, 1:k], resp)
}) |> bind_rows() |> data.frame(p) |> melt(id.vars = "p")
```

```
ggplot(mcl, aes(x = p, y = value, col = variable)) + geom_line() +
  geom_point() + ylab("Misclassification rate")
```

[![](12-chap_files/figure-html/fig-mclcv-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-mclcv-1-1.png "Figure 12.12: Cross-validation: the misclassification rate of LDA applied to random data, when evaluated on test data that were not used for learning, hovers around 0.5 independent of p. The misclassification rate on the training data is also shown. It behaves similar to what we already saw in Figure fig-learnbyheart-1.")

Figure 12.12: Cross-validation: the misclassification rate of LDA applied to random data, when evaluated on test data that were not used for learning, hovers around 0.5 independent of `p`. The misclassification rate on the training data is also shown. It behaves similar to what we already saw in [Figure 12.11](https://www.huber.embl.de/msmb/12-chap.html#fig-learnbyheart-1).

The result is show in [Figure 12.12](https://www.huber.embl.de/msmb/12-chap.html#fig-mclcv-1).

Question 12.11

Why are the curves in [Figure 12.12](https://www.huber.embl.de/msmb/12-chap.html#fig-mclcv-1) more variable (“wiggly”) than in [Figure 12.11](https://www.huber.embl.de/msmb/12-chap.html#fig-learnbyheart-1)? How can you overcome this?

Solution

Only one dataset (`xmat`, `resp`) was used to calculate [Figure 12.12](https://www.huber.embl.de/msmb/12-chap.html#fig-mclcv-1), whereas for [Figure 12.11](https://www.huber.embl.de/msmb/12-chap.html#fig-learnbyheart-1), we had the data generated within a `replicate` loop. You could similarly extend the above code to average the misclassification rate curves over many replicate simulated datasets.

### 12.4.2 The curse of dimensionality

In [Section 12.4.1](https://www.huber.embl.de/msmb/12-chap.html#sec-supervised-xval) we have seen overfitting and cross-validation on random data, but how does it look if there is in fact a relevant class separation?

```
p   = 2:20
mcl = replicate(100, {
  xmat = matrix(rnorm(n * last(p)), nrow = n)
  resp = sample(c("apple", "orange"), n, replace = TRUE)
  xmat[, 1:6] = xmat[, 1:6] + as.integer(factor(resp))

  lapply(p, function(k) {
    estimate_mcl_loocv(xmat[, 1:k], resp)
  }) |> bind_rows() |> cbind(p = p) |> melt(id.vars = "p")
}, simplify = FALSE) |> bind_rows()
```

```
mcl = group_by(mcl, p, variable) |> summarise(value = mean(value))

ggplot(mcl, aes(x = p, y = value, col = variable)) + geom_line() +
   geom_point() + ylab("Misclassification rate")
```

[![](12-chap_files/figure-html/fig-curseofdim-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-curseofdim-1.png "Figure 12.13: As we increase the number of features included in the model, the misclassification rate initially improves; as we start including more and more irrelevant features, it increases again, as we are fitting noise.")

Figure 12.13: As we increase the number of features included in the model, the misclassification rate initially improves; as we start including more and more irrelevant features, it increases again, as we are fitting noise.

[![](imgs/BiasVarianceTradeoff.png)](https://www.huber.embl.de/msmb/imgs/BiasVarianceTradeoff.png "Figure 12.14: Idealized version of Figure fig-curseofdim, from @HastieTibshiraniFriedman. A recurrent goal in machine learning is finding the sweet spot in the variance <- bias trade-off.")

Figure 12.14: Idealized version of [Figure 12.13](https://www.huber.embl.de/msmb/12-chap.html#fig-curseofdim), from Hastie, Tibshirani, and Friedman ([2008](https://www.huber.embl.de/msmb/16-chap.html#ref-HastieTibshiraniFriedman)). A recurrent goal in machine learning is finding the sweet spot in the variance <- bias trade-off.

The result is shown in [Figure 12.13](https://www.huber.embl.de/msmb/12-chap.html#fig-curseofdim). The group centers are the vectors (in \(\mathbb{R}^{20}\)) given by the coordinates \((1, 1, 1, 1, 1, 1, 0, 0, 0, ...)\) (apples) and \((2, 2, 2, 2, 2, 2, 0, 0, 0, ...)\) (oranges), and the optimal decision boundary is the hyperplane orthogonal to the line between them. For \(p\) smaller than \(6\), the decision rule cannot reach this hyperplane – it is biased. As a result, the misclassification rate is suboptimal, and it decreases with \(p\). But what happens for \(p\) larger than \(6\)? The algorithm is, in principle, able to model the optimal hyperplane, and it should not be distracted by the additional features. The problem is that it is. The more additional features enter the dataset, the higher the probability that one or more of them happen to fall in a way that they *look like* good, discriminating features in the training data – only to mislead the classifier and degrade its performance in the test data. Shortly we’ll see how to use penalization to (try to) control this problem.

The term **curse of dimensionality** was coined by Bellman ([1961](https://www.huber.embl.de/msmb/16-chap.html#ref-Bellman:1961)). It refers to the fact that high-dimensional spaces are very hard, if not impossible, to sample thoroughly: for instance, to cover a 2-dimensional square of side length 1 with grid points that are 0.1 apart, we need \(10^2=100\) points. In 100 dimensions, we need \(10^{100}\) – which is more than [the number of protons in the universe](https://en.wikipedia.org/wiki/Eddington_number). In genomics, we often aim to fit models to data with thousands of features. Also our intuitions about distances between points or about the relationship between a volume and its surface break down in a high-dimensional settings. We’ll explore some of the weirdnesses of high-dimensional spaces in the next few questions.

Question 12.12

Assume you have a dataset with 1 000 000 data points in \(p\) dimensions. The data are uniformly distributed in the unit hybercube (i.e., all features lie in the interval \([0,1]\)). What’s the side length of a hybercube that can be expected to contain just 10 of the points, as a function of \(p\)?

Solution

See [Figure 12.15](https://www.huber.embl.de/msmb/12-chap.html#fig-cursedimans1-1).

```
sideLength = function(p, pointDensity = 1e6, pointsNeeded = 10)
  (pointsNeeded / pointDensity) ^ (1 / p)
ggplot(tibble(p = 1:400, sideLength = sideLength(p)),
       aes(x = p, y = sideLength)) + geom_line(col = "red") +
  geom_hline(aes(yintercept = 1), linetype = 2)
```

[![](12-chap_files/figure-html/fig-cursedimans1-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-cursedimans1-1-1.png "Figure 12.15: Side length of a p-dimensional hybercube expected to contain 10 points out of 1 million uniformly distributed ones, as a function of the p. While for p=1, this length is conveniently small, namely 10/10^6=10^{-5}, for larger p it approaches 1, i.,e., becomes the same as the range of each the features. This means that a “local neighborhood” of 10 points encompasses almost the same data range as the whole dataset.")

Figure 12.15: Side length of a \(p\)-dimensional hybercube expected to contain 10 points out of 1 million uniformly distributed ones, as a function of the \(p\). While for \(p=1\), this length is conveniently small, namely \(10/10^6=10^{-5}\), for larger \(p\) it approaches 1, i.,e., becomes the same as the range of each the features. This means that a “local neighborhood” of 10 points encompasses almost the same data range as the whole dataset.

Next, let’s look at the relation between inner regions of the feature space versus its boundary regions. Generally speaking, prediction at the boundaries of feature space is more difficult than in its interior, as it tends to involve extrapolation, rather than interpolation. In the next question you’ll see how this difficulty explodes with feature space dimension.

Question 12.13

What fraction of a unit cube’s total volume is closer than 0.01 to any of its surfaces, as a function of the dimension?

Solution

See code below and [Figure 12.16](https://www.huber.embl.de/msmb/12-chap.html#fig-cursedimans2-1).

```
tibble(
  p = 1:400,
  volOuterCube = 1 ^ p,
  volInnerCube = 0.98 ^ p,  # 0.98 = 1 - 2 * 0.01
  `V(shell)` = volOuterCube - volInnerCube) |>
ggplot(aes(x = p, y =`V(shell)`)) + geom_line(col = "blue")
```

[![](12-chap_files/figure-html/fig-cursedimans2-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-cursedimans2-1-1.png "Figure 12.16: Fraction of a unit cube’s total volume that is in its “shell” (here operationalised as those points that are closer than 0.01 to its surface) as a function of the dimension p.")

Figure 12.16: Fraction of a unit cube’s total volume that is in its “shell” (here operationalised as those points that are closer than 0.01 to its surface) as a function of the dimension \(p\).

Question 12.14

What is the coefficient of variation (ratio of standard deviation over average) of the distance between two randomly picked points in the unit hypercube, as a function of the dimension?

Solution

We solve this one by simulation. We generate `n` pairs of random points in the hypercube (`x1`, `x2`) and compute their Euclidean distances. See [Figure 12.17](https://www.huber.embl.de/msmb/12-chap.html#fig-cursedimans3-1). This result can also be predicted from the central limit theorem.

```
n = 1000
df = tibble(
  p = round(10 ^ seq(0, 4, by = 0.25)),
  cv = vapply(p, function(k) {
    x1 = matrix(runif(k * n), nrow = n)
    x2 = matrix(runif(k * n), nrow = n)
    d = sqrt(rowSums((x1 - x2)^2))
    sd(d) / mean(d)
  }, FUN.VALUE = numeric(1)))
ggplot(df, aes(x = log10(p), y = cv)) + geom_line(col = "orange") +
  geom_point()
```

[![](12-chap_files/figure-html/fig-cursedimans3-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-cursedimans3-1-1.png "Figure 12.17: Coefficient of variation (CV) of the distance between randomly picked points in the unit hypercube, as a function of the dimension. As the dimension increases, everybody is equally far away from everyone else: there is almost no variation in the distances any more.")

Figure 12.17: Coefficient of variation (CV) of the distance between randomly picked points in the unit hypercube, as a function of the dimension. As the dimension increases, everybody is equally far away from everyone else: there is almost no variation in the distances any more.

## 12.5 Objective functions

We’ve already seen the **misclassification rate** (MCR) used to assess our classification performance in Figures [12.11](https://www.huber.embl.de/msmb/12-chap.html#fig-learnbyheart-1)–[12.13](https://www.huber.embl.de/msmb/12-chap.html#fig-curseofdim). Its population version is defined as

\[
\text{MCR} = \text{E}\left[ 𝟙\_{\hat{y} \neq y} \right],
\text{MCR} = \text{E}\left[ 𝟙\_{\hat{y} \neq y} \right],
\tag{12.1}\]

and for a finite sample

\[
\widehat{\text{MCR}} = \frac{1}{n}\sum\_{i=1}^n 𝟙\_{\hat{y\_i} \neq y\_i}.
\tag{12.2}\]

This is not the only choice we could make. Perhaps we care more about the misclassification of apples as oranges than vice versa, and we can reflect this by introducing weights that depend on the type of error made into the sum of [Equation 12.2](https://www.huber.embl.de/msmb/12-chap.html#eq-sec-super-obsmcl) (or the integral of [Equation 12.1](https://www.huber.embl.de/msmb/12-chap.html#eq-sec-super-expmcl)). This can get even more elaborate if we have more than two classes. Often we want to see the whole **confusion table**, which we can get via

```
table(truth, response)
```

An important special case is binary classification with asymmetric costs – think about, say, a medical test. Here, the **sensitivity** (a.k.a. **true positive rate** or **recall**) is related to the misclassification of healthy as ill, and the **specificity** (or **true negative rate**) depends on the probability of misclassification of ill as healthy. Often, there is a single parameter (e.g., a threshold) that can be moved up and down, allowing a trade-off between sensitivity and specificity (and thus, equivalently, between the two types of misclassification). In those cases, we usually are not content to know the classifier performance at one single choice of threshold, but at many (or all) of them. This leads to **receiver operating characteristic** (**ROC**) or **precision-recall** curves.

Question 12.15

What are the exact relationships between the per-class misclassification rates and sensitivity and specificity?

Solution

The sensitivity or true positive rate is

\[
\text{TPR} = \frac{\text{TP}}{\text{P}},
\]

where \(\text{TP}\) is the number of true positives and \(\text{P}\) the number of all positives. The specificity or true negative rate is

\[
\text{SPC} = \frac{\text{TN}}{\text{N}},
\]

where \(\text{TN}\) is the number of true negatives and \(\text{N}\) the number of all negatives. See also <https://en.wikipedia.org/wiki/Sensitivity_and_specificity>

Another cost function can be computed from the **Jaccard index**, which we already saw in [Chapter 5](https://www.huber.embl.de/msmb/05-chap.html).

\[
J(A,B) = \frac{|\,A\cap B\,|}{|\,A\cup B\,|},
\tag{12.3}\]

where \(A\) is the set of observations for which the true class is 1 (\(A=\{i\,|\,y\_i=1\}\)) and \(B\) is the set of observations for which the predicted class is 1. The number \(J\) is between 0 and 1, and when \(J\) is large, it indicates high overlap of the two sets. Note that \(J\) does not depend on the number of observations for which both true and predicted class is 0 – so it is particularly suitable for measuring the performance of methods that try to find rare events.

We can also consider probabilistic class predictions, which come in the form \(\hat{P}(Y\,|\,X)\). In this case, a possible risk function would be obtained by looking at distances between the true probability distribution and the estimated probability distributions. For two classes, the finite sample version of the \(\log \text{loss}\) is

\[
\log \text{loss} = -\frac{1}{n}\sum\_{i=1}^n y\_i\log(\hat{p}\_i) + (1 - y\_i)\log(1 - \hat{p}\_i),
\tag{12.4}\]

where \(\hat{p}\_i \in [0,1]\) is the prediction, and \(y\_i\in\{0,1\}\) is the truth.

[![Note that the \log\text{loss} will be infinite if a prediction is totally confident (\hat{p}_i is exactly 0 or 1) but wrong.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "Note that the \log\text{loss} will be infinite if a prediction is totally confident (\hat{p}_i is exactly 0 or 1) but wrong.")

Note that the \(\log\text{loss}\) will be infinite if a prediction is totally confident (\(\hat{p}\_i\) is exactly \(0\) or \(1\)) but wrong.

For continuous continuous response variables (regression), a natural choice is the **mean squared error** (**MSE**). It is the average squared error,

\[
\widehat{\text{MSE}} = \frac{1}{n}\sum\_{i=1}^n ( \hat{Y}\_i - Y\_i )^2.
\tag{12.5}\]

The population version is defined analogously, by turning the summation into an integral as in Equations [12.1](https://www.huber.embl.de/msmb/12-chap.html#eq-sec-super-expmcl) and [12.2](https://www.huber.embl.de/msmb/12-chap.html#eq-sec-super-obsmcl).

Statisticians call functions like Equations [12.1](https://www.huber.embl.de/msmb/12-chap.html#eq-sec-super-expmcl)—[12.5](https://www.huber.embl.de/msmb/12-chap.html#eq-sec-super-obsmse) variously (and depending on context and predisposition) **risk function**, **cost function**, **objective function**[7](https://www.huber.embl.de/msmb/12-chap.html#fn7).

7 There is even an R package dedicated to evaluation of statistical learners called **[metrics](https://cran.r-project.org/web/packages/metrics/)**.

## 12.6 Variance–bias trade-off

[![](imgs/TargetBias.png)](https://www.huber.embl.de/msmb/imgs/TargetBias.png "Figure 12.18: In bull’s eye (a), the estimates are systematically off target, but in a quite reproducible manner. The green segment represents the bias. In bull’s eye (b), the estimates are not biased, as they are centered in the right place, however they have high variance. We can distinguish the two scenarios since we see the result from many shots. If we only had one shot and missed the bull’s eye, we could not easily know whether that’s because of bias or variance.")

[![](imgs/TargetVariance.png)](https://www.huber.embl.de/msmb/imgs/TargetVariance.png "Figure 12.18: In bull’s eye (a), the estimates are systematically off target, but in a quite reproducible manner. The green segment represents the bias. In bull’s eye (b), the estimates are not biased, as they are centered in the right place, however they have high variance. We can distinguish the two scenarios since we see the result from many shots. If we only had one shot and missed the bull’s eye, we could not easily know whether that’s because of bias or variance.")

Figure 12.18: In bull’s eye (a), the estimates are systematically off target, but in a quite reproducible manner. The green segment represents the bias. In bull’s eye (b), the estimates are not biased, as they are centered in the right place, however they have high variance. We can distinguish the two scenarios since we see the result from many shots. If we only had one shot and missed the bull’s eye, we could not easily know whether that’s because of bias or variance.

An important fact that helps us understand the tradeoffs when picking a statistical learning model is that the MSE is the sum of two terms, and often the choices we can make are such that one of those terms goes down while the other one goes up. The bias measures how different the average of all the different estimates is from the truth, and variance, how much an individual one might scatter from the average value ([Figure 12.18](https://www.huber.embl.de/msmb/12-chap.html#fig-supervised-bullseye)). In applications, we often only get one shot, therefore being reliably almost on target can beat being right on the long term average but really off today. The decomposition

\[
\text{MSE} = \underbrace{\text{Var}(\hat{Y})}\_{\text{variance}}
+ \underbrace{\mathbb{E}[\hat{Y}-Y]^2}\_{\text{bias}}
\tag{12.6}\]

follows by straightforward algebra.

When trying to minimize the MSE, it is important to realize that sometimes we can pay the price of a small bias to greatly reduce variance, and thus overall improve MSE. We already encountered shrinkage estimation in [Chapter 8](https://www.huber.embl.de/msmb/08-chap.html). In classification (i.e., when we have categorical response variables), different objective functions than the MSE are used, and there is usually no such straightforward decomposition as in [Equation 12.6](https://www.huber.embl.de/msmb/12-chap.html#eq-sec-super-biasvariancedecompose). The good news is that we can usually go even much further than in the case of continuous responses with our trading biases for variance. This is because the discreteness of the response absorbs certain biases ([Friedman 1997](https://www.huber.embl.de/msmb/16-chap.html#ref-friedmanbiasvariance01)), so that the cost of higher bias is almost zero, while we still get the benefit of better (smaller) variance.

### 12.6.1 Penalization

In high-dimensional statistics, we are constantly plagued by variance: there is just not enough data to fit all the possible parameters. One of the most fruitful ideas in high-dimensional statistics is **penalization**: a tool to actively control and exploit the variance-bias tradeoff. Penalization is part of a larger class of regularization methods that are used to ensure stable estimates.

Although generalization of LDA to high-dimensional settings is possible ([Clemmensen et al. 2011](https://www.huber.embl.de/msmb/16-chap.html#ref-clemmensen2012sparse); [Witten and Tibshirani 2011](https://www.huber.embl.de/msmb/16-chap.html#ref-witten2011penalized)), it turns out that logistic regression is a more general approach[8](https://www.huber.embl.de/msmb/12-chap.html#fn8), and therefore we’ll now switch to that, using the **[glmnet](https://cran.r-project.org/web/packages/glmnet/)** package.

8 It fits into the framework of generalized linear models, which we encountered in [Chapter 8](https://www.huber.embl.de/msmb/08-chap.html).

For multinomial—or, for the special case of two classes, binomial—logistic regression models, the posterior log-odds between \(k\) classes and can be written in the form (see the section on *Logistic Regression* in the book by Hastie, Tibshirani, and Friedman ([2008](https://www.huber.embl.de/msmb/16-chap.html#ref-HastieTibshiraniFriedman)) for a more complete presentation):

\[
\log \frac{P(Y=i\,|\,X=x)}{P(Y=k\,|\,X=x)} = \beta^0\_i + \beta\_i x,
\tag{12.7}\]

where \(i=1,...,k-1\) enumerates the different classes and the \(k\)-th class is chosen as a reference. The data matrix \(x\) has dimensions \(n\times p\), where \(n\) is number of observations and \(p\) the number of features. The \(p\)-dimensional vector \(\beta\_i\) determines how the classification odds for class \(i\) versus class \(k\) depend on \(x\). The numbers \(\beta^0\_i\) are intercepts and depend, among other things, on the classes’ prior probabilities. Instead of the log odds [12.7](https://www.huber.embl.de/msmb/12-chap.html#eq-super-deflogregress) (i.e., ratios of class probabilities), we can also write down an equivalent model for the class probabilities themselves, and the fact that we here used the \(k\)-th class as a reference is an arbitrary choice, as the model estimates are equivariant under this choice ([Hastie, Tibshirani, and Friedman 2008](https://www.huber.embl.de/msmb/16-chap.html#ref-HastieTibshiraniFriedman)). The model is fit by maximising the log-likelihood \(\mathcal{l}(\beta, \beta^0; x)\), where \(\beta=(\beta\_1,...,\beta\_{k-1})\) and analogously for \(\beta^0\).

So far, so good. But as \(p\) gets larger, there is an increasing chance that some of the estimates go wildly off the mark, due to random sampling happenstances in the data (remember [Figure 12.1](https://www.huber.embl.de/msmb/12-chap.html#fig-overfitting-1)). This is true even if for each individual coordinate of the vector \(\beta\_i\), the error distribution is bounded: the probabilty of there being one coordinate that is in the far tails increases the more coordiates there are, i.e., the larger \(p\) is.

A related problem can also occur, not in [12.7](https://www.huber.embl.de/msmb/12-chap.html#eq-super-deflogregress), but in other, non-linear models, as the model dimension \(p\) increases while the sample size \(n\) remains the same: the likelihood landscape around its maximum becomes increasingly flat, and the maximum-likelihood estimate of the model parameters becomes more and more variable. Eventually, the maximum is no longer a point, but a submanifold, and the maximum likelihood estimate is unidentifiable. Both of these limitations can be overcome with a modification of the objective: instead of maximising the bare log-likelihood, we maximise a penalized version of it,

\[
\hat{\beta}= \arg\max\_\beta \mathcal{l}(\beta, \beta^0; x) + \lambda \operatorname{pen}(\beta),
\tag{12.8}\]

where \(\lambda\ge0\) is a real number, and \(\operatorname{pen}\) is a convex function, called the **penalty function**. Popular choices are \(\operatorname{pen}(\beta)=|\beta|^2\) (**ridge regression**) and \(\operatorname{pen}(\beta)=|\beta|^1\) (**lasso**).

[![Here, |\beta|^\nu=\sum_i\beta_i^\nu is the L_\nu-norm of the vector \beta. Variations are possible, for instead we could include in this summation only some but not all of the elements of \beta; or we could scale different elements differently, for instance based on some prior belief of their scale and importance.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "Here, |\beta|^\nu=\sum_i\beta_i^\nu is the L_\nu-norm of the vector \beta. Variations are possible, for instead we could include in this summation only some but not all of the elements of \beta; or we could scale different elements differently, for instance based on some prior belief of their scale and importance.")

Here, \(|\beta|^\nu=\sum\_i\beta\_i^\nu\) is the \(L\_\nu\)-norm of the vector \(\beta\). Variations are possible, for instead we could include in this summation only some but not all of the elements of \(\beta\); or we could scale different elements differently, for instance based on some prior belief of their scale and importance.

In the **elastic net**, ridge and lasso are hybridized by using the penalty function \(\operatorname{pen}(\beta)=(1-\alpha)|\beta|^1+\alpha|\beta|^2\) with some further parameter \(\alpha\in[0,1]\). The crux is, of course, how to choose the right \(\lambda\), and we will discuss that in the following.

### 12.6.2 Example: predicting colon cancer from stool microbiome composition

Zeller et al. ([2014](https://www.huber.embl.de/msmb/16-chap.html#ref-Zeller:MSB:2014)) studied metagenome sequencing data from fecal samples of 156 humans that included colorectal cancer patients and tumor-free controls. Their aim was to see whether they could identify biomarkers (presence or abundance of certain taxa) that could help with early tumor detection. The data are available from [Bioconductor](https://www.bioconductor.org/) through its **ExperimentHub** service under the identifier EH361.

```
library("ExperimentHub")
eh = ExperimentHub()
zeller = eh[["EH361"]]
```

```
table(zeller$disease)
```

```
       cancer large_adenoma             n small_adenoma 
           53            15            61            27
```

Question 12.16

Explore the `eh` object to see what other datasets there are.

Solution

Type `eh` into the R prompt and study the output.

For the following, let’s focus on the normal and cancer samples and set the adenomas aside.

```
zellerNC = zeller[, zeller$disease %in% c("n", "cancer")]
```

Before jumping into model fitting, as always it’s a good idea to do some exploration of the data. First, let’s look at the sample annotations. The following code prints the data from three randomly picked samples. (Only looking at the first ones, say with the R function `head`, is also an option, but may not be representative of the whole dataset).

```
pData(zellerNC)[ sample(ncol(zellerNC), 3), ]
```

```
                   subjectID age gender bmi country disease tnm_stage
CCIS62794166ST-4-0    FR-780  68   male  23  france  cancer    t3n1m1
CCIS87167916ST-4-0    FR-223  79   male  30  france  cancer    t4n1m1
CCIS17669415ST-4-0    FR-125  72 female  37  france  cancer    t1n0m0
                   ajcc_stage localization     fobt wif-1_gene_methylation_test
CCIS62794166ST-4-0         iv           lc negative                        <NA>
CCIS87167916ST-4-0         iv       rectum positive                    negative
CCIS17669415ST-4-0          i           lc negative                    positive
                   group bodysite ethnicity number_reads
CCIS62794166ST-4-0   crc    stool     white     39333809
CCIS87167916ST-4-0   crc    stool     white     50204707
CCIS17669415ST-4-0   crc    stool     white     69983469
```

Next, let’s explore the feature names:

[![We define the helper function formatfn to line wrap these long character strings for the available space here.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "We define the helper function formatfn to line wrap these long character strings for the available space here.")

We define the helper function `formatfn` to line wrap these long character strings for the available space here.

```
formatfn = function(x)
   gsub("|", "| ", x, fixed = TRUE) |> lapply(strwrap)

rownames(zellerNC)[1:4]
```

```
[1] "k__Bacteria"                  "k__Viruses"                  
[3] "k__Bacteria|p__Firmicutes"    "k__Bacteria|p__Bacteroidetes"
```

```
rownames(zellerNC)[nrow(zellerNC) + (-2:0)] |> formatfn()
```

```
[[1]]
[1] "k__Bacteria| p__Proteobacteria| c__Deltaproteobacteria|"         
[2] "o__Desulfovibrionales| f__Desulfovibrionaceae| g__Desulfovibrio|"
[3] "s__Desulfovibrio_termitidis"                                     

[[2]]
[1] "k__Viruses| p__Viruses_noname| c__Viruses_noname| o__Viruses_noname|"
[2] "f__Baculoviridae| g__Alphabaculovirus|"                              
[3] "s__Bombyx_mori_nucleopolyhedrovirus|"                                
[4] "t__Bombyx_mori_nucleopolyhedrovirus_unclassified"                    

[[3]]
[1] "k__Bacteria| p__Proteobacteria| c__Deltaproteobacteria|"         
[2] "o__Desulfovibrionales| f__Desulfovibrionaceae| g__Desulfovibrio|"
[3] "s__Desulfovibrio_termitidis| t__GCF_000504305"
```

As you can see, the features are a mixture of abundance quantifications at different taxonomic levels, from *k*ingdom over *p*hylum to *s*pecies. We could select only some of these, but here we continue with all of them. Next, let’s look at the distribution of some of the features. Here, we show an arbitrary choice of two, number 510 and 527; in practice, it is helpful to scroll through many such plots quickly to get an impression ([Figure 12.19](https://www.huber.embl.de/msmb/12-chap.html#fig-zellerHist-1)).

```
ggplot(melt(Biobase::exprs(zellerNC)[c(510, 527), ]), aes(x = value)) +
    geom_histogram(bins = 25) +
    facet_wrap( ~ Var1, ncol = 1, scales = "free")
```

[![](12-chap_files/figure-html/fig-zellerHist-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-zellerHist-1-1.png "Figure 12.19: Histograms of the distributions for two randomly selected features. The distributions are highly skewed, with many zero values and a thin, long tail of non-zero values.")

Figure 12.19: Histograms of the distributions for two randomly selected features. The distributions are highly skewed, with many zero values and a thin, long tail of non-zero values.

In the simplest case, we fit model [12.7](https://www.huber.embl.de/msmb/12-chap.html#eq-super-deflogregress) as follows.

```
library("glmnet")
glmfit = glmnet(x = t(Biobase::exprs(zellerNC)),
                y = factor(zellerNC$disease),
                family = "binomial")
```

A remarkable feature of the `glmnet` function is that it fits [12.7](https://www.huber.embl.de/msmb/12-chap.html#eq-super-deflogregress) not only for one choice of \(\lambda\), but for all possible \(\lambda\)s at once. For now, let’s look at the prediction performance for, say, \(\lambda=0.04\). The name of the function parameter is `s`:

```
predTrsf = predict(glmfit, newx = t(Biobase::exprs(zellerNC)),
                   type = "class", s = 0.04)
table(predTrsf, zellerNC$disease)
```

```
predTrsf cancer  n
  cancer     51  0
  n           2 61
```

Not bad – but remember that this is on the training data, without cross-validation. Let’s have a closer look at `glmfit`. The **[glmnet](https://cran.r-project.org/web/packages/glmnet/)** package offers a a diagnostic plot that is worth looking at ([Figure 12.20](https://www.huber.embl.de/msmb/12-chap.html#fig-plotglmfit-1)).

```
plot(glmfit, col = RColorBrewer::brewer.pal(8, "Dark2"), lwd = sqrt(3), ylab = "")
```

[![](12-chap_files/figure-html/fig-plotglmfit-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-plotglmfit-1-1.png "Figure 12.20: Regularization paths for glmfit.")

Figure 12.20: Regularization paths for `glmfit`.

Question 12.17

What are the \(x\)- and \(y\)-axes in [Figure 12.20](https://www.huber.embl.de/msmb/12-chap.html#fig-plotglmfit-1)? What are the different lines?

Solution

Consult the manual page of the function `plot.glmnet` in the **[glmnet](https://cran.r-project.org/web/packages/glmnet/)** package.

Let’s get back to the question of how to choose the parameter \(\lambda\). We could try many different choices –and indeed, all possible choices– of \(\lambda\), assess classification performance in each case using cross-validation, and then choose the best \(\lambda\).

[![You’ll already realize from the description of this strategy that if we optimize \lambda in this way, the resulting apparent classification performance will likely be exaggerated. We need a truly independent dataset, or at least another, outer cross-validation loop to get a more realistic impression of the generalizability. We will get back to this question at the end of the chapter.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "You’ll already realize from the description of this strategy that if we optimize \lambda in this way, the resulting apparent classification performance will likely be exaggerated. We need a truly independent dataset, or at least another, outer cross-validation loop to get a more realistic impression of the generalizability. We will get back to this question at the end of the chapter.")

You’ll already realize from the description of this strategy that if we optimize \(\lambda\) in this way, the resulting apparent classification performance will likely be exaggerated. We need a truly independent dataset, or at least another, outer cross-validation loop to get a more realistic impression of the generalizability. We will get back to this question at the end of the chapter.

We could do so by writing a loop as we did in the `estimate_mcl_loocv` function in [Section 12.4.1](https://www.huber.embl.de/msmb/12-chap.html#sec-supervised-xval). It turns out that the **[glmnet](https://cran.r-project.org/web/packages/glmnet/)** package already has built-in functionality for that, with the function `cv.glmnet`, which we can use instead.

```
cvglmfit = cv.glmnet(x = t(Biobase::exprs(zellerNC)),
                     y = factor(zellerNC$disease),
                     family = "binomial")
plot(cvglmfit)
```

[![](12-chap_files/figure-html/fig-colonCV-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-colonCV-1-1.png "Figure 12.21: Diagnostic plot for cv.glmnet: shown is a measure of cross-validated prediction performance, the deviance, as a function of \lambda. The dashed vertical lines show lambda.min and lambda.1se.")

Figure 12.21: Diagnostic plot for `cv.glmnet`: shown is a measure of cross-validated prediction performance, the deviance, as a function of \(\lambda\). The dashed vertical lines show `lambda.min` and `lambda.1se`.

The diagnostic plot is shown in [Figure 12.21](https://www.huber.embl.de/msmb/12-chap.html#fig-colonCV-1). We can access the optimal value with

```
cvglmfit$lambda.min
```

```
[1] 0.06998238
```

As this value results from finding a minimum in an estimated curve, it turns out that it is often too small, i.e., that the implied penalization is too weak. A heuristic recommended by the authors of the **[glmnet](https://cran.r-project.org/web/packages/glmnet/)** package is to use a somewhat larger value instead, namely the largest value of \(\lambda\) such that the performance measure is within 1 standard error of the minimum.

```
cvglmfit$lambda.1se
```

```
[1] 0.09691764
```

Question 12.18

How does the confusion table look like for \(\lambda=\;\)`lambda.1se`?

Solution

```
s0 = cvglmfit$lambda.1se
predict(glmfit, newx = t(Biobase::exprs(zellerNC)),type = "class", s = s0) |>
    table(zellerNC$disease)
```

```
         cancer  n
  cancer     37  6
  n          16 55
```

Question 12.19

What features drive the classification?

Solution

```
coefs = coef(glmfit)[, which.min(abs(glmfit$lambda - s0))]
topthree = order(abs(coefs), decreasing = TRUE)[1:3]
as.vector(coefs[topthree])
```

```
[1] -41.884114  -6.396031  -1.234294
```

```
formatfn(names(coefs)[topthree])
```

```
[[1]]
[1] "k__Bacteria| p__Candidatus_Saccharibacteria|"      
[2] "c__Candidatus_Saccharibacteria_noname|"            
[3] "o__Candidatus_Saccharibacteria_noname|"            
[4] "f__Candidatus_Saccharibacteria_noname|"            
[5] "g__Candidatus_Saccharibacteria_noname|"            
[6] "s__candidate_division_TM7_single_cell_isolate_TM7b"

[[2]]
[1] "k__Bacteria| p__Firmicutes| c__Clostridia| o__Clostridiales|"        
[2] "f__Ruminococcaceae| g__Subdoligranulum| s__Subdoligranulum_variabile"

[[3]]
[1] "k__Bacteria| p__Firmicutes| c__Clostridia| o__Clostridiales|"
[2] "f__Lachnospiraceae| g__Lachnospiraceae_noname|"              
[3] "s__Lachnospiraceae_bacterium_7_1_58FAA"
```

Question 12.20

How do the results change if we transform the data, say, with the `asinh` transformation as we saw in [Chapter 5](https://www.huber.embl.de/msmb/05-chap.html)?

Solution

See [Figure 12.22](https://www.huber.embl.de/msmb/12-chap.html#fig-colonCVTrsf-1).

```
cv.glmnet(x = t(asinh(Biobase::exprs(zellerNC))),
          y = factor(zellerNC$disease),
          family = "binomial") |> plot()
```

[![](12-chap_files/figure-html/fig-colonCVTrsf-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-colonCVTrsf-1-1.png "Figure 12.22: like Figure fig-colonCV-1, but using an \text{asinh} transformation of the data.")

Figure 12.22: like [Figure 12.21](https://www.huber.embl.de/msmb/12-chap.html#fig-colonCV-1), but using an \(\text{asinh}\) transformation of the data.

Question 12.21

Would a good classification performance on these data mean that this assay is ready for screening and early cancer detection?

Solution

No. The performance here is measured on a set of samples in which the cases have similar prevalence as the controls. This serves well enough to explore the biology. However, in a real-life application, the cases will be much less frequent. To be practically useful, the assay must have a much higher specificity, i.e., rarely diagnose disease where there is none. To establish specificity, a much larger set of normal samples need to be tested.

### 12.6.3 Example: classifying mouse cells from their expression profiles

Figures [12.21](https://www.huber.embl.de/msmb/12-chap.html#fig-colonCV-1) and [12.22](https://www.huber.embl.de/msmb/12-chap.html#fig-colonCVTrsf-1) are textbook examples of how we expect the dependence of (cross-validated) classification performance versus model complexity (\(\lambda\)) to look. Now let’s get back to the mouse embryo cells data. We’ll try to classify the cells from embryonic day `E3.25` with respect to their genotype.

```
sx = x[, x$Embryonic.day == "E3.25"]
embryoCellsClassifier = cv.glmnet(t(Biobase::exprs(sx)), sx$genotype,
                family = "binomial", type.measure = "class")
plot(embryoCellsClassifier)
```

[![](12-chap_files/figure-html/fig-mousecvglmfit-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-mousecvglmfit-1-1.png "Figure 12.23: Cross-validated misclassification error versus penalty parameter for the mouse cells data.")

Figure 12.23: Cross-validated misclassification error versus penalty parameter for the mouse cells data.

In [Figure 12.23](https://www.huber.embl.de/msmb/12-chap.html#fig-mousecvglmfit-1) we see that the misclassification error is (essentially) monotonously increasing with \(\lambda\), and is smallest for \(\lambda\to 0\), i.e., if we apply no penalization at all.

Question 12.22

What is going on with these data?

Solution

It looks that inclusion of more, and even of all features, does not harm the classification performance. In a way, these data are “too easy”. Let’s do a \(t\)-test for all features:

```
mouse_de = rowttests(sx, "genotype")
ggplot(mouse_de, aes(x = p.value)) +
  geom_histogram(boundary = 0, breaks = seq(0, 1, by = 0.01))
```

[![](12-chap_files/figure-html/fig-mousecellsrowttst-1-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-mousecellsrowttst-1-1.png "Figure 12.24: Histogram of p-values for the per-feature t-tests between genotypes in the E3.25 cells.")

Figure 12.24: Histogram of p-values for the per-feature \(t\)-tests between genotypes in the E3.25 cells.

The result, shown in [Figure 12.24](https://www.huber.embl.de/msmb/12-chap.html#fig-mousecellsrowttst-1), shows that large number of genes are differentially expressed, and thus informative for the class distinction. We can also compute the pairwise distances between all cells, using all features.

```
dists = as.matrix(dist(scale(t(Biobase::exprs(x)))))
diag(dists) = +Inf
```

and then for each cell determine the class of its nearest neighbor

```
nn = sapply(seq_len(ncol(dists)), function(i) which.min(dists[, i]))
table(x$sampleGroup, x$sampleGroup[nn]) |> `colnames<-`(NULL)
```

```
                  [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
  E3.25             33    0    0    0    3    0    0    0
  E3.25 (FGF4-KO)    1   15    0    1    0    0    0    0
  E3.5 (EPI)         2    0    3    0    6    0    0    0
  E3.5 (FGF4-KO)     0    0    0    8    0    0    0    0
  E3.5 (PE)          0    0    0    0   11    0    0    0
  E4.5 (EPI)         0    0    0    0    2    2    0    0
  E4.5 (FGF4-KO)     1    0    0    0    0    0    9    0
  E4.5 (PE)          0    0    0    0    2    0    0    2
```

Using all features, the 1 nearest-neighbor classifier is correct in almost all cases, including for the E3.25 wildtype vs FGF4-KO distinction. This means that for these data, there is no apparent benefit in regularization or feature selection. Limitations of using all features might become apparent with truly new data, but that is out of reach for cross-validation.

## 12.7 A large choice of methods

We have now seen three classification methods: linear discriminant analysis (`lda`), quadratic discriminant analysis (`qda`) and logistic regression using elastic net penalization (`glmnet`). In fact, there are hundreds of different learning algorithms[9](https://www.huber.embl.de/msmb/12-chap.html#fn9) available in R and its add-on packages. You can get an overview in the CRAN task view [Machine Learning & Statistical Learning](https://cran.r-project.org/web/views/MachineLearning.html). Some examples are:

9 For an introduction to the subject that uses R and provides many examples and exercises, we recommend ([James et al. 2013](https://www.huber.embl.de/msmb/16-chap.html#ref-James:2013)).

* Support vector machines: the function `svm` in the package **[e1071](https://cran.r-project.org/web/packages/e1071/)**; `ksvm` in **[kernlab](https://cran.r-project.org/web/packages/kernlab/)**
* Tree based methods in the packages **[rpart](https://cran.r-project.org/web/packages/rpart/)**, **[tree](https://cran.r-project.org/web/packages/tree/)**, **[randomForest](https://cran.r-project.org/web/packages/randomForest/)**
* Boosting methods: the functions `glmboost` and `gamboost` in package **[mboost](https://cran.r-project.org/web/packages/mboost/)**
* `PenalizedLDA` in the package **[PenalizedLDA](https://cran.r-project.org/web/packages/PenalizedLDA/)**, `dudi.discr` and `dist.pcaiv` in **[ade4](https://cran.r-project.org/web/packages/ade4/)**).

The complexity and heterogeneity of choices of learning strategies, tuning parameters and evaluation criteria in each of these packages can be confusing. You will already have noted differences in the interfaces of the `lda`, `qda` and `glmnet` functions, i.e., in how they expect their input data to presented and what they return. There is even greater diversity across all the other packages and functions. At the same time, there are common tasks such as cross-validation, parameter tuning and performance assessment that are more or less the same no matter what specific method is used. As you have seen, e.g., in our `estimate_mcl_loocv` function, the looping and data shuffling involved led to rather verbose code.

So what to do if you want to try out and explore different learning algorithms? Fortunately, there are several projects that provide unified interfaces to the large number of different machine learning interfaces in R, and also try to provide “best practice” implementations of the common tasks such as parameter tuning and performance assessment. The two most well-known ones are the packages **[caret](https://cran.r-project.org/web/packages/caret/)** and **[mlr](https://cran.r-project.org/web/packages/mlr/)**. Here were have a look at **[caret](https://cran.r-project.org/web/packages/caret/)**. You can get a list of supported methods through its `getModelInfo` function. There are quite a few, here we just show the first 8.

```
library("caret")
caretMethods = names(getModelInfo())
head(caretMethods, 8)
```

```
[1] "ada"         "AdaBag"      "AdaBoost.M1" "adaboost"    "amdai"      
[6] "ANFIS"       "avNNet"      "awnb"
```

```
length(caretMethods)
```

```
[1] 239
```

We will check out a neural network method, the `nnet` function from the eponymous package. The `parameter` slot informs us on the the available tuning parameters[10](https://www.huber.embl.de/msmb/12-chap.html#fn10).

10 They are described in the manual of the `nnet` function.

```
getModelInfo("nnet", regex = FALSE)[[1]]$parameters
```

```
  parameter   class         label
1      size numeric #Hidden Units
2     decay numeric  Weight Decay
```

Let’s try it out.

```
trnCtrl = trainControl(
  method = "repeatedcv",
  repeats = 3,
  classProbs = TRUE)
tuneGrid = expand.grid(
  size = c(2, 4, 8),
  decay = c(0, 1e-2, 1e-1))
nnfit = train(
  Embryonic.day ~ Fn1 + Timd2 + Gata4 + Sox7,
  data = embryoCells,
  method = "nnet",
  tuneGrid  = tuneGrid,
  trControl = trnCtrl,
  metric = "Accuracy")
```

That’s quite a mouthful, but the nice thing is that this syntax is standardized and applies across many different methods. All you need to do specify the name of the method and the grid of tuning parameters that should be explored via the `tuneGrid` argument.

Now we can have a look at the output ([Figure 12.25](https://www.huber.embl.de/msmb/12-chap.html#fig-ML-nnfit)).

```
nnfit
```

```
Neural Network 

66 samples
 4 predictor
 3 classes: 'E3.25', 'E3.5', 'E4.5' 

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 3 times) 
Summary of sample sizes: 59, 59, 60, 60, 60, 60, ... 
Resampling results across tuning parameters:

  size  decay  Accuracy   Kappa    
  2     0.00   0.7153571  0.4763093
  2     0.01   0.7808333  0.6185837
  2     0.10   0.7623810  0.5804533
  4     0.00   0.7736905  0.6059126
  4     0.01   0.7877778  0.6282527
  4     0.10   0.7601984  0.5719479
  8     0.00   0.7469048  0.5700943
  8     0.01   0.7701190  0.5954734
  8     0.10   0.7457143  0.5471589

Accuracy was used to select the optimal model using the largest value.
The final values used for the model were size = 4 and decay = 0.01.
```

```
plot(nnfit)
predict(nnfit) |> head(10)
```

```
 [1] E3.25 E3.25 E3.25 E3.25 E3.25 E3.25 E3.25 E3.25 E3.25 E3.25
Levels: E3.25 E3.5 E4.5
```

[![](12-chap_files/figure-html/fig-ML-nnfit-1.png)](https://www.huber.embl.de/msmb/12-chap_files/figure-html/fig-ML-nnfit-1.png "Figure 12.25: Parameter tuning of the neural net by cross-validation.")

Figure 12.25: Parameter tuning of the neural net by cross-validation.

Question 12.23

Will the accuracy that we obtained above for the optimal tuning parameters generalize to a new dataset? What could you do to address that?

Solution

No, it is likely to be too optimistic, as we have picked the optimum. To get a somewhat more realistic estimate of prediction performance when generalized, we could formalize (into computer code) all our data preprocessing choices and the above parameter tuning procedure, and embed this in another, outer cross-validation loop ([Ambroise and McLachlan 2002](https://www.huber.embl.de/msmb/16-chap.html#ref-ambroise2002selection)). However, this is likely still not enough, as we discuss in the next section.

### 12.7.1 Method hacking

In [Chapter 6](https://www.huber.embl.de/msmb/06-chap.html) we encountered *p-value hacking*. A similar phenomenon exists in statistical learning: given a dataset, we explore various different methods of preprocessing (such as normalization, outlier detection, transformation, feature selection), try out different machine learning algorithms and tune their parameters until we are content with the result. The measured accuracy is likely to be too optimistic, i.e., will not generalize to a new dataset. Embedding as many of our methodical choices into a computational formalism and having an outer cross-validation loop (not to be confused with the inner loop that does the parameter tuning) will ameliorate the problem. But is unlikely to address it completely, since not all our choices can be formalized.

The gold standard remains validation on truly unseen data. In addition, it is never a bad thing if the classifier is not a black box but can be interpreted in terms of domain knowledge. Finally, report not just summary statistics, such as misclassification rates, but lay open the complete computational workflow, so that anyone (including your future self) can convince themselves of the robustness of the result or of the influence of the preprocessing, model selection and tuning choices ([Holmes 2018](https://www.huber.embl.de/msmb/16-chap.html#ref-Holmes2017)).

## 12.8 Summary of this chapter

We have seen examples of machine learning applications; we have focused on predicting categorical variables (like diabetes type or cell class). Predicting continuous outcomes is also part of machine learning, although we have not considered it here. There are many parallels and overlaps between *machine learning* and *statistical regression* (which we studied in [Chapter 8](https://www.huber.embl.de/msmb/08-chap.html)). One can consider them two different names for pretty much the same activity, although each has its own flavors: in machine learning, the emphasis is on the prediction of the outcome variables, whereas in regression we often care at least as much about the role of the covariates – which of them have an effect on the outcome, and what is the nature of these effects? In other words, we do not only want predictions, we also want to understand them.

We saw linear and quadratic discriminant analysis, two intuitive methods for partitioning a two-dimensional data plane (or a \(p\)-dimensional space) into regions using either linear or quadratic separation lines (or hypersurfaces). We also saw logistic regression, which takes a slightly different approach but is more amenable to operating in higher dimensions and to regularization.

We encountered the main challenge of machine learning: how to avoid overfitting? We explored why overfitting happens in the context of the so-called curse of dimensionality, and we learned how it may be overcome using regularization.

In other words, machine learning would be easy if we had infinite amounts of data representatively covering the whole space of possible inputs and outputs[11](https://www.huber.embl.de/msmb/12-chap.html#fn11). The challenge is to make the best out of a finite amount of training data, and to generalize these to new, unseen inputs. There is a vigorous trade-off between the amount, resolution and coverage of training data and the complexity of the model. Many models have continuous parameters that enable us to “tune” their complexity or the strength of their regularization. Cross-validation can help us with such tuning, although it is not a panacea, and caveats apply, as we saw in [Section 12.6.3](https://www.huber.embl.de/msmb/12-chap.html#sec-supervised-classifymousecells).

11 It would “just” be a formidable database / data management problem.

## 12.9 Further reading

* An introduction to statistical learning that employs many concrete data examples and uses little mathematical formalism is given by James et al. ([2013](https://www.huber.embl.de/msmb/16-chap.html#ref-James:2013)). An extension, with more mathematical background, is the textbook by Hastie, Tibshirani, and Friedman ([2008](https://www.huber.embl.de/msmb/16-chap.html#ref-HastieTibshiraniFriedman)).
* The [CRAN task view on machine learning](https://cran.r-project.org/web/views/MachineLearning.html) gives an overview over machine learning software in R.
* [RStudio’s API for the “deep learning” platforms Keras and TensorFlow](https://tensorflow.rstudio.com/keras) and the associated teaching materials and demos are a good place to try out some of the recent developments in this field.