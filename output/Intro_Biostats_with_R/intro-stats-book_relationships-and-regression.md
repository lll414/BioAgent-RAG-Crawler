Source URL: https://tuos-bio-data-skills.github.io/intro-stats-book/relationships-and-regression.html
Date Scraped: 2025-12-12

---

# Chapter 12 Relationships and regression

## 12.1 Introduction

Much of biology is concerned with relationships between numeric variables. For example:

* We sample fish and measure their length and weight because we want to understand how weight changes with length.
* We survey grassland plots and measure soil pH and species diversity to understand how species diversity depends on soil pH.
* We manipulate temperature and measure fitness in insects because we want to describe their thermal tolerance.

In the previous chapter, we learnt about one technique for analysing associations between numeric variables (correlation). A correlation coefficient quantifies the strength and direction of an association between two variables: it will be close to zero if there is no association between the variables, whereas a strong association is implied if the coefficient is near -1 or +1. A correlation coefficient tells us nothing about the form of a relationship. Nor does it allow us to predict the value of one variable from the value of a second variable.

In contrast to correlation, a regression analysis allows us to make precise statements about how one numeric variable depends on the values of another. Graphically, we can evaluate such dependencies using a scatter plot. We may be interested in knowing:

1. Are the variables related or not? There’s not much point in studying a relationship that isn’t there:

![](intro-bio-stats-book_files/figure-html/reg-eg-related-1.png)

2. Is the relationship positive or negative? Sometimes we can answer a scientific question just by knowing the direction of a relationship:

![](intro-bio-stats-book_files/figure-html/reg-eg-posneg-1.png)

3. Is the relationship a straight line or a curve? It is important to know the form of a relationship if we want to make predictions:

![](intro-bio-stats-book_files/figure-html/reg-eg-linornot-1.png)

Although sometimes it may be obvious that there is a relationship between two variables from a plot of one against the other, at other times it may not. Take a look at the following:

![](intro-bio-stats-book_files/figure-html/reg-eg-confidence-1.png)

We might not be very confident in judging which, if either, of these plots provides evidence of a positive relationship between the two variables. Maybe the pattern we perceive can just be explained by sampling variation, or perhaps it can’t. We need a procedure—a statistical test—to evaluate how likely the relationship could have arisen as a result of sampling variation. In addition to judging the statistical significance of a relationship, we may also be interested in describing the relationship mathematically – i.e. finding the equation of the best fitting line through the data.

A linear regression analysis allows us to do all this.

## 12.2 What does linear regression do?

**Simple linear regression** allows us to *predict* how one variable (the **response variable**) *responds* or *depends* on to another (the **predictor variable**), assuming a straight-line relationship.

* What does the word ‘simple’ mean here? Simple linear regression is a regression model which only accounts for one predictor variable. If more than one predictor variable is considered, the correct term to describe the resulting model is ‘multiple regression’. Multiple regression is a very useful tool, but we’re only going to study simple regression in this book.
* What does the word ‘linear’ mean here? In statistics, the word linear is used in two slightly different but closely related ways. When discussing simple linear regression, the term linear is often understood to mean that the relationship follows a straight line. That’s all. The more technical definition concerns the relationship between the parameters of a statistical model. We don’t need to worry about that one here.

Writing ‘simple linear regression’ all the time becomes tedious, so we’ll often write ‘linear regression’ or ‘regression’. Just keep in mind that we’re always referring to simple linear regression in this book. These regression models account for a straight-line relationship between two numeric variables, i.e. they describe how the response variable changes in response to the values of the predictor variable. It is conventional to label the response variable as ‘yy’ and the predictor variable as ‘xx’. When we present such data graphically, the response variable always goes on the yy-axis and the predictor variable on the xx-axis. Try not to forget this convention!

#### ‘response / predictor’ or ‘dependent / independent’?

Another way to describe linear regression is that it finds the straight-line relationship which best describes the dependence of one variable (the **dependent variable**) on the other (the **independent variable**). The dependent vs. independent and response vs. predictor conventions for variables in a regression are equivalent. They only differ in the words they use to describe the variables involved.

To avoid confusion, we will stick with response / predictor naming convention in this book.

How do we decide how to select which is to be used as the response variable and which as the predictor variable? The decision is straightforward in an experimental setting: the manipulated variable is the predictor variable, and the measured outcome is the response variable. Consider the thermal tolerance example from earlier. The temperature was manipulated in that experiment, so it must be designated the predictor variable. Moreover, *a priori* (before conducting the experiment), we can reasonably suppose that temperature changes may cause changes in enzyme activity, but the reverse seems pretty unlikely.

Things are not so clear cut when working with data from an observational study. Indeed, the ‘response vs predictor’ naming convention can lead to confusion because the term ‘response variable’ tends to make people think in terms of causal relationships, i.e. that variation in the the predictor somehow *causes* changes in the response variable. Sometimes that is true, such as in the experimental setting described above. However, deciding to call one variable a response and the other a predictor should not be taken to automatically imply a casual relationship between them.

Finally, it matters which way round we designate the response and predictor variable in a regression analysis. Suppose we have two variables A and B. In that case, the relationship we find from a regression will not be the same for A against B as for B against A. Choosing which variable to designate as the response boils down to working out which of the variables needs to be explained (response) in terms of the other (predictor).

## 12.3 How does simple linear regression work?

### 12.3.1 Finding the best fit line

If we draw a straight line through a set of points on a graph then, unless they form a perfect straight line, some points will lie close to the line and others further away. The vertical distances between the line and each point (i.e. measured parallel to the yy-axis) have a special name. They are called the *residuals*. Here’s a visual example:

![Example of data (blue points) used in a simple regression. A fitted line and the associated residuals (vertical lines) are also shown](intro-bio-stats-book_files/figure-html/reg-eg-with-resids1-1.png)

Figure 12.1: Example of data (blue points) used in a simple regression. A fitted line and the associated residuals (vertical lines) are also shown

In this plot the blue points are the data and the vertical lines represent the residuals. The residuals represent the ‘left over’ variation after the line has been fitted through the data. They indicate how well the line fits the data. If all the points lay close to the line, the variability of the residuals would be low relative to the overall variation in the response variable, yy. When the observations are more scattered around the line, the variability of the residuals would be large relative to the variation in the response variable, yy.

Regression works by finding the line which minimises the size of the residuals in some sense. We’ll explain exactly how in a moment. The following illustration indicates the principle of this process:

![](intro-bio-stats-book_files/figure-html/reg-eg-four-plots1-1.png)

The data are identical in all four graphs, but in the top left-hand graph a horizontal line (i.e. no effect of xx on yy) has been fitted, while on the remaining three graphs sloping lines of different magnitude have been fitted.

#### Which line is best?

One of the four lines is the ‘line of best’ fit from a regression analysis. Spend a few moments looking at the four figures. Which line seems to fit the data best? Why do you think this line is ‘best’?

Let’s visualise the data, the candidate lines and the residuals:

![](intro-bio-stats-book_files/figure-html/reg-eg-four-plots2-1.png)

We said that regression works by finding the intercept and slope that minimises the vertical distances between the line and each observation in some way[7](#fn7). In fact, it minimises something called the ‘sum of squares’ of these distances: we calculate a sum of squares for a particular set of observations and a fitted line by squaring the residual distances and adding all of these up. This quantity is called the **residual sum of squares**. The line with the *lowest* residual sum of squares is the best line because it ‘explains’ the most variation in the response variable.

You should be able to see that, for the horizontal line (‘A’), the residual sum of squares is larger than any of the other three plots with the sloping lines. This suggests that the sloping lines fit the data better. Which one is best among the three we’ve plotted? To get at this, we need to calculate the residual sum of squares for each line:

```
##   Line    Residual Sum of Squares
## 1    A                   17.55067
## 2    B                   11.97966
## 3    C                   10.12265
## 4    D                   12.79674
```

So it looks like the line in panel C is the best fitting line among the candidates. In fact, it is the best fit line among all possible candidates. Did you manage to guess this by looking at the lines and the raw data? If not, think about why you got the answer wrong. Did you consider the vertical distances or the perpendicular distances?

#### Know your residuals

It is important to understand what a residual represents. Why? Because they pop up all the time when working with statistical models (not only regression, in fact). It is hard to understand what R is telling you about a model without knowing about these residual things.

## 12.4 What do you get out of a regression?

A regression analysis involves two activities:

* **Interpretation.** When we ‘fit’ a regression model to data we estimate the coefficients of a best-fit straight line through the data. This is the equation that best describes how the yy (response) variable *responds to* the xx (predictor) variable. To put it in slightly more technical terms, it describes the yy variable as a function of the xx variable. This model may be used to understand how the variables are related or make predictions.
* **Inference.** It is not enough to just estimate the regression equation. Before we can use it, we need to determine whether there is a statistically significant relationship between the xx and yy variables. That is, the analysis will tell us whether an apparent association is likely to be real or just a chance outcome resulting from sampling variation.

### 12.4.1 Interpreting a regression

What is the form of the relationship? The equation for a straight-line relationship is y=a+b×xy=a+b×x, where

* yy is the response variable,
* xx is the predictor variable,
* aa is the intercept (i.e. where the line crosses the yy axis), and
* bb is the slope of the line.

The aa and the bb are referred to as the *coefficients* (or *parameters*) of the line. The slope of the line is often the coefficient we care about most. It tells us the amount by which yy changes for a change of one unit in xx. If the value of bb is positive (i.e. a plus sign in the above equation) this means the line slopes upwards to the right. A negative slope (y=a−bxy=a−bx) means the line slopes downwards to the right. The diagram below shows the derivation of an equation for a straight line.

![](intro-bio-stats-book_files/figure-html/reg-line-explain-1.png)

Having the equation for a relationship allows us to predict the value of the yy variable for any value of xx. For example, in the thermal tolerance example, we want an equation that will allow us to work out how fitness changes with temperature. Such predictions can be made by hand (see below) or using R (details later).

In the above diagram, the regression equation is: y=1+0.66xy=1+0.66x. So to find the value of yy at x=2x=2 we use: y=1+(0.667×2)=2.32y=1+(0.667×2)=2.32. Obviously, by finding yy values for 2 (or preferably 3) different xx values from the equation, the actual line can easily be plotted on a graph manually if required—plot the values and join the dots! It’s much easier to use R to do this kind of thing though.

#### Regression involves a statistical model

A simple linear regression is underpinned by a statistical model. If you skim back through the [parametric statistics](parametric-statistics.html#parametric-statistics) chapter, you will see that the equation y=a+b×xy=a+b×x represents the ‘systematic component’ of the regression model. This bit describes the component of variation in yy that is explained by the model for the dependence of yy on xx. The residuals correspond to the ‘random component’ of the model. These represent the component of variation in the yy variable that our regression model fails to describe.

### 12.4.2 Evaluating hypotheses (‘inference’)

More than one kind of significance test can be carried out with a simple linear regression. We’re going to focus on the most common test: an *F* test of whether the slope coefficient is significantly different from 0. This addresses the important question, “Is there a relationship?”

How do we do this? We play exactly the same kind of gambit we used to develop the earlier significance tests:

1. We start with a null hypothesis of ‘no effect’. This corresponds to the hypothesis that the slope of the regression is zero.
2. We then work out what the distribution of some kind of test statistic should look like under the null hypothesis. The test statistic in this case is called the *F*-ratio.
3. We then calculate a *p*-value by asking how likely it is that we would see the observed test statistic, or a more extreme value, if the null hypothesis were really true.

It’s not so critical that someone understands the mechanics of an *F*-test to use it. However, knowing a bit about where is comes from does help to demystify the output produced by R. To that end, let’s step through the calculations involved in the *F* test using the example data shown in the above four-panel plot.

#### Total variation

We first need to calculate something called the **total sum of squares**. The figure below shows the raw data (blue points) and the grand mean (i.e. the sample mean).

![](intro-bio-stats-book_files/figure-html/reg-eg-total-1.png)

The vertical lines show the distance between each observation and the grand mean. These vertical lines are just the residuals from a model where the slope of the line is set to zero. What we need to do is quantify the variability of these residuals. We can’t just add them up, because by definition, they have to sum to zero, i.e. they are calculated relative to the grand mean.

Instead, we calculate the total sum of squares by taking each residual in turn, squaring it, and then adding up all the squared values. We call this the total sum of squares because it is a measure of the total variability in the response variable, yy. This number is 17.55 for the data in the figure above.

#### Residual variation

Next we need to calculate the **residual sum of squares**. We have already seen how this calculation works because it is used in the calculation of the best fit line—the best fit line is the one that minimises the residual sum of squares. Let’s plot this line along with the associated residuals of this line again:

![](intro-bio-stats-book_files/figure-html/reg-eg-with-resids2-1.png)

The vertical lines show the distance between each observation and the best fit line. We need to quantify the variability of these residuals. Again, we can’t just add up the deviations because they have to sum to zero as a result of how the best fit line is found. Instead we calculate the residual sum of squares by taking each residual in turn, squaring it, and then adding up all the squared values. This number is 10.12 for the figure above. We call this the residual, or error, sum of squares because it is a measure of the variation in yy that is ‘left over’ after accounting for the influence of the predictor variable xx.

#### Explained variation

Once the total sum of squares and the residual sum of squares are known, we can calculate the quantity we really want: the **explained sum of squares**. This is a measure of the variation in yy that is explained by the influence of the predictor variable xx. We calculate this by subtracting the residual sum of squares from the total sum of squares. This makes intuitive sense: if we subtract the variation in yy we can’t explain (residual) from all the variation in yy (total), we end up with the amount ‘explained’ by the regression. This number is 7.43 for the example.

#### Degrees of freedom, mean squares and *F* tests

The problem with sums of squares is that they are a function of sample size. The more data we have, the larger our sum of squares will get. The solution to this problem is to convert them into a measure of variability that doesn’t scale with sample size. We need to calculate **degrees of freedom** (written as df, or d.f.) to do this. We came across the concept of degrees of freedom when we studied the *t*-test. The idea is closely related to sample size. It is difficult to give a precise definition, but roughly speaking, the degrees of freedom of a statistic is a measures of how much information it is based on (bigger is better).

Each of the measures of variability we just calculated for the simple linear regression has a degrees of freedom associated with it. We need the explained and error degrees of freedom:

* Explained d.f. = 1
* Error d.f. = (Number of observations - 2)

Don’t worry if those seem a little cryptic. We don’t need to carry out degrees of freedom calculations by hand because R will do them for us. We’ll think about degrees of freedom a bit more when we start to learn about ANOVA models.

The reason degrees of freedom matter is because we can use them to standardise the sum of squares to account for sample size. The calculations are very simple. We take each sum of squares and divide it by its associated degrees of freedom. The resulting quantity is called a **mean square** (it’s the mean of squared deviations): Mean Square=Sum of SquaresDegrees of FreedomMean Square=Sum of SquaresDegrees of Freedom A mean square is actually an estimate of variance. Remember the variance? It is one of the standard measures of a distribution’s dispersion, or spread.

Now for the important bit. The two mean squares can be compared by calculating the ratio between them, which is designated by the letter *F*:

F=Variance Ratio=Explained Mean SquareResidual Mean SquareF=Variance Ratio=Explained Mean SquareResidual Mean Square

This is called the *F* ratio, or sometimes, the variance ratio. If the explained variation is large compared to the residual variation, then the *F* ratio will be large. Conversely, if the explained variation is relatively small the *F* ratio will be small. We can see where this is heading…

The *F* ratio is a type of test statistic—if the value of *F* is sufficiently large, we judge it to be statistically significant. For this judgement to be valid, we have to make one key assumption about the population from which the data has been sampled: we assume the residuals are drawn from a normal distribution. If this assumption is correct, it can be shown that the distribution of the *F* ratio under the null hypothesis (the ‘null distribution’) has a particular form: it follows a theoretical distribution called an *F*-distribution. And yes, that’s why the variance ratio is called ‘*F*’.

All this means we can assess statistical significance of the slope coefficient by comparing the *F* ratio calculated from a sample to this theoretical distribution. This procedure is called an *F* test. The *F* ratio is 7.34 in our example. This is quite high, which indicates that the slope is likely to be significantly different from 0.

We need more than an *F* ratio to calculate a *p*-value. We also need to consider the degrees of freedom of the test. Furthermore, because this test involves an *F* ratio, **there are two different degrees of freedom to consider**: the explained and residual df’s. Try to remember that fact—*F* ratio tests involve a pair of degrees of freedom.

We could calculate a *p*-value now by messing around with tables of *F* ratios and degrees of freedom. However, it is much easier to let R do this for us when working with a regression model. That’s what we’re going to do in the next chapter, so for now, we’ll leave significance tests alone.

## 12.5 Correlation or regression?

Before we finish up, it is worth pausing to review the difference between regression and correlation analyses. Whilst regression and correlation are both concerned with associations between numeric variables, they are different techniques and each is appropriate under distinct circumstances. This is a frequent source of confusion. Which technique is required for a particular analysis depends on

* the way the data were collected and
* the goal of the analysis.

There are two broad questions to consider—

#### Where do the data come from?

Think about how the data have been collected. If they are from a study where one of the variables has been experimentally manipulated, then choosing the best analysis is easy. We should use a regression analysis, in which the predictor variable is the experimentally manipulated variable and the response variable is the measured outcome. The fitted line from a regression analysis describes how the outcome variable *depends on* the manipulated variable—it describes the causal relationship between them.

It is generally inappropriate to use correlation to analyse data from an experimental setting. A correlation analysis examines association but does not imply the dependence of one variable on another. Since there is no distinction of response or predictor variables, it doesn’t matter which way round we do a correlation. In fact, the phrase ‘which way round’ doesn’t even make sense in the context of a correlation.

If the data are from a so-called ‘observational study’ where someone took measurements, but nothing was actually manipulated experimentally, then either method may be appropriate.

#### What is the goal of the analysis?

Think about what question is being asked. A correlation coefficient only quantifies the strength and direction of an association between variables—it tells us nothing about the form of a relationship. Nor does it allow us to make predictions about the value of one variable from the values of a second. A regression does allow this because it involves fitting a line through the data—i.e. regression involves **a model** for the relationship.

This means that if the goal is to understand the form of a relationship between two variables, or to use a fitted model to make predictions, we have to use regression. If we just want to know whether two variables are associated or not, the direction of the association, and whether the association is strong or weak, then a correlation analysis is sufficient. It is better to use a correlation analysis when the extra information produced by a regression is not needed, because the former will be simpler and potentially more robust.

---

7. Notice that it is the vertical distance that matters, not the perpendicular distance from the line.[↩︎](relationships-and-regression.html#fnref7)