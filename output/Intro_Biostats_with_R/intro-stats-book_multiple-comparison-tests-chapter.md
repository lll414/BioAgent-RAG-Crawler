Source URL: https://tuos-bio-data-skills.github.io/intro-stats-book/multiple-comparison-tests-chapter.html
Date Scraped: 2025-12-12

---

# Chapter 16 Multiple comparison tests

In the [One-way ANOVA in R](one-way-anova-in-R.html#one-way-anova-in-R) chapter, we learned how to examine the global hypothesis of no difference between means. That test does not evaluate *which* means might be driving a significant result. For example, in the corncrake example, we found evidence of a significant effect of dietary supplement on the mean hatchling growth rate. However, our analysis did not tell us which supplements are better or worse and did not make any statements about which supplements differ significantly from each other. The purpose of this chapter is to examine one method for assessing where the differences lie.

## 16.1 Post hoc multiple comparisons tests

The general method we will use is called a **post hoc multiple comparisons test**. The phrase ‘post hoc’ refers to the fact that these tests are conducted without any particular prior comparisons in mind. The words ‘multiple comparisons’ refer to the fact that they consider many different pairwise comparisons. There are quite a few multiple comparison tests—Scheffé’s test, the Student-Newman-Keuls test, Duncan’s new multiple range test, Dunnett’s test, … (the list goes on and on). Each one applies to particular circumstances, and none is universally accepted as ‘the best’. We are going to work with the most widely used test: the **Tukey multiple comparison test**. This test is also known as Tukey’s Honestly Significant Difference (Tukey HSD) test[11](#fn11).

People tend to favour Tukey’s HSD test because it is ‘conservative’: the test has a low false-positive rate compared to the alternatives. A false positive occurs when a test turns up a statistically significant result for an effect that is not there. A low false-positive rate is a good thing. It means that if we find a significant difference, we can be more confident it is ‘real’. The cost of using the Tukey HSD test is that it isn’t as powerful as alternatives: the test turns up a lot of false negatives. A false negative occurs when a test fails to produce a statistically significant result for an effect when it is present. A test with a high false-negative rate tends to miss effects.

One line of thinking says post hoc multiple comparisons tests of any kind should never be undertaken. We shouldn’t carry out an experiment without a prior prediction of what will happen—we should know which comparisons need to be made and should only undertake those particular comparisons rather than making every possible comparison. Nonetheless, post hoc multiple comparisons tests are easy to apply and widely used, so there is value in knowing how to use them. The Tukey HSD test at least tends to guard against picking up non-existent effects.

## 16.2 Tukey’s HSD test in R

We will use the corncrake example to illustrate how to use R to carry out and interpret Tukey’s HSD test.

#### 

The code below assumes CORN\_CRAKE.CSV has been read into a tibble called `corn_crake`. It also assumes a one-way ANOVA model for the effect of diet supplement on hatchling growth rate has already been fitted and given the name `corncrake_model`:

```
corncrake_model <- lm(WeightGain ~ Supplement, data = corn_crake)
```

In the [One-way ANOVA in R](one-way-anova-in-R.html#one-way-anova-in-R) chapter, we saw how to use `anova` on a model object to evaluate whether means differ, for example:

```
anova(corncrake_model)
```

```
## Analysis of Variance Table
## 
## Response: WeightGain
##            Df Sum Sq Mean Sq F value   Pr(>F)   
## Supplement  4 357.65  89.413  5.1281 0.002331 **
## Residuals  35 610.25  17.436                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

That gives us the result of the global significance test. Unfortunately, we need different functions to carry out the Tukey HSD test. First, we convert the linear model object into a different kind of model object using the `aov` function:

```
corncrake_aov <- aov(corncrake_model)
```

We don’t need to understand too much about what this is doing. In a nutshell, `aov` does some extra calculations needed to perform a Tukey HSD test. We gave the new object a name (`corncrake_aov`) because we need to use it in the next step.

It’s easy to perform a Tukey HSD test once we have the ‘aov’ version of our model. There are a few different options. Here is how to do this using the `TukeyHSD` function:

```
TukeyHSD(corncrake_aov, ordered = TRUE)
```

Pay attention! We applied the `TukeyHSD` function to the ‘aov’ object, *not* the original `lm` object. We have suppressed the output for now. Before we review it we need to get an idea of what it is going to show us.

The `ordered = TRUE` tells `TukeyHSD` that we want to order the treatment means from smallest to largest and then apply every pairwise comparison, starting with the smallest mean (‘Allvit’) and working up through the order. Here are the means ordered from smallest to largest, working left to right:

|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| Supplement | Allvit | None | Linseed | Sizefast | Earlybird |
| Mean | 14.3 | 15.6 | 18.5 | 19.0 | 22.9 |

So the `TukeyHSD` with `ordered = TRUE` will first compare ‘Allvit’ to ‘None’, then ‘Allvit’ to ‘Linseed’, then ‘Allvit’ to ‘Sizefast’, then ‘Allvit’ to ‘Earlybird’, then ‘None’ to ‘Linseed’, then ‘None’ to ‘Sizefast’, … and so on, until we get to ‘Sizefast’ vs. ‘Earlybird’. Let’s look at the output:

```
##   Tukey multiple comparisons of means
##     95% family-wise confidence level
##     factor levels have been ordered
## 
## Fit: aov(formula = corncrake_model)
## 
## $Supplement
##                     diff       lwr       upr     p adj
## None-Allvit        1.375 -4.627565  7.377565 0.9638453
## Linseed-Allvit     4.250 -1.752565 10.252565 0.2707790
## Sizefast-Allvit    4.750 -1.252565 10.752565 0.1771593
## Earlybird-Allvit   8.625  2.622435 14.627565 0.0018764
## Linseed-None       2.875 -3.127565  8.877565 0.6459410
## Sizefast-None      3.375 -2.627565  9.377565 0.4971994
## Earlybird-None     7.250  1.247435 13.252565 0.0113786
## Sizefast-Linseed   0.500 -5.502565  6.502565 0.9992352
## Earlybird-Linseed  4.375 -1.627565 10.377565 0.2447264
## Earlybird-Sizefast 3.875 -2.127565  9.877565 0.3592201
```

This table look confusing at first glance. It enables you to look up every pair of treatments, and see whether they are significantly different from each other. Lets see how this works… The first four lines compare the Allvit treatment (‘Allvit’) with each of the other treatments in turn:

```
##                     diff       lwr       upr     p adj 
## None-Allvit        1.375 -4.627565  7.377565 0.9638453 
## Linseed-Allvit     4.250 -1.752565 10.252565 0.2707790 
## Sizefast-Allvit    4.750 -1.252565 10.752565 0.1771593 
## Earlybird-Allvit   8.625  2.622435 14.627565 0.0018764
```

So to look up the difference between the control treatment and the ‘Allvit’ treatment, we read the first results row in the table. This says the means differ by 1.375, the confidence interval associated with this difference is [-4.63, 7.38], and that the comparison has a *p*-value of 0.96. So in this case we would conclude that there was a no significant difference between the control treatment and the ‘Allvit’ treatment. We could look up any comparison of the ‘Allvit’ treatment with a different treatment in the next three lines of this portion of the table.

This basic logic extends to the rest of the table. If we want to know whether the ‘Earlybird’ treatment is different from the control, we look up the ‘Earlybird-None’ row:

```
##                     diff       lwr       upr     p adj 
## Earlybird-None     7.250  1.247435 13.252565 0.0113786
```

It looks like the means of the ‘Earlybird’ and ‘None’ levels are significantly different at the *p* < 0.05 level.

Now we know how to look up any set of comparisons, we need to see whether the difference is significant. The next question is: How should we summarise such a table?

### 16.2.1 How to summarise multiple-comparison results

Summarising the results of multiple comparison tests can be a tricky business. The first rule is: don’t present the results like the `TukeyHSD` function does! A clear summary of the results will help us to interpret them correctly and makes it easier to explain them to others. How should we do this? Let’s list the means in order from smallest to largest again:

|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| Supplement | Allvit | None | Linseed | Sizefast | Earlybird |
| Mean | 14.3 | 15.6 | 18.5 | 19.0 | 22.9 |

Now, using the table of Tukey results, we can perform a sequence of pairwise comparisons between the supplements starting with the smallest pair… ‘Allvit’ and ‘None’. The appropriate test is in the first table:

```
##       diff        lwr        upr      p adj 
##  1.3750000 -4.6275646  7.3775646  0.9638453
```

The last column gives the *p*-value, which in this case is certainly not significant (it is much greater than 0.05), so we conclude there is no difference between the ‘Allvit’ and ‘None’ treatments. So now we continue with ‘Allvit’, but compare it to the next larger mean (‘Linseed’). In this case the values are:

```
##      diff       lwr       upr     p adj 
##  4.250000 -1.752565 10.252565  0.270779
```

The last column gives the *p*-value, which again is not significant, so we conclude there is no difference between the ‘Allvit’ and ‘Linseed’ treatments. So now we continue with ‘Allvit’, but compare it to the next larger mean (‘Sizefast’). In this case, the values are:

```
##       diff        lwr        upr      p adj 
##  4.7500000 -1.2525646 10.7525646  0.1771593
```

Once again, this difference is not significant, so we conclude there is no difference between the ‘Allvit’ and ‘Sizefast’ treatments either. So again, we continue with ‘Allvit’, which we compare to the next larger mean (‘Earlybird’).

```
##       diff        lwr        upr      p adj 
##  8.6250000  2.6224354 14.6275646  0.0018764
```

This time the *p*-value is clearly less than 0.05, so we conclude that this pair of treatments are significantly different. We record this by marking ‘Allvit’, ‘None’, ‘Linseed’ and ‘Sizefast’ to indicate that they don’t differ from each other. We’ll use the letter ‘b’ to do this.

|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| Supplement | Allvit | None | Linseed | Sizefast | Earlybird |
| Mean | 14.3 | 15.6 | 18.5 | 19.0 | 22.9 |
|  | b | b | b | b |  |

The ‘b’ defines a group of treatment means—‘Allvit’, ‘None’, ‘Linseed’ and ‘Sizefast’—which are not significantly different from one another. It doesn’t matter which letter we use by the way (the reason for using ‘b’ here will become apparent in a moment).

The means are ordered from smallest to largest, which means we can forget about ‘None’, ‘Linseed’ and ‘Sizefast’ treatments for a moment—if they are not significantly different from ‘Allvit’ they can’t be significantly different from one another.

We move on to ‘Earlybird’, but now, we *work back down* the treatments to see if we can define another overlapping group of means that are not significantly different from one another. When we do this, we find that ‘Earlybird’ is not significantly different from ‘Linseed’ and ‘Sizefast’, but that it is significantly different from ‘None’. This forms a second ‘not significantly different’ group. We will denote this with a new letter (‘a’) in our table:

|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| Supplement | Allvit | None | Linseed | Sizefast | Earlybird |
| Mean | 14.3 | 15.6 | 18.5 | 19.0 | 22.9 |
|  | b | b | b | b |  |
|  |  |  | a | a | a |

If there were additional treatments with a mean that was greater than ‘Earlybird’, we would have to carry on this process, *working back up* from ‘Earlybird’. Thankfully, there are no more treatments, so we are finished.

This leaves us with a concise and complete summary of where the differences between treatments are, which greatly simplifies the task of interpreting the results. Treatment labels that share a letter represent sets of means that are not different from each other. Treatments that are not linked by one or more shared letters are significantly different.

#### Significant ANOVA but no differences in a Tukey test?

ANOVA and the Tukey HSD test are different tests, with different goals. Because of this, it is possible to end up with a significant result from ANOVA, indicating at least one difference between means, but fail to get any differences detected by the Tukey test. This situation can arise when there are many group means being compared and the ANOVA result is marginal (close to *p* = 0.05).

When this happens, there isn’t much we can do except make the best interpretation we can from inspecting the data, and be suitably cautious in the conclusions we draw. It might be tempting to run a new *post hoc* analysis using a different kind of test. Don’t do this. It is a terrible strategy for doing statistics because this kind of practise is guaranteed to increase the chances of finding a false positive.

### 16.2.2 Doing it the easy way…

The results table we produced is concise and complete but no reasonable person would say it was easy to arrive at. Fortunately, someone has written an R function to do this for us. It isn’t part of ‘base R’ though, so we have to install a package to use it. The package we need is called `agricolae`:

```
install.packages("agricolae")
```

Once this has been installed and then loaded, we can use its `HSD.test` function to find the ‘not significantly different’ groups via a Tukey HSD test:

```
HSD.test(corncrake_aov, "Supplement", console=TRUE)
```

```
## 
## Study: corncrake_aov ~ "Supplement"
## 
## HSD Test for WeightGain 
## 
## Mean Square Error:  17.43571 
## 
## Supplement,  means
## 
##           WeightGain      std r       se Min Max   Q25  Q50   Q75
## Allvit        14.250 5.599745 8 1.476301   5  22 10.75 14.5 17.75
## Earlybird     22.875 3.482097 8 1.476301  17  27 20.50 24.5 25.00
## Linseed       18.500 3.023716 8 1.476301  15  24 16.00 18.5 20.00
## None          15.625 3.420004 8 1.476301  10  20 13.75 15.5 17.75
## Sizefast      19.000 4.780914 8 1.476301  10  25 17.25 19.5 22.25
## 
## Alpha: 0.05 ; DF Error: 35 
## Critical Value of Studentized Range: 4.065949 
## 
## Minimun Significant Difference: 6.002565 
## 
## Treatments with the same letter are not significantly different.
## 
##           WeightGain groups
## Earlybird     22.875      a
## Sizefast      19.000     ab
## Linseed       18.500     ab
## None          15.625      b
## Allvit        14.250      b
```

The `console = TRUE` argument tells the function to print the results for us. That’s a lot of output, but we can ignore most of it. The part that matters most is the table at the very end. This shows the group identities as letters, the treatment names, and the treatment means. If we compare that table with the one we just made, we can see they convey the same information. The package labels each group with a letter. For example, we can see that ‘Linseed’ and ‘SizeFast’ are both members of the ‘a’ and ‘b’ group.

So, there is no need to build these Tukey HSD tables by hand. Use the `HSD.test` function in the `agricolae` package instead. Why did we do it the long way first? Well, the usual reasoning applies: it helps us understand how the ‘letters notation’ works.

## 16.3 Summarising and presenting the results of a Tukey test

As with any statistical test it will usually be necessary to summarise the Tukey HSD test in a written form. This gets quite complicated when both the global significance test and the multiple comparisons tests need to be presented. In most cases, it is best to summarise the ANOVA results and concentrate on those comparisons that relate to the original hypothesis we were interested in. Then, refer to a table or figure for the additional detail. For example…

> There was a significant effect of supplement on the weight gain of hatchlings (ANOVA: F=5.1; d.f.= 4,35; p<0.01) (Figure 1). The only supplement that led to a significantly higher rate of weight gain than the control group was Earlybird (Tukey multiple comparisons test, p < 0.05).

When it comes to presenting the results in a report, we really need some way of presenting the means, and the results of the multiple comparison test, as the statement above cannot entirely capture the form of the results. The information can often be conveniently incorporated into a table or figure, using more or less the same format as the output from the `HSD.test` function in the `agricolae` package.

An example table might be:

| **Supplement** | Mean weight gain (g) |
| --- | --- |
| Earlybird | 22.9a |
| Sizefast | 19.0ab |
| Linseed | 18.5ab |
| None | 15.6b |
| Allvit | 14.3b |

However we present it, we need to provide some explanation to set out:

* what test we did,
* what the letter codes mean, and
* the critical threshold we used to judge significance.

In this case the information could be presented in a table legend:

> Table 1: Mean weight gain of hatchlings in the five supplement treatments. Means followed by the same letter did not differ significantly (Tukey test, p>0.05).

Letter coding can also be used effectively in a figure. Again, we must ensure all the relevant information is given in the figure legend.

---

11. N.B. Try to avoid a common mistake: it is *Tukey*, after its originator, the statistician Prof. John Tukey, not Turkey, a large domesticated bird which has made no useful contributions to statistical theory or practice.[↩︎](multiple-comparison-tests-chapter.html#fnref11)