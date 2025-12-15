Source URL: https://tuos-bio-data-skills.github.io/intro-stats-book/goodness-of-fit-tests.html
Date Scraped: 2025-12-12

---

# Chapter 26 Goodness of fit tests

## 26.1 When do we use a chi-square goodness of fit test?

A χ2 goodness of fit test is appropriate in situations where we are studying a categorical variable, and we want to compare the frequencies of each category to pre-specified, expected values. Here are a couple of examples:

* We’ve already seen one situation where a goodness of fit test might be useful: the analysis of the sex ratio among biology undergraduates. We have a prior prediction about what the sex ratio should be in the absence of bias (1:1) and we wanted to know if there was any evidence for sex-related bias in the decision to study biology. We can use a goodness of fit test to compare the number of males and females in a sample of students with the predicted values to determine whether the data are consistent with the equal sex ratio prediction.
* Red campion (*Silene dioica*) has separate male (stamen bearing) and female (ovary and stigma bearing) plants. Both sexes can be infected by the anther smut *Ustillago violacea*. This smut produces spores in the plant’s anthers, which are then transported to other host plants by insect vectors. *Ustillago* causes their sex to change on infecting the female flowers, triggering stamen development in the genetically female flowers. In populations of *Silene* in which there is no infection by *Ustillago* the ratio of male to female flowers is 1:1. Significant amounts of infection by the fungus may be indicated if there is an increase in the proportion of apparently male flowers relative to the expected 1:1 ratio.

The two examples considered above are as simple as things get: there are only two categories (Males and Females), and we expect a 1:1 ratio. However, the χ2 goodness of fit test can be employed in more complicated situations:

* The test can be applied to any number of categories. For example, we might have a diet choice experiment where squirrels are offered four different food types in equal proportions, and the food chosen by each squirrel is recorded. The study variable would then have four categories: one for each food type.
* The expected ratio need not be 1:1. For example, the principles of Mendelian genetics predict that the offspring of two plants which are heterozygous for flower colour (white recessive, pink dominant) will be either pink or white flowered, in the ratio 3:1. Plants from a breeding experiment could be tested against this expected ratio.

## 26.2 How does the chi-square goodness of fit test work?

The χ2 goodness of fit test uses raw counts to address questions about expected proportions or probabilities of events[24](#fn24). As always, we start by setting up the appropriate null hypothesis. This will be question specific, but it must always be framed in terms of ‘no effect’ or ‘no difference’. We then work out what a sampling distribution of some kind looks like under this null hypothesis and use this to assess how likely the observed result is (i.e. calculate a *p*-value).

We don’t need to work directly with the sampling distributions of counts in each category. Instead, we calculate an appropriate χ2 test statistic. The way to think about this is that the χ2 statistic reduces the information in the separate category counts down to a single number.

Let’s see how the χ2 goodness of fit test works using the *Silene* example discussed above. Imagine that we collected data on the frequency of plants bearing male and female flowers in a population of *Silene dioica*:

|  | Male | Female |
| --- | --- | --- |
| **Observed** | 105 | 87 |

We want to test whether the ratio of male to female flowers differs significantly from that expected in an uninfected population. The ‘expected in an uninfected population’ situation is the null hypothesis for the test.

**Step 1.** Calculate the counts expected when the null hypothesis is correct. This is the critical step. In the *Silene* example, we need to determine how many male and female plants we expected to sample if the sex ratio really was 1:1.

**Step 2.** Calculate the χ2 test statistic from the observed and expected counts. We will show this later. However, the calculation isn’t all that important, in the sense that we don’t learn much by doing it. The resulting χ2 statistic summarises—across all the categories—how likely the observed data are under the null hypothesis.

**Step 3.** Compare the χ2 statistic to the theoretical predictions of the χ2 distribution to assess the statistical significance of the difference between observed and expected counts.

The interpretation of this *p*-value in this test is the same as for any other kind of statistical test: it is the probability we would see the observed frequencies, or more extreme values, under the null hypothesis.

### 26.2.1 Assumptions of the chi-square goodness of fit test

Let’s remind ourselves about the assumptions of the χ2 goodness of fit test:

1. The data are independent counts of objects or events which can be classified into mutually exclusive categories. We shouldn’t aggregate *Silene* sex data from different surveys unless we were certain each survey had sampled different populations.
2. The expected counts are not too low. The rule of thumb is that the expected values (not the observed counts!) should be greater than 5. If any of the expected values are below 5, the *p*-values generated by the test start to become less reliable.

## 26.3 Carrying out a chi-square goodness of fit test in R

#### 

There is no need to download any data to work through this example. Although it is generally a good idea to keep data and code separate, the data used in a χ2 goodness of fit test are so simple we sometimes keep them in our R code. That’s what we’ll do here.

Let’s see how to use R to carry out a χ2 goodness of fit test with the *Silene* sex data. .

The first step is to construct a numeric vector containing the *observed frequencies* for each category. We’ll call this `observed_freqs`:

```
observed_freqs <- c(105, 87)
observed_freqs
```

```
## [1] 105  87
```

Note that `observed_freqs` is a numeric vector, not a data frame. It doesn’t matter what order the two category counts are supplied in.

Next, we need to calculate the *expected frequencies* of each category. Rather than expressing these as counts, R expects these to be proportions. We need to construct a second numeric vector containing this information. We’ll call this `expected_probs`:

```
# calculate the number of categories
n_cat <- length(observed_freqs)
# calculate the number in each category (equal frequencies in this e.g.)
expected_probs <- rep(1, n_cat) / n_cat
expected_probs
```

```
## [1] 0.5 0.5
```

Finally, use the `chisq.test` function to calculate the χ2 value, degrees of freedom and *p*-value for the test. The first argument is the numeric vector of observed counts (the data), and the second is the expected proportions in each category:

```
chisq.test(observed_freqs, p = expected_probs)
```

```
## 
##  Chi-squared test for given probabilities
## 
## data:  observed_freqs
## X-squared = 1.6875, df = 1, p-value = 0.1939
```

The vectors containing the data and expected proportions have to be the same length for this to work. R will complain and throw an error otherwise.

The output is straightforward to interpret. The first part (`Chi-squared test for given probabilities`) reminds us what kind of test we did. The phrase ‘for given probabilities’ is R-speak informing us that we have carried out a goodness of fit test. The next line (`data: observed_freqs`) reminds us what data we used. The final line is the one we care about: `X-squared = 1.6875, df = 1, p-value = 0.1939`. This output shows us the χ2 test statistic, the degrees of freedom associated with the test, and the *p*-value. Since *p* > 0.05, we conclude that the sex ratio is not significantly different from a 1:1 ratio.

#### Degrees of freedom for a χ2 goodness of fit test

We always need to know the degrees of freedom associated with our tests. This is n−1 in a χ2 goodness-of-fit test, where n is the number of categories. This comes from the fact that we have to calculate one expected frequency per category (= n frequencies), but since the frequencies have to add up to the total number of observations, once we know n−1 frequencies, the last one is fixed.

### 26.3.1 Summarising the result

Having obtained the result we need to write the conclusion. As always, we always go back to the original question to write the conclusion. In this case the appropriate conclusion is:

> The observed sex ratio of *Silene dioica* does not differ significantly from the expected 1:1 ratio (χ2 = 1.69, d.f. = 1, *p* = 0.19).

### 26.3.2 A bit more about goodness of fit tests in R

There is a useful short cut that we can employ when we expect equal numbers in every category (as above). There is no need to calculate the expected proportions in this situation because R will assume we meant to use the ‘equal frequencies’ null hypothesis. So the following…

```
chisq.test(observed_freqs)
```

```
## 
##  Chi-squared test for given probabilities
## 
## data:  observed_freqs
## X-squared = 1.6875, df = 1, p-value = 0.1939
```

…is exactly equivalent to the longer method we used first. We showed the first method because we do sometimes need to carry out a goodness of fit test assuming unequal expected values.

R will warn us if it thinks the data are not suitable for a χ2 test. Though it is often safe to ignore warnings produced by R, this is a situation where it is important to pay attention to it. We can see what the relevant warning looks like by using `chisq.test` with a fake data set:

```
chisq.test(c(2,5,7,5,5))
```

```
## Warning in chisq.test(c(2, 5, 7, 5, 5)): Chi-squared approximation may be
## incorrect
```

```
## 
##  Chi-squared test for given probabilities
## 
## data:  c(2, 5, 7, 5, 5)
## X-squared = 2.6667, df = 4, p-value = 0.6151
```

The warning `Chi-squared approximation may be incorrect` is telling us that there might be a problem with the test. What is the problem? The expected counts are below 5, which means the *p*-values produced by `chisq.test` will not be reliable.

### 26.3.3 Doing it the long way…

People tend not to carry out χ2 goodness of fit tests by hand these days. Why bother when a computer can do it for us? Nonetheless, we promised to show the calculation of the χ2 test statistic. It is found by taking the difference of each observed and expected count, squaring those differences, dividing each of these squared differences by the expected frequency, and finally, summing these numbers over the categories. That’s what this formula for the χ2 test statistic means:

χ2=∑(Oi−Ei)2Ei

In this formula the Oi are the observed frequencies, the Ei are the expected frequencies, the i in Oi and Ei refer to the different categories, and the Σ means summation (‘add up’).

Does that look formula familiar? It looks a lot like the quantity we calculated when comparing plant morph counts to their expected value the previous [categorical data](categorical-data-intro-chapter.html#categorical-data-intro-chapter) chapter. That’s no accident.

Here’s how the calculation of the χ2 test statistic works with our current example:

1. Calculate out how many male and female plants we would have sampled on average if the sex ratio really were 1:1. We already did this—the numbers are: (105 + 87)/2 = 192/2 = 96 of each sex.
2. Calculate the (Oi−Ei)2Ei term associated with each category. For the males, we have (105-96)2 / 96 = 0.844, and for females we (87-96)2 / 96 = 0.844[25](#fn25). The χ2 statistic is the sum of these two values.

Once we have obtained the χ2 value and the degrees of freedom for the test (d.f. = 1), we can calculate its *p*-value. We won’t show how to do this—it’s much simpler to get the whole job done in one step using `chisq.test`.

## 26.4 Determining appropriate expected values

A goodness of fit test can only be applied if we can determine some expected values with which to compare the observed counts. Equally obviously, almost anything can be made significant by using inappropriate expected values. This means we always need to have a justifiable basis for choosing the expected values. In many cases the experiment can be designed, or the data collected, in such a way that we would expect to find equal numbers in each category if whatever it is we are interested in is not having an effect. At other times, the expectation can be generated by knowledge, or prediction, of a biological process e.g. a 1:1 sex ratio, a 3:1 ratio of phenotypes. At other times the expectation may need a bit more working out.

---

24. We said it before, but it’s worth saying again. Do not apply the goodness of fit test to proportions or means.[↩︎](goodness-of-fit-tests.html#fnref24)
25. These are the same because there are only two categories in this example.[↩︎](goodness-of-fit-tests.html#fnref25)