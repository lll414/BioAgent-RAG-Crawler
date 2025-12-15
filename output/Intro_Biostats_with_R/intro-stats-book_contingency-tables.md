Source URL: https://tuos-bio-data-skills.github.io/intro-stats-book/contingency-tables.html
Date Scraped: 2025-12-12

---

# Chapter 27 Contingency tables

## 27.1 When do we use a chi-square contingency table test?

A χ2 contingency table test is appropriate in situations where we are studying *two or more* categorical variables—each object is classified according to more than one categorical variable—and we want to evaluate whether categories are associated. Here are a couple of examples:

* Returning to the data on biology students, we suggested that we might want to know if eye colour was related to sex in any way. That is, we want to know whether brown and blue eyes occur in different proportions among males and females. If we recorded the sex and eye colour of male and female students, we could use a χ2 contingency table test to evaluate whether eye colour is associated with sex.
* The two-spot ladybird (*Adalia bipunctata*) occurs in two forms: the typical form, which is red with black spots and the dark form, which has much of the elytral surface black, with the two spots red. We want to know whether the relative frequency of the colour morphs is different in industrial and rural habitats. We could address this question by applying a χ2 contingency table test to an aggregate sample of two-spot ladybirds taken from rural and industrial habitats.

Let’s think about how these kinds of data look. Here are the biology student sex and eye colour data again, organised into a table:

|  | Blue eyes | Brown eyes |
| --- | --- | --- |
| **Male** | 22 | 10 |
| **Female** | 29 | 21 |

This is called a two-way contingency table. It is a *two-way* contingency table because it summarises the frequency distribution of two categorical variables at the same time[26](#fn26). If we had measured three variables, we would have a *three-way* contingency table (e.g. 2 x 2 x 2).

A contingency table takes its name from the fact that it captures the ‘contingencies’ among the categorical variables: it summarises how the frequencies of one categorical variable are associated with the categories of another. The term association is used here to describe the non-independence of categories among categorical variables. Other terms used to refer to the same idea include ‘linkage’, ‘non-independence’, and ‘interaction’.

Associations are evident when the proportions of objects in one set of categories (e.g. R1 and R2) depends on a second set of categories (e.g. C1 and C2). Here are two possibilities:

Table 27.1: Contingency tables without an association (left table),
and with an association (right table), among two categorical variables

| |  | C1 | C2 | | --- | --- | --- | | **R1** | 10 | 20 | | **R2** | 40 | 80 | | |  | C1 | C2 | | --- | --- | --- | | **R1** | 10 | 80 | | **R2** | 40 | 20 | |

* There is no evidence of association in the first table (left): the numbers in category R1 are 1/4 of those in category R2, whether the observations are in category C1 or C2. Notice that this reasoning isn’t about the total numbers in each category—there are 100 category C2 cases and only 50 category C1 cases.
* There is evidence of an association in the second table (right): the proportion of observations in R1 changes markedly depending on whether we are looking at observations for category C1 or category C2. The R1 cases are less frequent in the C1 column relative to the C2 column. Again, it is the proportions that matter, not the raw numbers.

We can ask a variety of different questions of data in a contingency table, but they are usually used to test for associations between the variables they summarise. That’s the topic of this chapter. There are different ways to carry out such tests of association. We’ll focus on the most widely used tool—the ‘Pearson’s Chi Square’ ($\chi^{2}$) contingency table test[27](#fn27).

## 27.2 How does the chi-square contingency table test work?

The χ2 contingency table test uses data in the form of a contingency table to address questions about the dependence between two or more different *kinds* of outcomes, or events. We start to tackle this question by setting up the appropriate null hypothesis. The null hypothesis is always the same for the standard contingency table test of association: it proposes that different events are independent of one another. This means the occurrence of one kind of event does not depend on the other type of event, i.e. they are *not* associated.

Once the null hypothesis has been worked out, the remaining calculations are no different from those used in a goodness of fit test. We calculate the frequencies expected in each cell under the null hypothesis, we calculate a χ2 test statistic to summarise the mismatch between observed and expected values and then use this to assess how likely the observed result is under the null hypothesis, resulting in the *p*-value.

We’ll continue with the two-spot ladybird (*Adalia bipunctata*) example from the beginning of this chapter. The dark (melanic) form is under the control of a single gene. Melanic and red types occur at different frequencies in different areas. Two observations are pertinent to this study:

1. In London melanics comprise about 10%, whereas in rural towns in northern England the frequency is greater (e.g. Harrogate 63%, Hexham 75%).
2. The frequency of melanics has decreased in Birmingham since smoke control legislation was introduced.

It was thought that the different forms might be differentially susceptible to some toxic component of smoke, but this doesn’t explain the geographic variation in the proportions of melanics. It turns out that the effect is a subtle one in which melanic forms do better in conditions of lower sunshine than red forms due to their greater ability to absorb solar radiation. That means melanics are favoured where the climate is naturally less sunny, though there will also be smaller-scale variations due to local environmental conditions that affect solar radiation, such as smoke.

A survey was carried out of *Adalia bipunctata* in a large urban area and the more rural surrounding areas to test whether this effect still occurs in industrial areas. The following frequencies of different colour forms were obtained.

|  | Black | Red | Totals |
| --- | --- | --- | --- |
| **Rural** | 30 | 70 | *100* |
| **Industrial** | 115 | 85 | *200* |
| **Totals** | *145* | *155* | *300* |

We want to test whether the proportions of melanics are different between urban and rural areas. In order to make it a bit easier to discuss the calculations involved, we’ll refer to each cell in the table by a letter…

|  | Black | Red | Totals |
| --- | --- | --- | --- |
| **Rural** | a | b | e |
| **Industrial** | c | d | f |
| **Totals** | g | h | k |

We can now step through the steps involved in a χ2 contingency table test:

**Step 1.** We need to work out the expected numbers in cells a-d, under the null hypothesis that the two kinds of outcomes (colour and habitat type) are independent of one another. Let’s see how it works for the Black-Rural cell (‘a’):

* Calculate the probability that a randomly chosen individual in the sample is from a rural location (p(Rural)). This is the total ‘Rural’ count (e) divided by the grand total (k):

p(Rural)=ek

* Calculate the probability that a randomly chosen individual in the sample set is black (p(Black)). This is the total ‘Black’ count (g) divided by the grand total (k):

p(Black)=gk

* Calculate the probability that a randomly chosen individual is both ‘Rural’ and ‘Black’, *assuming* these events are ‘independent’ . This is given by the product of the probabilities from steps 1 and 2:

p(Rural, Black)=ek×gk

* Convert this to the expected number of individuals that are ‘Rural’ and ‘Black’, under the independence assumption. This is the probability from step 3 multiplied by the grand total (k):

E(Rural, Black)=ek×gk×k=g×ek

In general, the expected value for any particular cell in a contingency table is given by multiplying the associated row and column totals and then dividing by the grand total.

**Step 2** Calculate the χ2 test statistic from the four observed cell counts and their expected values. The resulting χ2 statistic summarises—across all the cells—how likely the observed frequencies are under the null hypothesis of no association (= independence).

**Step 3** Compare the χ2 statistic to the theoretical predictions of the χ2 distribution in order to assess the statistical significance of the difference between observed and expected counts. This *p*-value is the probability we would see the observed pattern of cell counts, or a more strongly associated pattern, under the null hypothesis of no association. A low *p*-value represents evidence against the null hypothesis.

### 27.2.1 Assumptions of the chi-square contingency table test

The assumptions of the χ2 goodness of fit test are essentially the same as those of the goodness of fit test:

1. The observations are independent counts. For example, it would not be appropriate to apply a χ2 goodness of fit test where observations are taken before and after an experimental intervention is applied to the same objects.
2. The expected counts are not too low. The rule of thumb is that the expected values (again, not the observed counts) should be greater than 5. If any of the expected values are below 5, the *p*-values generated by the test start to become less reliable.

The second assumption is not as critical as some people imagine. If we find ourselves in a situation where the expected value in one or more cells is less than 5 there is a simple fix—use a different test. We’ll highlight one way to achieve this below.

## 27.3 Carrying out a chi-square contingency table test in R

We’ll use the ladybird colour and habitat type example to demonstrate how to conduct a chi-square contingency table test in R.

#### 

Three versions of the ladybird data will be used: LADYBIRDS1.CSV, LADYBIRDS2.CSV, and LADYBIRDS3.CSV. These all contain the same information. They only differ in how it is organised. We assume that each of these has been read into tibbles called `lady_bird_df1`, `lady_bird_df2` and `lady_bird_df3`. For example:

```
lady_bird_df1 <- read_csv(file = "LADYBIRDS1.CSV")
lady_bird_df2 <- read_csv(file = "LADYBIRDS2.CSV")
lady_bird_df3 <- read_csv(file = "LADYBIRDS3.CSV")
```

Carrying out a χ2 contingency table test in R is very simple: we use the `chisq.test` function again. The only snag is that we have to ensure the data are correctly formatted. When we read data into R using a function like `read_csv` (or `read.csv`), we end up with a tibble (or ordinary data frame). Unfortunately, the `chisq.test` function is one of the few statistical functions not designed to work with data frames. We have to convert our data to a special contingency table object to use it with `chisq.test`.[28](#fn28)

How do we construct a contingency table object? We use a function called `xtabs`. The `xtabs` function does categorical ‘cross tabulation’—it sums up the number of occurrences of different combinations of categories among variables. We need to see how to use `xtabs` before running through the chi-square contingency table test.

### 27.3.1 Step 1. Get the data into the correct format

The usage of `xtabs` depends upon how the raw data are organised. There are possible cases.

#### 27.3.1.1 **Case 1: Raw data**

Data suitable for analysis with a χ2 contingency table test are often represented in a data set with one column per categorical variable and one row per observation. This is the version of the ladybird data in the `lady_bird_df1`:

```
glimpse(lady_bird_df1)
```

```
## Rows: 300
## Columns: 2
## $ Habitat <chr> "Rural", "Rural", "Rural", "Rural", "Rural", "Rural", "Rural",…
## $ Colour  <chr> "Black", "Black", "Black", "Black", "Black", "Black", "Black",…
```

This data frame contains 300 rows (one for each ladybirds) and two columns (`Habitat` and `Colour`). Each column is a categorical variable containing information about the ladybirds. Each one has two unique values:

```
unique(lady_bird_df1$Habitat)
```

```
## [1] "Rural"      "Industrial"
```

```
unique(lady_bird_df1$Colour)
```

```
## [1] "Black" "Red"
```

We require a two-way table that contains the total counts in each combination of categories. This is what `xtabs` does. It takes two arguments: the first is a formula (involving the `~` symbol) that specifies the required contingency table; the second is the name of the data frame containing the raw data. When working with data in the above format—one observation per row—we use a formula that only contains the names of the categorical variables on the right hand side of the `~` (`~ Habitat + Colour`):

```
lady_bird_table <- xtabs(~ Habitat + Colour, data = lady_bird_df1)
```

When used like this, `xtabs` sums up the number of observations with each combination of `Habitat` and `Colour`.

We called the output `lady_bird_table` to emphasise that the data from `xtabs` are now stored in a contingency table. When we print this to the console we see that `lady_bird_table` does indeed refer to something that looks like a 2 x 2 contingency table of counts:

```
lady_bird_table
```

```
##             Colour
## Habitat      Black Red
##   Industrial   115  85
##   Rural         30  70
```

#### 27.3.1.2 **Case 2: Partial counts**

Sometimes data suitable for analysis with a χ2 contingency table test are partially summarised into counts. For example, imagine that we visited five rural and five urban sites and recorded the numbers of red and black colour forms found at each site. This is the version of the ladybird data in the `lady_bird_df2`:

```
glimpse(lady_bird_df2)
```

```
## Rows: 20
## Columns: 4
## $ Habitat <chr> "Rural", "Rural", "Rural", "Rural", "Rural", "Rural", "Rural",…
## $ Site    <chr> "R1", "R2", "R3", "R4", "R5", "R1", "R2", "R3", "R4", "R5", "U…
## $ Colour  <chr> "Black", "Black", "Black", "Black", "Black", "Red", "Red", "Re…
## $ Number  <dbl> 10, 3, 4, 7, 6, 15, 18, 9, 12, 16, 32, 25, 25, 17, 16, 17, 23,…
```

It’s a bit easier to see how the data are organised if we print the `lady_bird_df2` object:

```
lady_bird_df2
```

```
## # A tibble: 20 × 4
##    Habitat    Site  Colour Number
##    <chr>      <chr> <chr>   <dbl>
##  1 Rural      R1    Black      10
##  2 Rural      R2    Black       3
##  3 Rural      R3    Black       4
##  4 Rural      R4    Black       7
##  5 Rural      R5    Black       6
##  6 Rural      R1    Red        15
##  7 Rural      R2    Red        18
##  8 Rural      R3    Red         9
##  9 Rural      R4    Red        12
## 10 Rural      R5    Red        16
## 11 Industrial U1    Black      32
## 12 Industrial U2    Black      25
## 13 Industrial U3    Black      25
## 14 Industrial U4    Black      17
## 15 Industrial U5    Black      16
## 16 Industrial U1    Red        17
## 17 Industrial U2    Red        23
## 18 Industrial U3    Red        21
## 19 Industrial U4    Red         9
## 20 Industrial U5    Red        15
```

The counts at each site are in the `Number` variable, and the site identities are in the `Site` variable. We need to sum over the sites to get the total number within each combination of `Habitat` and `Colour`. We use `xtabs` again, but this time we also tell it which variable to sum over (`Number ~ Habitat + Colour`):

```
lady_bird_table <- xtabs(Number ~ Habitat + Colour, data = lady_bird_df2)
```

When working with partially summed counts data we have to use a formula where the name of the variable containing the counts is on left hand side of the `~`, and the names of the categorical variables to sum over are on the right hand side of the `~`. The `lady_bird_table` object produced by `xtabs` is no different than before:

```
lady_bird_table
```

```
##             Colour
## Habitat      Black Red
##   Industrial   115  85
##   Rural         30  70
```

Good! These are the same data.

#### 27.3.1.3 **Case 3: Complete counts**

Data suitable for analysis with a χ2 contingency table test are sometimes already summed into the required total counts. This is the version of the ladybird data in the `lady_bird_df3`:

```
glimpse(lady_bird_df3)
```

```
## Rows: 4
## Columns: 3
## $ Habitat <chr> "Industrial", "Industrial", "Rural", "Rural"
## $ Colour  <chr> "Black", "Red", "Black", "Red"
## $ Number  <dbl> 115, 85, 30, 70
```

Again, we can see how the data are organised by printing the `lady_bird_df3` object:

```
lady_bird_df3
```

```
## # A tibble: 4 × 3
##   Habitat    Colour Number
##   <chr>      <chr>   <dbl>
## 1 Industrial Black     115
## 2 Industrial Red        85
## 3 Rural      Black      30
## 4 Rural      Red        70
```

The total counts are already in the `Number` variable, meaning there is no need to sum over anything to get the total for each combination of `Habitat` and `Colour`.

However, we still need to convert the data to a contingency table object so we can feed it to `chisq.test`. One way to do this is by using `xtabs` again, exactly as before:

```
lady_bird_table <- xtabs(Number ~ Habitat + Colour, data = lady_bird_df3)
lady_bird_table
```

```
##             Colour
## Habitat      Black Red
##   Industrial   115  85
##   Rural         30  70
```

In this case `xtabs` doesn’t change the data because it’s just ‘summing’ over one value in each combination of categories. We’re only using `xtabs` to convert it to a contingency table object. The resulting `lady_bird_table` object is the same as before.

### 27.3.2 Step 2. Do the test

Once we have the data in the form of a contingency table the associated χ2 test of independence between the two categorical variables is easy to carry out:

```
chisq.test(lady_bird_table)
```

```
## 
##  Pearson's Chi-squared test with Yates' continuity correction
## 
## data:  lady_bird_table
## X-squared = 19.103, df = 1, p-value = 1.239e-05
```

That’s it. Just pass one argument to `chisq.test`: the contingency table.

This output should make sense in the light of what we saw in the previous chapter. R first prints a reminder of the test employed (`Pearson's Chi-squared test with Yates' continuity correction`) and the data used (`data: lady_bird_table`). We’ll come back to the “Yates’ continuity correction” bit in a moment. R then summarises the χ2 value, the degrees of freedom, and the *p*-value: `X-squared = 19.103, df = 1, p-value = 1.239e-05` The *p*-value is highly significant (*p*<0.001) indicating that the colour type frequency varies among the two kinds of habitats.[29](#fn29)

#### Degrees of freedom for a χ2 contingency table test

We need to know the degrees of freedom associated with the test: in a two-way χ2 contingency table test these are (nA−1)×(nB−1), where nA−1 is the number of categories in the first variable, and nB−1 is the number of categories in the second variable. So if we’re working with a 2 x 2 table the d.f. are 1, if we’re working with a 2 x 3 table the d.f. are 2, if we’re working with a 3 x 3 table the d.f. are 4, and so on.

What was that “Yates’ continuity correction” all about? The reasoning behind using this correction is a bit beyond this book, but in a nutshell, it generates more reliable *p*-values under certain circumstances. By default, the `chisq.test` function applies this correction to all 2 x 2 contingency tables. We can force R to use the standard calculation by setting `correct = FALSE` if we want to:

```
chisq.test(lady_bird_table, correct = FALSE)
```

```
## 
##  Pearson's Chi-squared test
## 
## data:  lady_bird_table
## X-squared = 20.189, df = 1, p-value = 7.015e-06
```

Both methods give similar results in this example, though they aren’t exactly the same—the χ2 value calculated when `correct = FALSE` is very slightly higher than the value found when using the correction.

Don’t use the `correct = FALSE` option! The default correction is a safer option for 2 x 2 tables.

#### What should we do when the expecetd values are below 5?

The expected cell counts all need to be greater than 5 when using the standard χ2 contingency table test. If thsi assumption is not met, the χ2 test statistic will deviate from its theoretical distribution and the test becomes unreliable.

Fortunately, if we find ourselves in a situation with low cell counts, `chisq.test` can be used to carry out an alternative test. We can access this alternative by setting the `simulate.p.value` argument to `TRUE`:

```
chisq.test(lady_bird_table, simulate.p.value = TRUE)
```

```
## 
##  Pearson's Chi-squared test with simulated p-value (based on 2000
##  replicates)
## 
## data:  lady_bird_table
## X-squared = 20.189, df = NA, p-value = 0.0004998
```

What just happened? R used a simulation approach to construct the null distribution for the test. The details don’t matter too much, but logic of the approach is a lot like bootstrapping or permutation tests. R generates many hypothetical sets of cell counts based on the null hypothesis of no association, and then compares the distribution of their χ2 statistics to the value from the true data to calculate a *p*-value.

### 27.3.3 Summarising the result

We have obtained the result, so now we need to write the conclusion. As always, we go back to the original question to write the conclusion. In this case the appropriate conclusion is:

> There is a significant association between the colour of *Adalia bipunctata* individuals and habitat, such that black individuals are more likely to be found in industrial areas (χ2 = 19.1, d.f. = 1, *p* < 0.001).

Notice that we summarised the nature of the association alongside the statistical result. This is easy to do in the text when describing the results of a 2 x 2 contingency table test. It’s much harder to summarise the association in written form when working with larger tables. Instead, we often present a table or a bar chart showing the observed counts.

## 27.4 Working with larger tables

Interpreting the results of a contingency table test is fairly straightforward in the case of a 2 x 2 table. We can certainly use contingency table tests with larger two-way tables (e.g. 3 x 2, 3 x 3, 3 x 4, etc), and higher-dimensional tables (e.g. 2 x 2 x 2), but the results become progressively harder to understand as the tables increase in size. If you get a significant result, it is often best to compare the observed and expected counts for each cell and look for the highest differences to try and establish what is driving the significant association. Visualising the data with a bar chart can also help with interpretation.

There are ways of subdividing large tables to make subsequent χ2 tests on individual parts of a table to establish specific effects, but these are not detailed here. Note that, unless we had planned to make more detailed comparisons before collecting the data, it is hard to justify such *post hoc* analysis.

---

26. This is called their ‘joint distribution’, in case you were wondering.[↩︎](contingency-tables.html#fnref26)
27. This is the same Pearson who invented the correlation coefficient for measuring linear associations by the way.[↩︎](contingency-tables.html#fnref27)
28. The table objects produced by `xtabs` are **not** the same as the **dplyr** table-like objects: ‘tibbles’. This is one reason why **dplyr** adopted the name ‘tibble’. The overlap in names is unfortunate, but we’ll have to live with it—there are only so many ways to name things that look like tables.[↩︎](contingency-tables.html#fnref28)
29. We could have summarised the result as: habitat type varies among the two colour types. This way of explaining the result seems odd though. Ladybirds are found within habitats, not the other way around. Just keep in mind that this is a semantic issue. The contingency table test doesn’t make a distinction between directions of effects.[↩︎](contingency-tables.html#fnref29)