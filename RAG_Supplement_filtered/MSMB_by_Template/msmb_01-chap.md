Source URL: https://www.huber.embl.de/msmb/01-chap.html
Date Scraped: 2025-12-12

---

# 1  Generative Models for Discrete Data

[![](imgs/Pile_ou_face.png)](https://www.huber.embl.de/msmb/imgs/Pile_ou_face.png)

In molecular biology, many situations involve counting events: how many codons use a certain spelling, how many reads of DNA match a reference, how many CG digrams are observed in a DNA sequence. These counts give us *discrete* variables, as opposed to quantities such as mass and intensity that are measured on *continuous* scales.

If we know the rules that the mechanisms under study follow, even if the outcomes are random, we can generate the probabilities of any events we are interested in by computations and standard probability laws. This is a *top-down* approach based on deduction and our knowledge of how to manipulate probabilities. In [Chapter 2](https://www.huber.embl.de/msmb/02-chap.html), you will see how to combine this with data-driven (*bottom-up*) statistical modeling.

## 1.1 Goals for this chapter

In this chapter we will:

* Learn how to obtain the probabilities of all possible outcomes from a given model and see how we can compare the theoretical frequencies with those observed in real data.
* Explore a complete example of how to use the Poisson distribution to analyse data on epitope detection.
* See how we can experiment with the most useful generative models for discrete data: Poisson, binomial, multinomial.
* Use the R functions for computing probabilities and counting rare events.
* Generate random numbers from specified distributions.

## 1.2 A real example

Let’s dive into an example where we have a probability model for the data generating process. Our model says that mutations along the genome of HIV (Human Immunodeficiency Virus) occur with a rate of \(5 \times 10^{-4}\) per nucleotide per replication cycle. The rate is the same at each nucleotide position, and mutations at one position happen independently of what happens at other positions[1](https://www.huber.embl.de/msmb/01-chap.html#fn1). The genome size of HIV is about \(10^4=10,000\) nucleotides, thus, after one cycle, the total number of mutations will follow a **Poisson** distribution[2](https://www.huber.embl.de/msmb/01-chap.html#fn2) with rate \(5 \times 10^{-4} \times 10^4 = 5\). What does that tell us?

1 In practice, and strictly speaking, complete and utter independence will rarely hold in reality, if you look close enough. Thus, what modellers usually mean with such assertions is that any possible correlations or dependencies are so weak and rare that ignoring them is a good enough approximation.

2 We will give more details later about this type of probability distribution

This probability model predicts that the number of mutations over one replication cycle will be close to 5, and that the variability of this estimate is \(\sqrt{5}\) (the standard error). We now have baseline reference values for both the number of mutations we expect to see in a typical HIV strain and its variability.

In fact, we can deduce even more detailed information. If we want to know how often 3 mutations could occur under the Poisson(5) model, we can use an R function to generate the probability of seeing \(x=3\) events, taking the value of the **rate parameter** of the Poisson distribution, called lambda (\(\lambda\)), to be \(5\).

Greek letters such as \(\lambda\) and \(\mu\) often denote important parameters that characterize the probability distributions we use.

```
dpois(x = 3, lambda = 5)
```

```
[1] 0.1403739
```

This says the chance of seeing exactly three events is around 0.14, or about 1 in 7.

If we want to generate the probabilities of all values from 0 to 12, we do not need to write a loop. We can simply set the first argument to be the **vector** of these 13 values, using R’s sequence operator, the colon “`:`”. We can see the probabilities by plotting them ([Figure 1.1](https://www.huber.embl.de/msmb/01-chap.html#fig-Poisson5)). As with this figure, most figures in the margins of this book are created by the code shown in the text.

[![Note how the output from R is formatted: the first line begins with the first item in the vector, hence the [1], and the second line begins with the 9th item, hence the [9]. This helps you keep track of elements in long vectors. The term vector is R parlance for an ordered list of elements of the same type (in this case, numbers).](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "Note how the output from R is formatted: the first line begins with the first item in the vector, hence the [1], and the second line begins with the 9th item, hence the [9]. This helps you keep track of elements in long vectors. The term vector is R parlance for an ordered list of elements of the same type (in this case, numbers).")

Note how the output from R is formatted: the first line begins with the first item in the vector, hence the [1], and the second line begins with the 9th item, hence the [9]. This helps you keep track of elements in long vectors. The term *vector* is R parlance for an ordered list of elements of the same type (in this case, numbers).

```
0:12
```

```
 [1]  0  1  2  3  4  5  6  7  8  9 10 11 12
```

```
dpois(x = 0:12, lambda = 5)
```

```
 [1] 0.0067 0.0337 0.0842 0.1404 0.1755 0.1755 0.1462 0.1044 0.0653 0.0363
[11] 0.0181 0.0082 0.0034
```

```
barplot(dpois(0:12, 5), names.arg = 0:12, col = "red")
```

[![](01-chap_files/figure-html/fig-Poisson5-1.png)](https://www.huber.embl.de/msmb/01-chap_files/figure-html/fig-Poisson5-1.png "Figure 1.1: Probabilities of seeing 0,1,2,…,12 mutations, as modeled by the Poisson(5) distribution. The plot shows that we will often see 4 or 5 mutations but rarely as many as 12. The distribution continues to higher numbers (13,...), but the probabilities will be successively smaller, and here we don’t visualize them.")

Figure 1.1: Probabilities of seeing 0,1,2,…,12 mutations, as modeled by the Poisson(5) distribution. The plot shows that we will often see 4 or 5 mutations but rarely as many as 12. The distribution continues to higher numbers (\(13,...\)), but the probabilities will be successively smaller, and here we don’t visualize them.

Mathematical theory tells us that the Poisson probability of seeing the value \(x\) is given by the formula \(e^{-\lambda} \lambda^x / x!\). In this book, we’ll discuss theory from time to time, but give preference to displaying concrete numeric examples and visualizations like [Figure 1.1](https://www.huber.embl.de/msmb/01-chap.html#fig-Poisson5).

The Poisson distribution is a good model for rare events such as mutations. Other useful probability models for **discrete events** are the Bernoulli, binomial and multinomial distributions. We will explore these models in this chapter.

## 1.3 Using discrete probability models

[![Think of a categorical variable as having different alternative values. These are the levels, similar to the different alternatives at a gene locus: alleles.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "Think of a categorical variable as having different alternative values. These are the levels, similar to the different alternatives at a gene locus: alleles.")

Think of a categorical variable as having different alternative values. These are the levels, similar to the different alternatives at a gene locus: *alleles*.

A point mutation can either occur or not; it is a binary event. The two possible outcomes (yes, no) are called the **levels** of the categorical variable.

Not all events are binary. For example, the genotypes in a diploid organism can take three levels (AA, Aa, aa).

Sometimes the number of levels in a categorical variable is very large; examples include the number of different types of bacteria in a biological sample (hundreds or thousands) and the number of codons formed of 3 nucleotides (64 levels).

When we measure a categorical variable on a sample, we often want to tally the frequencies of the different levels in a vector of counts. R has a special encoding for categorical variables and calls them **factors**[3](https://www.huber.embl.de/msmb/01-chap.html#fn3). Here we capture the different blood genotypes for 19 subjects in a vector which we tabulate.

3 R makes sure that the factor variable will accept no other, “illegal” values, and this is useful for keeping your calculations safe.

[![c() is one of the most basic functions. It collates elements of the same type into a vector. In the code shown here, the elements of genotype are character strings.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "c() is one of the most basic functions. It collates elements of the same type into a vector. In the code shown here, the elements of genotype are character strings.")

`c()` is one of the most basic functions. It collates elements of the same type into a vector. In the code shown here, the elements of `genotype` are character strings.

```
genotype = c("AA","AO","BB","AO","OO","AO","AA","BO","BO",
             "AO","BB","AO","BO","AB","OO","AB","BB","AO","AO")
table(genotype)
```

```
genotype
AA AB AO BB BO OO 
 2  2  7  3  3  2
```

On creating a *factor*, R automatically detects the levels. You can access the levels with the `levels` function.

[![It is not obvious from the output of the table function that the input was a factor; however if there had been another level with no instances, the table would also have contained that level, with a zero count.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "It is not obvious from the output of the table function that the input was a factor; however if there had been another level with no instances, the table would also have contained that level, with a zero count.")

It is not obvious from the output of the `table` function that the input was a factor; however if there had been another level with no instances, the table would also have contained that level, with a zero count.

```
genotypeF = factor(genotype)
levels(genotypeF)
```

```
[1] "AA" "AB" "AO" "BB" "BO" "OO"
```

```
table(genotypeF)
```

```
genotypeF
AA AB AO BB BO OO 
 2  2  7  3  3  2
```

Question 1.1

What if you want to create a *factor* that has some levels not yet in your data?

Solution

Look at the manual page of the `factor` function.

If the order in which the data are observed doesn’t matter, we call the random variable **exchangeable**. In that case, all the information available in the factor is summarized by the counts of the factor levels. We then say that the vector of frequencies is **sufficient** to capture all the relevant information in the data, thus providing an effective way of compressing the data.

### 1.3.1 Bernoulli trials

[![](imgs/BallsinBoxes2.png)](https://www.huber.embl.de/msmb/imgs/BallsinBoxes2.png "Figure 1.2: Two possible events with unequal probabilities. We model this by a Bernoulli distribution with probability parameter p=2/3.")

Figure 1.2: Two possible events with unequal probabilities. We model this by a Bernoulli distribution with probability parameter \(p=2/3\).

Tossing a coin has two possible outcomes. This simple experiment, called a Bernoulli trial, is modeled using a so-called Bernoulli random variable. Understanding this building block will take you surprisingly far. We can use it to build more complex models.

Let’s try a few experiments to see what some of these random variables look like. We use special R functions tailored to generate outcomes for each type of distribution. They all start with the letter `r`, followed by a specification of the model, here `rbinom`, where `binom` is the abbreviation used for binomial.

Suppose we want to simulate a sequence of 15 fair coin tosses. To get the outcome of 15 Bernoulli trials with a probability of success equal to 0.5 (a fair coin), we write

```
rbinom(15, prob = 0.5, size = 1)
```

```
 [1] 0 0 0 1 1 0 0 0 0 1 0 1 1 0 0
```

We use the `rbinom` function with a specific set of **parameters**[4](https://www.huber.embl.de/msmb/01-chap.html#fn4): the first parameter is the number of trials we want to observe; here we chose 15. We designate by `prob` the probability of success. By `size=1` we declare that each individual trial consists of just one single coin toss.

4 For R functions, parameters are also called **argument**s.

Question 1.2

Repeat this function call a number of times. Why isn’t the answer always the same?

Success and failure can have unequal probabilities in a Bernoulli trial, as long as the probabilities sum to one[5](https://www.huber.embl.de/msmb/01-chap.html#fn5). To simulate twelve trials of throwing a ball into the two boxes as shown in [Figure 1.2](https://www.huber.embl.de/msmb/01-chap.html#fig-twoballs), with probability of falling in the right-hand box \(\frac{2}{3}\) and in the left-hand box \(\frac{1}{3}\), we write

5 We call such events **complementary**.

```
rbinom(12, prob = 2/3, size = 1)
```

```
 [1] 1 1 1 0 1 0 1 1 0 1 1 1
```

The 1 indicates success, meaning that the ball fell in the right-hand box, 0 means the ball fell in the left-hand box.

### 1.3.2 Binomial success counts

If we only care how many balls go in the right-hand box, then the order of the throws doesn’t matter[6](https://www.huber.embl.de/msmb/01-chap.html#fn6), and we can get this number by just taking the sum of the cells in the output vector. Therefore, instead of the binary vector we saw above, we only need to report a single number. In R, we can do this using one call to the `rbinom` function with the parameter `size` set to 12.

Two outcomes and a size of 1 or more makes it a binomial trial. If the size is 1, then this is the special case of the Bernoulli trial.

6 The exchangeability property.

```
rbinom(1, prob = 2/3, size = 12)
```

```
[1] 9
```

This output tells us how many of the twelve balls fell into the right-hand box (the outcome that has probability 2/3). We use a random two-box model when we have only two possible outcomes such as heads or tails, success or failure, CpG or non-CpG, M or F, Y = pyrimidine or R = purine, diseased or healthy, true or false. We only need the probability of “success” \(p\), because “failure” (the *complementary* event) will occur with probability \(1-p\). When looking at the result of several such trials, if they are exchangeable[7](https://www.huber.embl.de/msmb/01-chap.html#fn7), we record only the number of successes. Therefore, SSSSSFSSSSFFFSF is summarized as (#Successes=10, #Failures=5), or as \(x=10\), \(n=15\).

7 One situation in which trials are exchangeable is if they are **independent** of each other.

The number of successes in 15 Bernoulli trials with a probability of success of 0.3 is called a **binomial** random variable or a random variable that follows the \(B(15,0.3)\) distribution. To generate samples, we use a call to the `rbinom` function with the number of trials set to 15:

[![What does set.seed do here?](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "What does set.seed do here?")

What does `set.seed` do here?

```
set.seed(235569515)
rbinom(1, prob = 0.3, size = 15)
```

```
[1] 5
```

Question 1.3

Repeat this function call ten times. What seems to be the most common outcome?

Solution

The most frequent value is 4. In fact, the theoretical proportion of times that we expect 4 to appear is the value of the probability that \(X=4\) if \(X\) follows \(B(15, 0.3)\).

The complete **probability mass distribution** is available by typing:

We use the function `round` to keep the number of printed decimal digits down to 2.

```
probabilities = dbinom(0:15, prob = 0.3, size = 15)
round(probabilities, 2)
```

```
 [1] 0.00 0.03 0.09 0.17 0.22 0.21 0.15 0.08 0.03 0.01 0.00 0.00 0.00 0.00 0.00
[16] 0.00
```

We can produce a bar plot of this distribution, shown in [Figure 1.3](https://www.huber.embl.de/msmb/01-chap.html#fig-binombarplot).

```
barplot(probabilities, names.arg = 0:15, col = "red")
```

[![](01-chap_files/figure-html/fig-binombarplot-1.png)](https://www.huber.embl.de/msmb/01-chap_files/figure-html/fig-binombarplot-1.png "Figure 1.3: Theoretical distribution of B(15,0.3) . The highest bar is at x=4. We have chosen to represent theoretical values in red throughout.")

Figure 1.3: Theoretical distribution of \(B(15,0.3)\) . The highest bar is at \(x=4\). We have chosen to represent theoretical values in red throughout.

The number of trials is the number we input in R as the `size` parameter and is often written \(n\), while the probability of success is \(p\). Mathematical theory tells us that for \(X\) distributed as a binomial distribution with parameters \((n,p)\) written \(X \sim B(n,p)\), the probability of seeing \(X=k\) successes is

[![Instead of \frac{n!}{(n-k)!k!} we can use the special notation {n \choose k} as a shortcut.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "Instead of \frac{n!}{(n-k)!k!} we can use the special notation {n \choose k} as a shortcut.")

Instead of \(\frac{n!}{(n-k)!k!}\) we can use the special notation \({n \choose k}\) as a shortcut.

\[
\begin{aligned}
P(X=k)
&=&\frac{n\times (n-1)... (n-k+1)}{k\times(k-1)... 1}\; p^k\, (1-p)^{n-k}\\
&=&\frac{n!}{(n-k)!k!}\;p^k\, (1-p)^{n-k}\\
&=&{ n \choose k}\; p^k\, (1-p)^{n-k}.
\end{aligned}
\]

Question 1.4

What is the output of the formula for \(k=3\), \(p=2/3\), \(n=4\)?

### 1.3.3 Poisson distributions

[![](imgs/Simeon_Poisson.jpg)](https://www.huber.embl.de/msmb/imgs/Simeon_Poisson.jpg "Figure 1.4: Simeon Poisson, after whom the Poisson distribution is named (this is why it always has a capital letter, except in our R code). ")

Figure 1.4: Simeon Poisson, after whom the Poisson distribution is named (this is why it always has a capital letter, except in our R code).

When the probability of success \(p\) is small and the number of trials \(n\) large, the binomial distribution \(B(n, p)\) can be faithfully approximated by a simpler distribution, the **Poisson distribution** with rate parameter \(\lambda=np\). We already used this fact, and this distribution, in the HIV example ([Figure 1.1](https://www.huber.embl.de/msmb/01-chap.html#fig-Poisson5)).

Question 1.5

What is the probability mass distribution of observing `0:12` mutations in a genome of \(n = 10^4\) nucleotides, when the probability is \(p = 5 \times 10^{-4}\) per nucleotide? Is it similar when modeled by the binomial \(B(n,p)\) distribution and by the Poisson\((\lambda=np)\) distribution?

Solution

Note that, unlike the binomial distribution, the Poisson no longer depends on two separate parameters \(n\) and \(p\), but only on their product \(np\). As in the case of the binomial distribution, we also have a mathematical formula for computing Poisson probabilities:

\[
P(X=k)= \frac{\lambda^k\;e^{-\lambda}}{k!}.
\]

For instance, let’s take \(\lambda=5\) and compute \(P(X=3)\):

```
5^3 * exp(-5) / factorial(3)
```

```
[1] 0.1403739
```

which we can compare with what we computed above using `dpois`.

Task

Simulate a mutation process along 10,000 positions with a mutation rate of \(5\times10^{-4}\) and count the number of mutations. Repeat this many times and plot the distribution with the `barplot` function (see [Figure 1.5](https://www.huber.embl.de/msmb/01-chap.html#fig-gen-simpoisson)).

```
rbinom(1, prob = 5e-4, size = 10000)
```

```
[1] 6
```

```
simulations = rbinom(n = 300000, prob = 5e-4, size = 10000)
barplot(table(simulations), col = "lavender")
```

[![](01-chap_files/figure-html/fig-gen-simpoisson-1.png)](https://www.huber.embl.de/msmb/01-chap_files/figure-html/fig-gen-simpoisson-1.png "Figure 1.5: Simulated distribution of B(10000, 10^{-4}) for 300000 simulations.")

Figure 1.5: Simulated distribution of B(10000, \(10^{-4}\)) for 300000 simulations.

Now we are ready to use probability calculations in a case study.

### 1.3.4 A generative model for epitope detection

When testing certain pharmaceutical compounds, it is important to detect proteins that provoke an allergic reaction. The molecular sites that are responsible for such reactions are called **epitopes**. The technical definition of an epitope is:

> A specific portion of a macromolecular antigen to which an antibody binds. In the case of a protein antigen recognized by a T-cell, the epitope or determinant is the peptide portion or site that binds to a Major Histocompatibility Complex (MHC) molecule for recognition by the T cell receptor (TCR).

And in case you’re not so familiar with immunology: an **antibody** (as schematized in [Figure 1.6](https://www.huber.embl.de/msmb/01-chap.html#fig-antibody)) is a type of protein made by certain white blood cells in response to a foreign substance in the body, which is called the **antigen**.

[![](imgs/Antibody_IgG2.png)](https://www.huber.embl.de/msmb/imgs/Antibody_IgG2.png "Figure 1.6: A diagram of an antibody showing several immunoglobulin domains in color.")

Figure 1.6: A diagram of an antibody showing several immunoglobulin domains in color.

An antibody binds (with more or less specificity) to its antigen. The purpose of the binding is to help destroy the antigen. Antibodies can work in several ways, depending on the nature of the antigen. Some antibodies destroy antigens directly. Others help recruit white blood cells to destroy the antigen. An epitope, also known as antigenic determinant, is the part of an antigen that is recognized by the immune system, specifically by antibodies, B cells or T cells.

#### ELISA error model with known parameters

ELISA[8](https://www.huber.embl.de/msmb/01-chap.html#fn8) assays are used to detect specific epitopes at different positions along a protein. Suppose the following facts hold for an ELISA array we are using:

8 **E**nzyme-**L**inked **I**mmuno**S**orbent **A**ssay ([Wikipedia link ELISA](http://en.wikipedia.org/wiki/ELISA)).

* The baseline noise level per position, or more precisely the false positive rate, is 1%. This is the probability of declaring a hit – we think we have an epitope – when there is none. We write this \(P(\text{declare epitope}|\text{no epitope})\)[9](https://www.huber.embl.de/msmb/01-chap.html#fn9).
* The protein is tested at 100 different positions, supposed to be independent.
* We are going to examine a collection of 50 patient samples.

9 The vertical bar in expressions such as \(X|Y\) means “\(X\) happens *conditional on* \(Y\) being the case”.

#### One patient’s data

The data for one patient’s assay look like this:

```
  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 [38] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 [75] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

where the 1 signifies a hit (and thus the potential for an allergic reaction), and the zeros signify no reaction at that position.

Task

Verify by simulation that the sum of 50 independent Bernoulli variables with \(p=0.01\) is –to good enough approximation– the same as a Poisson(\(0.5\)) random variable.

#### Results from the 50 assays

We’re going to study the data for all 50 patients tallied at each of the 100 positions. If there are no allergic reactions, the false positive rate means that for one patient, each individual position has a probability of 1 in 100 of being a 1. So, after tallying 50 patients, we expect at any given position the sum of the 50 observed \((0,1)\) variables to have a Poisson distribution with parameter 0.5. A typical result may look like [Figure 1.7](https://www.huber.embl.de/msmb/01-chap.html#fig-typicalP). Now suppose we see actual data as shown in [Figure 1.8](https://www.huber.embl.de/msmb/01-chap.html#fig-epitopedata), loaded as R object `e100` from the data file `e100.RData`.

[![](01-chap_files/figure-html/fig-typicalP-1.png)](https://www.huber.embl.de/msmb/01-chap_files/figure-html/fig-typicalP-1.png "Figure 1.7: Plot of typical data from our generative model for the background, i.,e., for the false positive hits: 100 positions along the protein, at each position the count is drawn from a Poisson(0.5) random variable.")

Figure 1.7: Plot of typical data from our generative model for the background, i.,e., for the false positive hits: 100 positions along the protein, at each position the count is drawn from a Poisson(0.5) random variable.

```
load("../data/e100.RData")
barplot(e100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5, 100.5),
  names.arg = seq_along(e100), col = "darkolivegreen")
```

[![](01-chap_files/figure-html/fig-epitopedata-1.png)](https://www.huber.embl.de/msmb/01-chap_files/figure-html/fig-epitopedata-1.png "Figure 1.8: Output of the ELISA array results for 50 patients in the 100 positions.")

Figure 1.8: Output of the ELISA array results for 50 patients in the 100 positions.

The spike in [Figure 1.8](https://www.huber.embl.de/msmb/01-chap.html#fig-epitopedata) is striking. *What are the chances of seeing a value as large as 7, if no epitope is present?*  
If we look for the probability of seeing a number as big as 7 (or larger) when considering one Poisson(\(0.5\)) random variable, the answer can be calculated in closed form as

\[
P(X\geq 7)= \sum\_{k=7}^\infty P(X=k).
\]

This is, of course, the same as \(1-P(X\leq 6)\). The probability \(P(X\leq 6)\) is the so-called **cumulative distribution** function at 6, and R has the function `ppois` for computing it, which we can use in either of the following two ways:[10](https://www.huber.embl.de/msmb/01-chap.html#fn10)

10 Besides the convenience of not having to do the subtraction from one, the second of these computations also tends to be more accurate when the probability is small. This has to do with limitations of floating point arithmetic.

```
1 - ppois(6, 0.5)
```

```
[1] 1.00238e-06
```

```
ppois(6, 0.5, lower.tail = FALSE)
```

```
[1] 1.00238e-06
```

Task

Check the manual page of `ppois` for the meaning of the `lower.tail` argument.

We denote this number by \(\epsilon\), the Greek letter epsilon[11](https://www.huber.embl.de/msmb/01-chap.html#fn11). We have shown that the probability of seeing a count as large as \(7\), assuming no epitope reactions, is:

11 Mathematicians often call small numbers (and children) \(\epsilon\)s.

\[
\epsilon=P(X\geq 7)=1-P(X\leq 6)\simeq10^{-6}.
\]

#### Extreme value analysis for the Poisson distribution

Stop! The above calculation is *not* the correct computation in this case.

Question 1.6

Can you spot the flaw in our reasoning if we want to compute the probability that we observe these data if there is no epitope?

Solution

We looked at all 100 positions, looked for the largest value and found that it was 7. Due to this selection, a value as large as 7 is more likely to occur than if we only looked at one position.

So instead of asking what the chances are of seeing a Poisson(0.5) as large as 7, we should ask ourselves, what are the chances that the maximum of 100 Poisson(0.5) trials is as large as 7? We use **extreme value** analysis here[12](https://www.huber.embl.de/msmb/01-chap.html#fn12). We order the data values \(x\_1,x\_2,... ,x\_{100}\) and rename them \(x\_{(1)},x\_{(2)},x\_{(3)},... ,x\_{(100)}\), so that \(x\_{(1)}\) denotes the smallest and \(x\_{(100)}\) the largest of the counts over the 100 positions. Together, \(x\_{(1)},... x\_{(100)}\) are called the **rank statistic** of this sample of 100 values.

12 Meaning that we’re interested in the behavior of the very large or very small values of a random distribution, for instance the maximum or the minimum.

13 The notation with the \(\prod\) is just a compact way to write the product of a series of terms, analogous to the \(\sum\) for sums.

The maximum value being as large as 7 is the **complementary event** of having all 100 counts be smaller than or equal to 6. Two complementary events have probabilities that sum to 1. Because the positions are supposed to be independent, we can now do the computation[13](https://www.huber.embl.de/msmb/01-chap.html#fn13):

\[
\begin{aligned}
P(x\_{(100)}\geq7)
&=&1-P(x\_{(100)}\leq6)\\
&=&1-P(x\_{(1)}\leq6)\times P(x\_{(2)}\leq6)\times \cdots \times P(x\_{(100)}\leq6)\\
&=&1-P(x\_1\leq6)\times P(x\_2\leq6)\times \cdots \times P(x\_{100}\leq6)\\
&=&1-\prod\_{i=1}^{100} P(x\_i\leq6).
\end{aligned}
\]

Because we suppose each of these 100 events are independent, we can use our result from above:

\[
\prod\_{i=1}^{100} P(x\_i \leq 6)=
\left(P(x\_i \leq 6)\right)^{100}=
\left(1-\epsilon\right)^{100}.
\]

#### Actually computing the numbers

We could just let R compute the value of this number, \(\left(1-\epsilon\right)^{100}\). For those interested in how such calculations can be shortcut through approximation, we give some details. These can be skipped on a first reading.

We recall from above that \(\epsilon\simeq 10^{-6}\) is much smaller than 1. To compute the value of \(\left(1-\epsilon\right)^{100}\) approximately, we can use the binomial theorem and drop all “higher order” terms of \(\epsilon\), i.e., all terms with \(\epsilon^2, \epsilon^3, ...\), because they are negligibly small compared to the remaining (“leading”) terms.

\[
(1-\epsilon)^n = \sum\_{k=0}^n {n\choose k} \, 1^{n-k} \, (-\epsilon)^k =
1-n\epsilon+{n\choose 2} \epsilon^2 - {n\choose 3} \epsilon^3 + ... \simeq
1-n\epsilon \simeq 1 - 10^{-4}
\]

Another, equivalent, route goes by using the approximation \(e^{-\epsilon} \simeq 1-\epsilon\), which is the same as \(\log(1-\epsilon)\simeq -\epsilon\). Hence

\[
(1-\epsilon)^{100} = e^{\log\left((1-\epsilon)^{100}\right)} = e^{ 100 \log (1-\epsilon)} \simeq e^{-100 \epsilon}
\simeq e^{-10^{-4}} \simeq 1 - 10^{-4}.
\]

Thus the correct probability of seeing a number of hits as large or larger than 7 in the 100 positions, if there is no epitope, is about 100 times the probability we wrongly calculated previously.

Both computed probabilities \(10^{-6}\) and \(10^{-4}\) are smaller than standard significance thresholds (say, \(0.05, 0.01\) or \(0.001\)). The decision to reject the null of no epitope would have been the same. However if one has to stand up in court and defend the p-value to 8 significant digits as in some forensic court cases:[14](https://www.huber.embl.de/msmb/01-chap.html#fn14) that is another matter. The adjusted p-value that takes into account the multiplicity of the test is the one that should be reported, and we will return to this important issue in [Chapter 6](https://www.huber.embl.de/msmb/06-chap.html).

14 This occurred in the examination of the forensic evidence in the OJ Simpson case.

#### Computing probabilities by simulation

In the case we just saw, the theoretical probability calculation was quite simple and we could figure out the result by an explicit calculation. In practice, things tend to be more complicated, and we are better to compute our probabilities using the **Monte Carlo** method: a computer simulation based on our generative model that finds the probabilities of the events we’re interested in. Below, we generate 100,000 instances of picking the maximum from 100 Poisson distributed numbers.

```
maxes = replicate(100000, {
  max(rpois(100, 0.5))
})
table(maxes)
```

```
maxes
    1     2     3     4     5     6     7     9 
    7 23028 60840 14364  1604   141    15     1
```

In 16 of 100000 trials, the maximum was 7 or larger. This gives the following approximation for \(P(X\_{\text{max}}\geq 7)\)[15](https://www.huber.embl.de/msmb/01-chap.html#fn15):

15 In R, the expression `maxes >= 7` evaluates into a logical vector of the same length as `maxes`, but with values of `TRUE` and `FALSE`. If we apply the function `mean` to it, that vector is converted into 0s and 1s, and the result of the computation is the fraction of 1s, which is the same as the fraction of `TRUE`s.

```
mean( maxes >= 7 )
```

```
[1] 0.00016
```

which more or less agrees with our theoretical calculation. We already see one of the potential limitations of Monte Carlo simulations: the “granularity” of the simulation result is determined by the inverse of the number of simulations (100000) and so will be around 10^{-5}. Any estimated probability cannot be more precise than this granularity, and indeed the precision of our estimate will be a few multiples of that. Everything we have done up to now is only possible because we know the false positive rate per position, we know the number of patients assayed and the length of the protein, we suppose we have identically distributed independent draws from the model, and there are no unknown parameters. This is an example of **probability or generative modeling**: all the parameters are known and the mathematical theory allows us to work by **deduction** in a **top-down** fashion.

We postulated the Poisson distribution for the noise, pretending we knew all the parameters and were able to conclude through mathematical deduction.

If instead we are in the more realistic situation of knowing the number of patients and the length of the proteins, but don’t know the distribution of the data, then we have to use **statistical modeling**. This approach will be developed in [Chapter 2](https://www.huber.embl.de/msmb/02-chap.html). We will see that if we have only the data to start with, we first need to **fit** a reasonable distribution to describe it. However, before we get to this harder problem, let’s extend our knowledge of discrete distributions to more than binary, success-or-failure outcomes.

## 1.4 Multinomial distributions: the case of DNA

#### More than two outcomes.

When modeling four possible outcomes, as for instance the boxes in [Figure 1.9](https://www.huber.embl.de/msmb/01-chap.html#fig-Boxes4) or when studying counts of the four nucleotides [A,C,G] and [T], we need to extend the [binomial] model.

[![](imgs/BallsinBoxes4.png)](https://www.huber.embl.de/msmb/imgs/BallsinBoxes4.png "Figure 1.9: The boxes represent four outcomes or levels of a discrete categorical variable. The box on the right represents the more likely outcome.")

Figure 1.9: The boxes represent four outcomes or levels of a discrete **categorical** variable. The box on the right represents the more likely outcome.

Recall that when using the binomial, we can consider unequal probabilities for the two outcomes by assigning a probability \(p=P(1)=p\_1\) to the outcome 1 and \(1-p=p(0)=p\_0\) to the outcome 0. When there are more than two possible outcomes, say [A,C,G] and [T], we can think of throwing balls into boxes of differing sizes corresponding to different probabilities, and we can label these probabilities \(p\_A,p\_C,p\_G,p\_T\). Just as in the binomial case the sum of the probabilities of all possible outcomes is 1, \(p\_A+p\_C+p\_G+p\_T=1\).

[![You are secretly meeting a continuous distribution here, the uniform distribution: runif.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "You are secretly meeting a continuous distribution here, the uniform distribution: runif.")

You are secretly meeting a continuous distribution here, the uniform distribution: `runif`.

Task

Experiment with the random number generator that generates all possible numbers between \(0\) and \(1\) through the function called `runif`. Use it to generate a random variable with 4 levels (A, C, G, T) where \(p\_{\text{A}}=\frac{1}{8}, p\_{\text{C}}=\frac{3}{8}, p\_{\text{G}}=\frac{3}{8}, p\_{\text{T}}=\frac{1}{8}\).

**Mathematical formulation.** Multinomial distributions are the most important model for tallying counts and R uses a general formula to compute the probability of a **multinomial** vector of counts \((x\_1,...,x\_m)\) for outcomes of \(n\) draws from \(m\) boxes with probabilities \(p\_1,...,p\_m\):

The first term reads: the joint probability of observing count \(x\_1\) in box 1 and \(x\_2\) in 2 and … \(x\_m\) in box m, given that box 1 has probability \(p\_1\), box 2 has probability \(p\_2\), … and box \(m\) has probability \(p\_m\).

\[\begin{align}
P(x\_1,x\_2,...,x\_m)
&=\frac{n!}{\prod\_{i=1}^m x\_i!} \prod\_{i=1}^m p\_i^{x\_i}\\
&={{n}\choose{x\_1,x\_2,...,x\_m}} \; p\_1^{x\_1}\,p\_2^{x\_2}\cdots p\_m^{x\_m}.
\end{align}\]

The term in brackets is called the multinomial coefficient and is an abbreviation for \[{n\choose x\_1,x\_2,...,x\_m}=\frac{n!}{x\_1!x\_2!\cdots x\_m!}.\] So this is a generalization of the binomial coefficient – for \(m=2\) it is the same as the binomial coefficient.

Question 1.7

Suppose we have four boxes that are equally likely. Using the formula, what is the probability of observing 4 in the first box, 2 in the second box, and none in the two other boxes?

Solution

\[
P(4,2,0,0)=\frac{6\times 5\times 4\times 3\times 2}{4\times 3\times 2 \times 2}
\frac{1}{4^6}
=\frac{15}{4^6}\simeq 0.0037.
\]

```
dmultinom(c(4, 2, 0, 0), prob = rep(1/4, 4))
```

```
[1] 0.003662109
```

We often run simulation experiments to check whether the data we see are consistent with the simplest possible four-box model where each box has the same probability 1/4. In some sense it is the strawman (nothing interesting is happening). We’ll see more examples of this in [Chapter 2](https://www.huber.embl.de/msmb/02-chap.html). Here we use a few R commands to generate such vectors of counts. First suppose we have 8 characters of four different, equally likely types:

```
pvec = rep(1/4, 4)
t(rmultinom(1, prob = pvec, size = 8))
```

```
     [,1] [,2] [,3] [,4]
[1,]    1    3    1    3
```

Question 1.8

Try the code without using the `t()` function; what does `t` stand for?

Question 1.9

How do you interpret the difference between `rmultinom(n = 8, prob = pvec, size = 1)` and `rmultinom(n = 1, prob = pvec, size = 8)`? Hint: remember what we did in Sections [1.3.1](https://www.huber.embl.de/msmb/01-chap.html#sec-generative-bernoulli) and [1.3.2](https://www.huber.embl.de/msmb/01-chap.html#sec-generative-binomial).

### 1.4.1 Simulating for power

Let’s see an example of using Monte Carlo for the **multinomial** in a way which is related to a problem scientists often have to solve when planning their experiments: how big a sample size do I need?

[![](imgs/SampleSize.png)](https://www.huber.embl.de/msmb/imgs/SampleSize.png)

Ask a statistician about sample size, they will always tell you they need more data. The larger the sample size, the more sensitive the results. However lab work is expensive, so there is a tricky cost-benefit tradeoff to be considered. This is such an important problem that we have dedicated a whole chapter to it at the end of the book ([Chapter 13](https://www.huber.embl.de/msmb/13-chap.html)).

The term **power** has a special meaning in statistics. It is the probability of detecting something if it *is* there, also called the **true positive rate**.

Conventionally, experimentalists aim for a power of 80% (or more) when planning experiments. This means that if the same experiment is run many times, about 20% of the time it will fail to yield significant results even though it should.

Let’s call \(H\_0\) the null hypothesis that the DNA data we have collected comes from a *fair* process, where each of the 4 nucleotides is equally likely \((p\_A,p\_C,p\_G,p\_T)=(0.25,0.25,0.25,0.25)\). Null here just means: the baseline, where nothing interesting is going on. It’s the strawman that we are trying to disprove (or “reject”, in the lingo of statisticians), so the null hypothesis should be such that deviations from it are interesting[16](https://www.huber.embl.de/msmb/01-chap.html#fn16).

16 If you know a little biology, you will know that DNA of living organisms rarely follows that null hypothesis – so disproving it may not be all that interesting. Here we proceed with this null hypothesis because it allows us to illustrate the calculations, but it can also serve us as a reminder that the choice of a good null hypothesis (one whose rejection is interesting) requires scientific input.

As you saw by running the R commands for 8 characters and 4 equally likely outcomes, represented by equal-sized boxes, we do not always get 2 in each box. It is impossible to say, from looking at just 8 characters, whether the nucleotides come from a fair process or not.

Let’s determine if, by looking at a sequence of length \(n=20\), we can detect whether the original distribution of nucleotides is fair or whether it comes from some other (“alternative”) process.

We generate 1000 simulations from the null hypothesis using the `rmultinom` function. We display only the first 11 columns to save space.

```
obsunder0 = rmultinom(1000, prob = pvec, size = 20)
dim(obsunder0)
```

```
[1]    4 1000
```

```
obsunder0[, 1:11]
```

```
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
[1,]    6    5    6    8    4    6    2    7    5     4     4
[2,]    6    6    3    7    3    3    8    4    3     3     5
[3,]    3    3    6    2    8    3    5    7    4     7     6
[4,]    5    6    5    3    5    8    5    2    8     6     5
```

[![Notice that the top of every column there is an index of the form [,1][,2]... These are the column indices. The rows are labeled [1,][2,].... The object obsunder0 is not a simple vector as those we have seen before, but an array of numbers in a matrix.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "Notice that the top of every column there is an index of the form [,1][,2]... These are the column indices. The rows are labeled [1,][2,].... The object obsunder0 is not a simple vector as those we have seen before, but an array of numbers in a matrix.")

Notice that the top of every column there is an index of the form `[,1][,2]...` These are the column indices. The rows are labeled `[1,][2,]...`. The object `obsunder0` is not a simple vector as those we have seen before, but an array of numbers in a matrix.

Each column in the matrix `obsunder0` is a simulated instance. You can see that the numbers in the boxes vary a lot: some are as big as 8, whereas the expected value is 5=20/4.

#### Creating a test

Remember: we know these values come from a fair process. Clearly, knowing the process’ expected values isn’t enough. We also need a measure of variability that will enable us to describe how much variability is expected and how much is too much. We use as our measure the following statistic, which is computed as the sum of the squares of the differences between the observed values and the expected values relative to the expected values. Thus, for each instance,

This measure weights each of the square residuals relative to their expected values.

\[
{\tt stat}=\frac{(E\_A-x\_A)^2}{E\_A}+
\frac{(E\_C-x\_C)^2}{E\_C}+
\frac{(E\_G-x\_G)^2}{E\_G}+
\frac{(E\_T-x\_T)^2}{E\_T}
=\sum\_i\frac{(E\_i-x\_i)^2}{E\_i}
\tag{1.1}\]

How much do the first three columns of the generated data differ from what we expect? We get:

```
expected0 = pvec * 20
sum((obsunder0[, 1] - expected0)^2 / expected0)
```

```
[1] 1.2
```

```
sum((obsunder0[, 2] - expected0)^2 / expected0)
```

```
[1] 1.2
```

```
sum((obsunder0[, 3] - expected0)^2 / expected0)
```

```
[1] 1.2
```

The values of the measure can differ- you can look at a few more than 3 columns, we’re going to see how to study all 1,000 of them. To avoid repetitive typing, we encapsulate the formula for `stat`, [Equation 1.1](https://www.huber.embl.de/msmb/01-chap.html#eq-generative-stat), in a function:

```
stat = function(obsvd, exptd = 20 * pvec) {
  sum((obsvd - exptd)^2 / exptd)
}
stat(obsunder0[, 1])
```

```
[1] 1.2
```

To get a more complete picture of this variation, we compute the measure for all 1000 instances and store these values in a vector we call `S0`: it contains values generated under \(H\_0\). We can consider the histogram of the `S0` values shown in [Figure 1.10](https://www.huber.embl.de/msmb/01-chap.html#fig-histS0) an estimate of our **null distribution**.

```
S0 = apply(obsunder0, 2, stat)
summary(S0)
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.000   1.200   2.800   3.126   4.400  17.600
```

```
hist(S0, breaks = 25, col = "lavender", main = "")
```

[![](01-chap_files/figure-html/fig-histS0-1.png)](https://www.huber.embl.de/msmb/01-chap_files/figure-html/fig-histS0-1.png "Figure 1.10: The histogram of simulated values S0 of the statistic stat under the null (fair) distribution provides an approximation of the sampling distribution of the statistic stat.")

Figure 1.10: The histogram of simulated values `S0` of the statistic `stat` under the null (fair) distribution provides an approximation of the **sampling distribution** of the statistic `stat`.

[![The apply function is shorthand for a loop over the rows or columns of an array. Here the second argument, 2, indicates looping over the columns.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "The apply function is shorthand for a loop over the rows or columns of an array. Here the second argument, 2, indicates looping over the columns.")

The `apply` function is shorthand for a loop over the rows or columns of an array. Here the second argument, 2, indicates looping over the columns.

The summary function shows us that `S0` takes on a spread of different values. From the simulated data we can approximate, for instance, the 95% quantile (a value that separates the smaller 95% of the values from the 5% largest values).

```
q95 = quantile(S0, probs = 0.95)
q95
```

```
95% 
7.6
```

So we see that 5% of the `S0` values are larger than 7.6. We’ll propose this as our critical value for testing data and will reject the hypothesis that the data come from a fair process, with equally likely nucleotides, if the weighted sum of squares `stat` is larger than 7.6.

#### Determining our test’s power

We must compute the probability that our test—based on the weighted sum-of-square differences—will detect that the data in fact do not come from the null hypothesis. We compute the probability of rejecting by simulation. We generate 1000 simulated instances from an alternative process, parameterized by `pvecA`.

[![](imgs/roulette.png)](https://www.huber.embl.de/msmb/imgs/roulette.png)

```
pvecA = c(3/8, 1/4, 1/4, 1/8)
observed = rmultinom(1000, prob = pvecA, size = 20)
dim(observed)
```

```
[1]    4 1000
```

```
observed[, 1:7]
```

```
     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
[1,]   10    4    8    8    4    7    7
[2,]    3   10    5    6    6    7    2
[3,]    5    3    5    6    4    2    6
[4,]    2    3    2    0    6    4    5
```

```
apply(observed, 1, mean)
```

```
[1] 7.469 4.974 5.085 2.472
```

```
expectedA = pvecA * 20
expectedA
```

```
[1] 7.5 5.0 5.0 2.5
```

As with the simulation from the null hypothesis, the observed values vary considerably. The question is: how often (out of 1000 instances) will our test detect that the data depart from the null?

The test doesn’t reject the first observation, (10, 3, 5, 2), because the value of the statistic is within the 95th percentile.

```
stat(observed[, 1])
```

```
[1] 7.6
```

```
S1 = apply(observed, 2, stat)
q95
```

```
95% 
7.6
```

```
sum(S1 > q95)
```

```
[1] 199
```

```
power = mean(S1 > q95)
power
```

```
[1] 0.199
```

[![We read the vertical line as given or conditional on.](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "We read the vertical line as given or conditional on.")

We read the vertical line as **given** or **conditional on**.

Run across 1000 simulations, the test identified 199 as coming from an alternative distribution. We have thus computed that the probability \(P(\text{reject }H\_0 \;|\; H\_A)\) is 0.199.

With a sequence length of \(n = 20\) we have a *power* of about 20% to detect the difference between the fair generating process and our **alternative**.

Task

In practice, as we mentioned, an acceptable value of power is \(0.8\) or more. Repeat the simulation experiments and suggest a new sequence length \(n\) that will ensure that the power is acceptable.

#### Classical statistics for classical data

We didn’t need to simulate the data using Monte Carlo to compute the 95th percentiles; there is an adequate theory to help us with the computations.

Our statistic `stat` actually has a well-known distribution called the chi-square distribution (with 3 degrees of freedom) and written \({\chi}^2\_3\).

[![](imgs/ProbaDiagram.png)](https://www.huber.embl.de/msmb/imgs/ProbaDiagram.png "Figure 1.11: We have studied how a probability model has a distribution, which we call F. F often depends on parameters, which are–by convention–denoted by Greek letters, such as \theta. The observed data are generated via the brown arrow and are represented by Latin letters, such as x. The vertical bar in the probability computation stands for supposing that or conditional on")

Figure 1.11: We have studied how a probability model has a distribution, which we call \(F\). \(F\) often depends on parameters, which are–by convention–denoted by Greek letters, such as \(\theta\). The observed data are generated via the brown arrow and are represented by Latin letters, such as \(x\). The vertical bar in the probability computation stands for **supposing that** or **conditional on**

We will see in [Chapter 2](https://www.huber.embl.de/msmb/02-chap.html) how to compare distributions using Q-Q plots (see [Figure 2.8](https://www.huber.embl.de/msmb/02-chap.html#fig-qqplot3-1)). We could have used a more standard test instead of running a hand-made simulation. However, the procedure we’ve learned extends to many situations in which the chi-square distribution doesn’t apply. For instance, when some of the boxes have extremely low probabilities and their counts are mostly zero.

## 1.5 Summary of this chapter

We have used mathematical formulæ and R to compute probabilities of various discrete *events*, using a few basic distributions:

**The *Bernoulli* distribution** was our most basic building block – it is used to represent a single binary trial such as a coin flip. We can code the outcomes as 0 and 1. We call \(p\) the probability of success (the 1 outcome).

**The *binomial* distribution** is used for the number of 1s in \(n\) binary trials, and we can compute the probabilities of seeing \(k\) successes using the R function `dbinom`. We also saw how to simulate a binomial eperiment with \(n\) trials using the function `rbinom`.

**The *Poisson* distribution** is most appropriate for cases when \(p\) is small (the 1s are rare). It has only one parameter \(\lambda\), and the Poisson distribution for \(\lambda=np\) is approximately the same as the binomial distribution for \((n,p)\) if \(p\) is small. We used the Poisson distribution to model the number of randomly occurring false positives in an assay that tested for epitopes along a sequence, presuming that the per-position false positive rate \(p\) was small. We saw how such a parametric model enabled us to compute the probabilities of extreme events, as long as we knew all the parameters.

**The *multinomial* distribution** is used for discrete events that have more than two possible outcomes or **levels**. The power example showed us how to use Monte Carlo simulations to decide how much data we need to collect if we want to test whether a multinomial model with equal probabilities is consistent with the data. We used probability distributions and probabilistic models to evaluate hypotheses about how our data were generated, by making assumptions about the generative models. We term the probability of seeing the data, given a hypothesis, a **p-value**. This is not the same as the probability that the hypothesis is true!

[![P(H_0\;|\;\text{data}) is not the same as a p-value P(\text{data}\;|\;H_0).](imgs/devil.png)](https://www.huber.embl.de/msmb/imgs/devil.png "P(H_0\;|\;\text{data}) is not the same as a p-value P(\text{data}\;|\;H_0).")

\(P(H\_0\;|\;\text{data})\) is not the same as a p-value \(P(\text{data}\;|\;H\_0)\).

## 1.6 Further reading

* The elementary book by Freedman, Pisani, and Purves ([1997](https://www.huber.embl.de/msmb/16-chap.html#ref-Freedman:1997)) provides the best introduction to probability through the type of box models we mention here.
* The book by Durbin et al. ([1998](https://www.huber.embl.de/msmb/16-chap.html#ref-DEKM)) covers many useful probability distributions and provides in its appendices a more complete view of the theoretical background in probability theory and its applications to sequences in biology.
* Monte Carlo methods are used extensively in modern statistics. Robert and Casella ([2009](https://www.huber.embl.de/msmb/16-chap.html#ref-Casella2009)) provides an introduction to these methods using R.
* [Chapter 6](https://www.huber.embl.de/msmb/06-chap.html) will cover the subject of hypothesis testing. We also suggest Rice ([2006](https://www.huber.embl.de/msmb/16-chap.html#ref-Rice:2007)) for more advanced material useful for the type of more advanced probability distributions, beta, gamma, exponentials we often use in data analyses.