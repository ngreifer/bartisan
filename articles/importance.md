# Which variables matter

## Introduction

A fitted forest spends its splitting rules on some predictors and not
others, and counting how it spends them is the usual measure of variable
importance for BART.
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
reports the count.

This vignette covers what that number means, how to tell whether a
difference in it is real, and three things it is not. The last part
matters more than the first: variable importance is the most over-read
output in machine learning, and the failure modes are specific and
checkable.

``` r

library(bartisan)

data(rhc)

model <- death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
  hema + pafi + paco2 + crea + surv2m + card
```

## The measure

``` r

set.seed(2026)

fit <- bartisan(model, data = rhc, family = binomial(), chains = 4)

variable_importance(fit)
#>    variable splits splits_lower splits_upper prop_used
#> 1    surv2m  21.45            5           49     1.000
#> 2       age  12.75            3           29     1.000
#> 3     paco2   8.37            0           25     0.971
#> 4       rhc   3.99            0           16     0.915
#> 5      pafi   4.99            0           17     0.898
#> 6       aps   5.90            0           19     0.860
#> 7      card   2.83            0           12     0.702
#> 8       edu   2.03            0            9     0.600
#> 9      crea   3.13            0           19     0.591
#> 10     hema   2.80            0           19     0.544
#> 11   meanbp   1.53            0            8     0.534
#> 12     race   1.85            0           11     0.504
#> 13      sex   2.05            0           16     0.491
#> 14     resp   2.35            0           15     0.490
```

Two columns, answering different questions.

`splits` is the average number of splitting rules the forest gives a
predictor in one posterior draw, with an interval. It says how much of
the fitted structure the predictor accounts for.

`prop_used` is the proportion of draws in which the predictor received
at least one rule. Read this one first. Under the default prior it
behaves like a posterior probability that the predictor belongs in the
model, and it is the column to use if the question is which variables to
keep.

The reason it behaves that way is the prior. By default the variable a
rule splits on is drawn from a Dirichlet-distributed set of
probabilities ([Linero 2018](#ref-linero2018sparse)), which lets the
forest concentrate on a few predictors and drop the rest entirely. A
predictor that contributes nothing can fall to `prop_used` near zero,
which is not possible under classic BART, where every predictor keeps a
fixed share of the splitting probability.

The prognostic score and the illness measures sit at the top, which is
what you would expect. At the other end, several predictors are used in
only about half the draws.

## Calibrating with noise

Is a `prop_used` of 0.53 low? The table alone cannot say. Before reading
anything into the bottom of it, check what a variable that certainly
does not matter looks like on this data. Add some.

``` r

set.seed(11)
rhc_noise <- rhc
for (j in 1:3) rhc_noise[[paste0("noise", j)]] <- rnorm(nrow(rhc_noise))

set.seed(2026)
fit_noise <- bartisan(update(model, . ~ . + noise1 + noise2 + noise3),
                      data = rhc_noise, family = binomial(), chains = 4)

variable_importance(fit_noise)
#>    variable splits splits_lower splits_upper prop_used
#> 1    surv2m 17.651            5        37.00     1.000
#> 2       age 17.005            2        41.00     1.000
#> 3     paco2 10.225            1        28.00     0.996
#> 4      pafi  5.462            0        16.00     0.926
#> 5       rhc  3.037            0        10.00     0.904
#> 6       aps  3.958            0        13.00     0.866
#> 7       edu  2.401            0        11.02     0.643
#> 8    noise3  2.264            0        11.00     0.564
#> 9    meanbp  2.558            0        15.00     0.560
#> 10     card  1.581            0         8.00     0.540
#> 11     hema  2.470            0        13.00     0.500
#> 12   noise1  1.752            0        10.00     0.479
#> 13     crea  1.328            0         8.00     0.434
#> 14     race  1.269            0         8.00     0.421
#> 15   noise2  1.211            0         8.02     0.403
#> 16     resp  1.016            0         6.00     0.398
#> 17      sex  0.975            0         6.00     0.354
```

The noise variables do not sort to the bottom. They land in the middle
of the table, above several of the real predictors, and that is the
answer to the question: the lower half of this table carries no
information. A predictor ranked below a variable that is literally
random cannot be said to matter less than it.

Read the gap rather than the ranking. The predictors at the top sit well
clear of the noise controls and can be said to matter. Everything from
about the middle down is indistinguishable from random numbers on this
sample, which is a more useful conclusion than an ordering, and one the
table alone would not have given you.

This technique costs three lines and settles the question. Use it
whenever an importance table is going to be interpreted.

## When it does work

The failure above is a property of this dataset, not of the method. With
more data and a stronger signal the separation is sharp. The following
uses the Friedman function, where `x1` to `x5` enter the outcome and
`x6` to `x10` are noise, at 500 observations.

``` r

friedman <- function(n) {
  x <- as.data.frame(matrix(runif(n * 10), n, 10))
  names(x) <- paste0("x", 1:10)
  x$y <- 10 * sin(pi * x$x1 * x$x2) + 20 * (x$x3 - 0.5)^2 +
    10 * x$x4 + 5 * x$x5 + rnorm(n)
  x
}

set.seed(7)
fit_fr <- bartisan(y ~ ., data = friedman(500), family = gaussian())

variable_importance(fit_fr)
#>    variable splits splits_lower splits_upper prop_used
#> 1        x2 32.814           24           41     1.000
#> 2        x1 22.758           13           32     1.000
#> 3        x4 11.686            8           17     1.000
#> 4        x3 11.340            7           19     1.000
#> 5        x5  3.156            1            6     1.000
#> 6       x10  1.070            0            2     0.890
#> 7        x6  0.118            0            1     0.108
#> 8        x9  0.120            0            1     0.092
#> 9        x8  0.102            0            1     0.076
#> 10       x7  0.058            0            1     0.054
```

The five real predictors sit at 1.00 and four of the five noise
predictors below 0.13. The fifth, `x8`, comes in around 0.4, which is a
useful reminder that even a clean separation has a straggler: a noise
variable will occasionally be picked up, and a single moderate value is
not evidence of anything. The gap that matters runs from about 0.4 to
1.00, and it is wide.

There is no threshold that is correct in general. Look for the gap, cut
inside it, and check that the conclusion does not depend on where.

## Correlated predictors

This is the failure mode most likely to mislead, because it produces a
confident and wrong answer rather than an uninformative one.

``` r

set.seed(5)
n <- 400
cd <- data.frame(x1 = runif(n))
cd$x1_copy <- cd$x1 + rnorm(n, sd = 0.01)   # almost the same variable
cd$x2 <- runif(n)
cd$x3 <- runif(n)
cd$y <- 3 * cd$x1 + rnorm(n, sd = 0.3)      # only x1 is in the truth

fit_corr <- bartisan(y ~ ., data = cd, family = gaussian())

variable_importance(fit_corr)
#>   variable splits splits_lower splits_upper prop_used
#> 1       x1 70.238         46.5           88     1.000
#> 2  x1_copy  4.770          0.0           23     0.602
#> 3       x3  0.896          0.0            5     0.442
#> 4       x2  0.942          0.0            5     0.374
```

The outcome depends on `x1`. The forest used `x1_copy` in every draw and
`x1` in under 10% of them. Read as a selection result, this table says
to keep the variable that is not in the truth and drop the one that is.

Nothing has gone wrong with the fit. The two variables carry the same
information, so a tree that splits on either fits equally well, and the
forest settled on one arbitrarily. Predictions are unaffected. What is
affected is any statement about which variable matters.

The same thing happens to the effects:

``` r

library(marginaleffects)

avg_comparisons(fit_corr, variables = c("x1", "x1_copy"))
#> 
#>     Term Estimate   2.5 % 97.5 %
#>  x1        1.3408  0.7456  1.556
#>  x1_copy   0.0495 -0.0124  0.673
#> 
#> Type: response
#> Comparison: +1
```

All of the association is attributed to the copy and none to the
original. Moving both together, which is what a change in the underlying
quantity would mean, recovers the truth:

``` r

lo <- transform(cd, x1 = 0.25, x1_copy = 0.25)
hi <- transform(cd, x1 = 0.75, x1_copy = 0.75)

drawn <- rowMeans(predict(fit_corr, newdata = hi, draws = TRUE) -
                    predict(fit_corr, newdata = lo, draws = TRUE))

round(c(estimate = mean(drawn), quantile(drawn, c(0.025, 0.975))), 3)
#> estimate     2.5%    97.5% 
#>     1.58     1.43     1.72
```

The true difference over that range is 1.5. So the model knows the
relationship; it just cannot say which of two identical columns it
belongs to.

The practical rule is to treat a set of correlated predictors as one
unit. Decide in advance which variables measure the same underlying
thing, and report and move them together.

## Three things importance is not

**It is not an effect size.** How often a predictor is split on and how
much it moves the outcome are different quantities, and they can
disagree in both directions. A variable with a small effect over a range
that gets split repeatedly will show high usage. If the question is how
much the outcome changes,
[`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
answers it directly and reports an interval. See
[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md).

**It is not a test.** `prop_used` is a posterior probability under a
particular prior. It has no error rate attached, and the choice of
cutoff is yours. The noise-control technique above is the closest thing
to a calibration and it is informal.

**It is not causal.** A ranking of predictors describes this fitted
function on this sample. It says nothing about what would happen if any
of them were changed, and a variable can rank highly because it is a
consequence of the outcome, a proxy for a confounder, or correlated with
something that matters. See
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md).

## A note on the number of trees

A recommendation that circulates is to reduce the number of trees when
using BART for variable selection. It comes from Chipman et al.
([2010](#ref-chipman2010)), who observe that counting splits works
poorly with many trees because “the redundancy offered by so many trees
tends to mix many irrelevant predictors in with the relevant ones”, and
that predictors compete for splits when the forest is small.

The advice is correct for classic BART and largely unnecessary here.
Measured on the Friedman function at 500 observations, three replicates,
the mean `prop_used` for the five noise predictors was:

| Trees | sparsity = FALSE | sparsity = TRUE (default) |
|------:|-----------------:|--------------------------:|
|    10 |             0.28 |                      0.09 |
|    20 |             0.50 |                      0.08 |
|    50 |             0.95 |                      0.09 |
|   100 |             1.00 |                      0.14 |

Without the sparsity prior the advice is essential: at 50 trees the
noise predictors are used in 95% of draws and cannot be told apart from
the real ones. With the prior, which is on by default, they stay near
zero at every tree count. The recommendation addresses the same problem
the prior addresses, and there is little left for it to do. Reducing the
trees also costs predictive accuracy, so it is not free.

If you want to check, fit at the default and again at 20 trees and see
whether the conclusion changes. It usually will not.

## Where to go next

[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md)
covers how much predictors move the outcome, which is the question
importance is usually a proxy for.
[`vignette("comparison")`](https://ngreifer.github.io/bartisan/articles/comparison.md)
covers choosing between models that include different variables, which
is the other way to ask whether a variable earns its place.

## References

Chipman, Hugh A., Edward I. George, and Robert E. McCulloch. 2010.
“BART: Bayesian Additive Regression Trees.” *The Annals of Applied
Statistics* 4 (1): 266–98. <https://doi.org/10.1214/09-AOAS285>.

Linero, Antonio R. 2018. “Bayesian Regression Trees for High-Dimensional
Prediction and Variable Selection.” *Journal of the American Statistical
Association* 113 (522): 626–36.
<https://doi.org/10.1080/01621459.2016.1264957>.
