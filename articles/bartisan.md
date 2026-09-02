# Getting started with bartisan

## Introduction

*bartisan* fits Bayesian additive regression trees (BART) using the same
interface as [`glm()`](https://rdrr.io/r/stats/glm.html). You supply a
formula, a data frame, and a family, and the model estimates the
relationship between the predictors and the outcome without you having
to say what shape that relationship takes. Nonlinearity and interactions
are found rather than specified.

This vignette walks through a complete analysis: fitting a model,
checking that it worked, seeing which predictors it uses, reading off
the effects, predicting for new observations, and comparing models. It
assumes you are comfortable with regression but does not assume
familiarity with machine learning or Bayesian methods.

The main thing to take from it is that the defaults are meant to be
used. The priors, the number of trees, and the sampler settings are
chosen to work across a wide range of problems, and tuning them is
rarely where the gains are. Almost everything below is a single function
call with no arguments beyond the formula and the data.

Each section ends with a pointer to a vignette that covers the same
ground in more depth.

``` r

library(bartisan)
```

## The data

`rhc` records 1500 critically ill patients from the SUPPORT study and
whether each received right heart catheterization, a monitoring
procedure, within a day of arriving in intensive care ([Connors et al.
1996](#ref-connors1996)). The question the study asked is whether the
procedure helps or harms.

``` r

data(rhc)

str(rhc)
#> 'data.frame':    1500 obs. of  16 variables:
#>  $ rhc   : int  0 0 0 0 0 1 0 0 0 1 ...
#>  $ death : int  1 0 0 1 0 0 1 1 0 1 ...
#>  $ days  : int  37 235 189 13 202 203 663 15 239 22 ...
#>  $ age   : num  75.3 55 34.4 42.2 41.4 ...
#>  $ sex   : Factor w/ 2 levels "female","male": 1 2 2 1 2 2 2 2 2 2 ...
#>  $ race  : Factor w/ 3 levels "white","black",..: 1 1 1 1 2 1 1 1 1 1 ...
#>  $ edu   : num  9 14 15 16 11 ...
#>  $ aps   : int  48 29 21 55 60 68 26 89 59 105 ...
#>  $ meanbp: num  55 67 66 77 53 47 63 44 50 33 ...
#>  $ resp  : num  26 10 30 40 12 40 22 0 33 44 ...
#>  $ hema  : num  26.3 29 23.8 53 31 ...
#>  $ pafi  : num  157 149 202 171 390 ...
#>  $ paco2 : num  30 45 37 25 31 28 40 39 36 31 ...
#>  $ crea  : num  1.7 1 0.5 1.9 15 ...
#>  $ surv2m: num  0.441 0.339 0.846 0.672 0.777 ...
#>  $ card  : Factor w/ 2 levels "no","yes": 1 1 1 1 1 1 2 2 1 1 ...
```

The outcome comes in two forms. `death` is whether the patient died
during follow-up, and `days` is how long that took. This vignette uses
the binary form;
[`vignette("survival")`](https://ngreifer.github.io/bartisan/articles/survival.md)
uses the other.

The remaining variables are the patient’s age, sex, race and years of
education, seven physiological measurements taken on the first day,
whether cardiovascular disease was among the diagnoses, and `surv2m`,
the study’s own estimate of the patient’s chance of surviving two
months. All were recorded before catheterization.

## Fitting the model

``` r

set.seed(2026)

# For parallelization; optional
future::plan(future::multisession)

fit <- bartisan(
  death ~ rhc + age + sex + race + edu + aps + meanbp + resp + hema +
    pafi + paco2 + crea + surv2m + card,
  data = rhc, family = binomial(), chains = 4
)

fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ rhc + age + sex + race + edu + aps + 
#>     meanbp + resp + hema + pafi + paco2 + crea + surv2m + card, 
#>     data = rhc, family = binomial(), chains = 4)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 1 forest of 50 trees, soft decision rules
#> Draws: 2000 kept across 4 chains after 500 warmup
```

That is the whole call. `family = binomial()` says the outcome is
binary, exactly as in [`glm()`](https://rdrr.io/r/stats/glm.html). The
family may be omitted, in which case it is read off the outcome and
reported; naming it is clearer and silences the message.

`chains = 4` runs the sampler four times from different starting points.
The default is one chain, which is enough to get estimates, but running
several is what makes the convergence diagnostics in the next section
available.

Choosing a family is the one modeling decision that usually matters more
than any sampler setting.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers the choice, and
[`vignette("survival")`](https://ngreifer.github.io/bartisan/articles/survival.md)
covers censored outcomes such as this one’s `days`.

## Checking the model

Two questions are worth separating: whether the sampler converged, and
whether the model fits.

### Convergence

[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
is the one call to make. It computes the convergence and mixing
statistics, applies the conventional thresholds, and says what to change
about whichever of them fall short.

``` r

diagnose(fit)
#> Convergence and mixing
#> 
#>                            quantity  rhat rhat_late ess_bulk ess_tail
#>                              loglik 1.178     1.319       16       76
#>                          splits.eta 1.058     1.055       62      258
#>  eta.eta (worst 5% of observations) 1.086     1.145       34      137
#> 
#> What to do
```

`rhat` compares variation between chains to variation within them;
values near 1 indicate the chains have settled on the same answer.
`rhat_late` is the same statistic on the second half of the draws alone,
which is what distinguishes a warmup that ended too early from chains
that have each settled somewhere different. `ess_bulk` and `ess_tail`
are effective sample sizes, and count how many independent draws the
correlated ones are worth, in the middle of the distribution and in the
tails.

Read the rows that correspond to quantities you will report. `eta.eta`
is the fitted function, summarized over its worst 5% of observations, so
it is the row that matters here. Forests mix slowly on their fitted
values, so a figure above 1.01 on that row is ordinary rather than
alarming;
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
explains what to do about it and when to worry.

The table is also in `fit$rhat` for a fit with more than one chain, if
what you want is the numbers rather than the report.

### Fit

A posterior predictive check simulates new outcomes from the fitted
model and compares their distribution to the observed one.

``` r

bayesplot::pp_check(fit)
```

For a binary outcome this is a weaker check than it is for a continuous
one: there are only two values to get right, so the replicates match
unless something has gone badly wrong.
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
covers what to do instead.

## Which predictors the model uses

``` r

variable_importance(fit)
#>    variable splits splits_lower splits_upper prop_used
#> 1    surv2m 21.455            5           49    1.0000
#> 2       age 12.748            3           29    1.0000
#> 3     paco2  8.373            0           25    0.9705
#> 4       rhc  3.987            0           16    0.9155
#> 5      pafi  4.992            0           17    0.8985
#> 6       aps  5.897            0           19    0.8595
#> 7      card  2.828            0           12    0.7015
#> 8       edu  2.030            0            9    0.6000
#> 9      crea  3.126            0           19    0.5905
#> 10     hema  2.802            0           19    0.5440
#> 11   meanbp  1.530            0            8    0.5340
#> 12     race  1.846            0           11    0.5040
#> 13      sex  2.053            0           16    0.4910
#> 14     resp  2.349            0           15    0.4900
```

`splits` is the average number of splitting rules the forest spends on
each predictor per draw, and `prop_used` is the proportion of draws in
which the predictor received any rule at all.

`surv2m` takes the most rules, which is unsurprising: it is a prognostic
score built to predict survival. Age follows. At the bottom, race and
sex are used in about half the draws, which says the model can often do
without them.

Two cautions. Usage is not effect size: a predictor can be split on
constantly and still move the prediction very little, and
[`marginaleffects::avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
in the next section is the better guide to that. And, when predictors
are correlated, the usage distributes among them more or less
arbitrarily.

[`vignette("importance")`](https://ngreifer.github.io/bartisan/articles/importance.md)
covers variable importance and selection, including how to tell whether
a difference in this table means anything.

## Interpreting the fit

A forest has no coefficients, so there is no table of slopes to read.
The question “what is the effect of catheterization” is answered by
asking the fitted model what it predicts when every patient receives it,
asking again when none does, and taking the difference.
*marginaleffects* does this.

``` r

library(marginaleffects)

avg_comparisons(fit, variables = "rhc")
#> 
#>  Estimate 2.5 % 97.5 %
#>    0.0578     0  0.107
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
```

Catheterization is associated with an increase of about six percentage
points in the probability of death. The interval runs from roughly zero
to eleven points, so the direction is reasonably clear and the size is
not.

The lower bound is exactly zero rather than merely close to it, and that
is worth knowing about. The default splitting prior can drop a predictor
from the forest entirely, and in a draw where it drops `rhc` the
contrast is exactly zero, so the posterior has a point mass there. The
default settings are not necessarily the best ones to use for causal
effect estimation; more specialized methods, like Bayesian causal
forests (BCF) and BART without sparsity-inducing priors, are described
at
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md).

Because the outcome is binary, this is a difference in probability,
which is interpretable without reference to the model. That is usually
the number to report. The underlying probabilities are also worth
showing:

``` r

avg_predictions(fit, variables = "rhc")
#> 
#>  rhc Estimate 2.5 % 97.5 %
#>    0    0.633 0.604  0.663
#>    1    0.690 0.644  0.729
#> 
#> Type: response
```

About 63% of patients would be expected to die without catheterization
and 69% with it, averaging over the covariates as they actually occur in
this sample.

### Looking at a relationship

Effects averaged over the sample hide the shape of the relationship. To
see the shape, plot the model’s predictions against one predictor.

``` r

plot_predictions(fit, condition = "surv2m") +
  ggplot2::labs(x = "estimated probability of surviving two months",
                y = "fitted probability of death") +
  ggplot2::theme_bw()
```

![](bartisan_files/figure-html/pdp-1.png)

The fitted probability of death falls from about 0.88 to about 0.45 as
the prognostic score rises, and the fall is not a straight line. A
logistic regression reports one slope on the log-odds scale for the
whole range. Nothing had to be specified to find the shape.

The band is a credible interval, and it widens at the top where few
patients were that healthy.

[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md)
covers effects, curves, and interactions.

## Predicting new observations

[`predict()`](https://rdrr.io/r/stats/predict.html) returns the
posterior mean prediction.

``` r

new_patient <- rhc[1, ]
new_patient$rhc <- 1

predict(fit, newdata = new_patient)
#> [1] 0.8301
```

For a prediction with an interval, use
[`marginaleffects::predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html)[^1]:

``` r

predictions(fit, newdata = new_patient)
#> 
#>  Estimate 2.5 % 97.5 %
#>     0.835 0.729  0.901
#> 
#> Type: response
```

This describes the probability that a patient with these characteristics
dies. It is an interval for that probability, not a statement about
which way any individual patient will go: the outcome itself is either 0
or 1, and a probability of 0.7 is entirely compatible with survival.

For a continuous outcome, the distinction between an interval for the
mean and an interval for a new observation matters a great deal, and
[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md)
covers it. For a binary outcome the second is rarely what anyone wants.

## Comparing models

Approximate leave-one-out cross-validation estimates how well a model
predicts data it has not seen.

``` r

library(loo)

loo(fit)
#> 
#> Computed from 2000 by 1500 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -849.4 17.4
#> p_loo        33.5  1.0
#> looic      1698.8 34.8
#> ------
#> MCSE of elpd_loo is 0.6.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.0, 0.3]).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
```

`elpd_loo` is the estimated log predictive density on held-out data,
where higher is better. `p_loo` is the effective number of parameters,
which is a measure of how much of the data the forest is actually using.
The Pareto k diagnostics are all good, meaning the approximation is
trustworthy for this fit.

Two models can be compared directly. Here we ask whether the
physiological measurements earn their keep over knowing the patient’s
demographics alone:

``` r

set.seed(2026)
demographics <- bartisan(death ~ rhc + age + sex + race + edu, data = rhc,
                         family = binomial(), chains = 4)

loo_compare(list(full = loo(fit), demographics = loo(demographics)))
#>         model elpd_diff se_diff p_worse diag_diff diag_elpd
#>          full       0.0     0.0      NA                    
#>  demographics     -67.4    11.2    1.00
```

The full model predicts better by around six times the standard error of
the difference, which is what you would expect: how sick a patient is on
arrival is the main thing that predicts whether they die.

[`vignette("comparison")`](https://ngreifer.github.io/bartisan/articles/comparison.md)
covers model comparison, and the cases where leave-one-out fails.

## What to be careful about

The model is flexible about the shape of the relationship and nothing
else.

It does not make an association causal. Patients were not randomized to
catheterization; sicker patients were more likely to receive it, which
is exactly the kind of confounding that can produce an apparent harm.
Whether the estimate above can be read as the effect of the procedure
depends on whether the covariates account for that selection, which is a
question about the study and not about the fit.
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md)
covers what is required, using this same data.

It does not extrapolate reliably. Predictions for predictor values
outside the range of the training data are shrunk toward the overall
mean rather than continuing any trend.

It does not fix a badly chosen family. Getting the outcome distribution
wrong matters more than any sampler setting.

## Where to go next

| Topic | Vignette |
|----|----|
| How BART works, and the sampler | [`vignette("implementation")`](https://ngreifer.github.io/bartisan/articles/implementation.md) |
| Choosing a family | [`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md) |
| Convergence and fit | [`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md) |
| Variable importance and selection | [`vignette("importance")`](https://ngreifer.github.io/bartisan/articles/importance.md) |
| Effects, curves, and interactions | [`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md) |
| Model comparison | [`vignette("comparison")`](https://ngreifer.github.io/bartisan/articles/comparison.md) |
| Causal inference | [`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md) |
| Censored and survival outcomes | [`vignette("survival")`](https://ngreifer.github.io/bartisan/articles/survival.md) |

`?bartisan-package` has a shorter version of the same map, organized by
task.

## References

Connors, Alfred F., Theodore Speroff, Neal V. Dawson, et al. 1996. “The
Effectiveness of Right Heart Catheterization in the Initial Care of
Critically Ill Patients.” *JAMA* 276 (11): 889–97.
<https://doi.org/10.1001/jama.1996.03540110043030>.

[^1]: Note that by default, *marginaleffects* uses the posterior median
    as the point estimate, whereas
    [`predict()`](https://rdrr.io/r/stats/predict.html) uses the
    posterior mean, so these values may differ slightly. Use
    `options("marginaleffects_posterior_center" = "mean")` prior to
    running
    [`predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html)
    to produce the posterior mean.
