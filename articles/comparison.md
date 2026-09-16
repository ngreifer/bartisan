# Choosing Between Models

## Introduction

Model selection means something different for a BART model than it does
for a linear model, and the difference is worth stating before any code.

With [`lm()`](https://rdrr.io/r/stats/lm.html) we choose which terms
enter, whether to add a squared term, and whether to include an
interaction; those choices are the model. With BART the forest makes
them, so they are not ours to make. What remains is a shorter list:

1.  Which variables the model is allowed to see.
2.  Which likelihood, meaning which family and link.
3.  Occasionally, a sampler setting (e.g., the number of trees).

In this guide we will compare those choices using leave-one-out
cross-validation. First we’ll fit a model to the `rhc` data and read its
[`loo()`](https://mc-stan.org/loo/reference/loo.html) output, including
the diagnostics that say when the approximation cannot be trusted, and
the two forms of cross-validation to fall back on when it cannot. Next
we’ll compare two sets of variables and then three links, work through
the two ways two log densities end up on different scales and what to do
about each, and put the forest up against a Bayesian logistic
regression. Finally we’ll cover how to tune a setting and how to select
variables, both of which have a right way and a more tempting wrong one.

``` r

library(bartisan)
library(loo)

data("rhc")

model <- death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
  hema + pafi + paco2 + crea + surv2m + card

set.seed(2026)

full <- bartisan(model, data = rhc, family = binomial())
```

The fits here use the default single chain rather than the four used
elsewhere; leave-one-out needs draws rather than chains, and four times
as many fits would make this vignette slow to build for no gain.
Convergence should still be checked separately, as in
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md),
before we trust any of these comparisons.

## Leave-One-Out Cross-Validation (`loo()`)

``` r

loo(full)
#> 
#> Computed from 800 by 1500 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -847.8 17.3
#> p_loo        32.9  0.9
#> looic      1695.5 34.6
#> ------
#> MCSE of elpd_loo is 0.7.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.0, 0.4]).
#> 
#> All Pareto k estimates are good (k < 0.66).
#> See help('pareto-k-diagnostic') for details.
```

`elpd_loo` estimates the log predictive density on data the model has
not seen, with higher values being better. It is computed by importance
sampling from the fitted posterior rather than by refitting, which is
why it is fast ([Vehtari et al. 2017](#ref-vehtari2017)).

`p_loo` is the effective number of parameters, a little over thirty
here; for a forest with hundreds of leaves across its trees that number
is small, because the prior shrinks most of them toward zero. It is a
useful measure of how much of the data the model is actually using.

The Pareto \\k\\ diagnostics are the thing to check: leave-one-out by
importance sampling is trustworthy only when the importance weights are
well behaved, and a \\k\\ above the threshold *loo* prints with them,
which is at most .7 and a little lower for a fit with this many draws,
says that for that observation they are not. Here they are all good.

### When the Approximation Fails

A forest is a flexible function, so a single observation can have a good
deal of influence on the leaves it falls into, and the worry is
therefore that high \\k\\ values will be more common than they are for a
parametric model. Measured, they usually are not: the leaf prior shrinks
every leaf toward zero and the fit is a sum over many trees, so no
single observation dominates the leaves it reaches. The exceptions that
do turn up are usually about the likelihood rather than the trees, which
makes the warning worth reading rather than expecting. If many
observations are flagged, the estimate is unreliable, and the remedy is
held-out data rather than a different diagnostic: fit to part of the
sample and score the part the model was not shown.

``` r

set.seed(2026)
train_id <- sample.int(nrow(rhc), 1200)

train <- rhc[train_id, ]
held  <- rhc[-train_id, ]

fit_train <- bartisan(model, data = train, family = binomial())

score <- predict(fit_train, newdata = held, type = "density", log = TRUE)

sum(score)
#> [1] -172.6
```

`type = "density"` evaluates the outcome under each posterior draw and
averages over the draws before taking the log, which is the same form as
one pointwise `elpd_loo` contribution. The two are therefore estimates
of the same thing, which the two numbers do not make obvious:

``` r

elpd_total <- loo(full)$estimates["elpd_loo", "Estimate"]

rbind(loo     = c(total = elpd_total, n = nrow(rhc),     per_obs = elpd_total / nrow(rhc)),
      heldout = c(total = sum(score), n = length(score), per_obs = mean(score)))
#>          total    n per_obs
#> loo     -847.8 1500 -0.5652
#> heldout -172.6  300 -0.5753
```

The totals differ by a factor of five because they sum different numbers
of terms: [`loo()`](https://mc-stan.org/loo/reference/loo.html) scores
every observation and the split scores only the ones held out. It is in
the per-observation column that the two can be read against each other,
and there they nearly agree. The gap that remains is in the direction to
expect, since leave-one-out trains on 1499 observations and the split
trains on 1200.

That is the whole of what a held-out score gives us on its own, which is
not much: a log score has no scale, and the only thing to do with one is
put it beside another. The next section does that, first through
[`loo()`](https://mc-stan.org/loo/reference/loo.html) and then through
this split.

## Comparing Two Models (`loo_compare()`)

We can compare the fit of two models by supplying their
[`loo()`](https://mc-stan.org/loo/reference/loo.html) output to
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).
Below, we fit a model that only includes demographic variables as
predictors to compare to our full model.

``` r

set.seed(2026)
demographics <- bartisan(death ~ rhc + age + sex + race + edu,
                         data = rhc, family = binomial())

loo_compare(list(full = loo(full),
                 demographics = loo(demographics)))
#>         model elpd_diff se_diff p_worse diag_diff diag_elpd
#>          full       0.0     0.0      NA                    
#>  demographics     -69.4    11.1    1.00
```

The full model predicts better by around six times the standard error of
the difference. The standard error is the important half: a difference
of many standard errors is clear, and a difference smaller than its own
standard error is not evidence of anything.

`p_worse` reports that reading as a probability. It is
`pnorm(0, elpd_diff, se_diff)`, the chance that the model on that row is
really the worse of the two given how far apart they came out and how
precisely the difference is known. The best-ranked model has nothing to
be compared against and gets `NA`. Two properties are worth knowing
before reading one. Because the models are sorted by `elpd_loo` before
it is computed, every reported value is at least .5 by construction, so
.5 does not mean “even odds after weighing the evidence” but “the
ranking is arbitrary and another sample could reverse it”; and because
it comes from the normal approximation behind `se_diff`, it inherits
that approximation’s failures, which is why *loo* flags them in
`diag_diff` and `diag_elpd` when it detects them. Here it is 1.00, and
the reading is that a sample like this one would essentially never put
the demographic model ahead.

### Totals, Averages, and Why the Convention Is a Total

`elpd_diff` is a difference of *totals*, and the reason is that a log
score is additive. Summing it over observations gives a log predictive
likelihood, so the difference between two models is a log likelihood
ratio: the evidence the sample carries about which model predicts
better, in nats. On that reading the total is the quantity with meaning,
and it should grow with the sample, because more data is more evidence.
It also keeps `elpd` a total over observations, as AIC, WAIC and DIC
are, though those three are quoted on the deviance scale, which is
\\-2\\ times a log score.

Nothing is lost by preferring the average, because the decision does not
depend on the choice. The mean difference and the total differ by a
factor of \\n\\, its standard error differs by the same factor, and the
ratio of the two, which is what says whether the difference is real, is
identical either way.

What the average is better for is a comparison whose two halves do not
sum over the same observations. The section above ran into exactly that:
`elpd_loo` totals 1500 terms and the held-out score totals 300, so the
totals are not comparable and the per-observation averages are. The same
is true of a score quoted across datasets, or across subgroups of
different sizes. Within one
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
call it never arises, since every model there is scored on the same
rows, which is why the convention can afford to be a total.

### The Same Comparison Through the Split

When the Pareto diagnostics say leave-one-out cannot be trusted, this
comparison is the one to rebuild on held-out data. Both candidates are
fitted to the training half and scored on the observations neither of
them saw:

``` r

demographics_train <- bartisan(death ~ rhc + age + sex + race + edu,
                               data = train, family = binomial())

score_demographics <- predict(demographics_train, newdata = held,
                              type = "density", log = TRUE)

d <- score - score_demographics

c(mean_diff = mean(d), se = sd(d) / sqrt(length(d)),
  ratio = mean(d) / (sd(d) / sqrt(length(d))))
#> mean_diff        se     ratio 
#>   0.03691   0.01645   2.24372
```

That is
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)’s
arithmetic done by hand and quoted per observation: the average
difference in log score, its standard error, and the ratio that decides
whether to believe it. The full model is ahead by a little over two
standard errors rather than six, which is what the split costs. It sees
300 observations where
[`loo()`](https://mc-stan.org/loo/reference/loo.html) sees 1500, and the
precision of a difference goes as the square root of the count.

### Every Observation Held Out Once (K-Fold)

A single split is honest but wasteful: it scores 300 observations and
throws the other 1200 into training, where they tell us nothing about
prediction. Splitting the sample into \\K\\ parts instead, and fitting
\\K\\ times so that each part is scored by a model that never saw it,
holds every observation out exactly once.
[`kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html) does
that:

``` r

set.seed(2026)
folds <- loo::kfold_split_random(K = 5, N = nrow(rhc))

kfold_full <- kfold(full, folds = folds)

kfold_full
#> 
#> Based on 5-fold cross-validation.
#> 
#>            Estimate   SE
#> elpd_kfold   -843.7 17.2
#> p_kfold        28.9  2.3
#> kfoldic      1687.5 34.3
```

`elpd_kfold` is the held-out log score, summed over every observation,
and it owes nothing to an importance-sampling approximation: the model
really was refitted without each fold. `p_kfold` is the gap between what
the model predicts for an observation it was fitted to and what it
predicts for the same one held out, which is the price of having used
it. The cost is five fits rather than one.

The folds are passed rather than drawn so that the second model is
scored on the same split, which is what makes the two comparable:

``` r

kfold_demographics <- kfold(demographics, folds = folds)

loo_compare(list(full = kfold_full, demographics = kfold_demographics))
#>         model elpd_diff se_diff p_worse diag_diff diag_elpd
#>          full       0.0     0.0      NA                    
#>  demographics     -72.2    11.1    1.00
```

The same reading as before, and at the same precision as the
leave-one-out comparison rather than the single split’s: every
observation is scored, so nothing is thrown away.

It also gives a way to check
[`loo()`](https://mc-stan.org/loo/reference/loo.html) rather than taking
the Pareto diagnostics’ word for it, since the two estimate the same
thing:

``` r

c(kfold = kfold_full$estimates["elpd_kfold", "Estimate"] / nrow(rhc),
  loo = elpd_total / nrow(rhc))
#>   kfold     loo 
#> -0.5625 -0.5652
```

They agree, which is what a clean Pareto \\k\\ column was saying in less
direct form. That is the division of labor between the three routes:
[`loo()`](https://mc-stan.org/loo/reference/loo.html) for one fit and a
diagnostic, \\K\\-fold for \\K\\ fits and no approximation, and a single
split when even that is too expensive. When the outcome is rare enough
that a random split could leave a fold with no events in it,
[`loo::kfold_split_stratified()`](https://mc-stan.org/loo/reference/kfold-helpers.html)
assigns the folds instead.

This is the right way to ask whether a set of variables earns its place,
and it is a better question than the one variable importance answers,
because it is about prediction rather than about how the forest happened
to spend its splits. Here the answer is not in doubt by any of the three
routes: how sick a patient is on arrival predicts whether they die, and
demographics alone do not.

## Comparing Links (`family`)

The other choice is the likelihood. For a binary outcome the family is
settled, and what remains is the link (i.e., the function that maps the
forest’s output onto a probability). Below we compare probit and
complementary log-log BART models to our original logistic BART model.

``` r

set.seed(2026)
probit <- bartisan(model, data = rhc, family = binomial("probit"))

set.seed(2026)
cloglog <- bartisan(model, data = rhc, family = binomial("cloglog"))

loo_compare(list(logit = loo(full),
                 probit = loo(probit),
                 cloglog = loo(cloglog)))
#>    model elpd_diff se_diff p_worse       diag_diff diag_elpd
#>  cloglog       0.0     0.0      NA                          
#>   probit      -0.3     1.8    0.57 |elpd_diff| < 4          
#>    logit      -1.5     2.1    0.77 |elpd_diff| < 4
```

The three are within a point or two of each other, and the differences
are smaller than their standard errors, which *loo* flags directly. The
reading is that the link does not matter here.

That is a useful negative result and worth reporting as one. It is also
the usual outcome: with a flexible function on the inside, the link has
little left to do, because the forest can absorb the difference between
one link and another. This is not true of a generalized linear model,
where the link carries the whole shape of the relationship.

## The Scale of the Log Density

[`loo()`](https://mc-stan.org/loo/reference/loo.html) is built on the
pointwise log density of each observation, so comparing two models by it
assumes both densities are taken with respect to the same measure. That
is usually automatic. Two different things break it: one is a change of
variable and is always repairable, the other is a point mass and usually
is not.

### A Change of Variable (`ph()` Against an AFT Family)

The accelerated failure time families report the density of \\\log T\\,
while
[`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
reports the density of \\T\\. Below we fit one of each to the survival
outcome in these data and compare them as they come out:

``` r

library(survival)

surv_model <- Surv(days, death) ~ age + sex + race + edu + aps + meanbp + surv2m

set.seed(2026)
aft <- bartisan(surv_model, data = rhc, family = weibull_aft())

set.seed(2026)
prop_haz <- bartisan(surv_model, data = rhc, family = ph())

loo_compare(list(aft = loo(aft), ph = loo(prop_haz)))
#>  model elpd_diff se_diff p_worse diag_diff       diag_elpd
#>    aft       0.0     0.0      NA           2 k_psis > 0.66
#>     ph   -3390.9    90.8    1.00
```

The accelerated failure time model appears to win by thousands of
points, which is implausible on its face: no choice of baseline hazard
is worth that much. What the comparison is measuring is the change of
variable.

Writing \\Y = \log T\\, the density of \\T\\ is \\f_T(t) = f_Y(\log t) /
t\\, so the log density on the \\T\\ scale is \\\log f_Y(\log t) - \log
t\\. Subtracting \\\log t\\ from each of the accelerated failure time
model’s pointwise contributions is the whole of the correction, and it
applies to events only: a censored observation contributes a survival
probability, which is a probability on either scale and carries no
measure to change.

`scale` asks [`loo()`](https://mc-stan.org/loo/reference/loo.html) for
it. Naming one measure for every model in the comparison is enough,
because a fit already on the scale named is returned untouched:

``` r

loo_compare(list(aft = loo(aft, scale = "time"),
                 ph = loo(prop_haz, scale = "time")))
#>  model elpd_diff se_diff p_worse diag_diff       diag_elpd
#>     ph       0.0     0.0      NA                          
#>    aft    -131.0    14.9    1.00           2 k_psis > 0.66
```

The ordering reverses. Read on the scale they share, the two are a
hundred-odd points apart rather than three thousand, and it is the
proportional hazards model that predicts these times better.

The correction is not applied without being asked for, because
[`loo()`](https://mc-stan.org/loo/reference/loo.html) would then stop
reporting the model’s own predictive density: it would no longer agree
with `log_lik()`, and a comparison against a proportional hazards fit
from some other package would quietly take on the error the correction
exists to remove.

Which of the two scales is named does not matter, since a constant per
observation cancels from a difference. Asking for `"log_time"` moves
both totals and leaves the gap where it was:

``` r

c(on_time = elpd(loo(prop_haz, scale = "time")) - elpd(loo(aft, scale = "time")),
  on_log_time = elpd(loo(prop_haz, scale = "log_time")) -
    elpd(loo(aft, scale = "log_time")))
#>     on_time on_log_time 
#>         131         131
```

The everyday version of this mistake takes the same repair by hand. A
model of `log(y)` reports \\\log f_Y(\log y)\\, so subtracting \\\log
y\\ from each of its pointwise contributions, with
`loo(sweep(rstantools::log_lik(fit), 2, log(y), "-"))`, puts it on the
scale of a model of `y`.

There is also a way to sidestep the correction rather than apply it.
Survival at a horizon is a probability under every family, so comparing
fits on `predict(type = "survival")` needs no Jacobian at all;
[`vignette("survival")`](https://ngreifer.github.io/bartisan/articles/survival.md)
gives that route, along with the same correction written for a held-out
log score instead of for
[`loo()`](https://mc-stan.org/loo/reference/loo.html).

### A Point Mass Against a Density (`gaussian()` Against `tweedie()`)

The second break is not a transformation, and unlike the first it
usually cannot be repaired.
[`tweedie()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
places a point mass at zero and a density on the positive half-line;
[`gaussian()`](https://rdrr.io/r/stats/family.html) places a density on
the whole line. Earnings in the `lalonde` data are where this comes up,
since a quarter of the sample earned exactly nothing.

``` r

data("lalonde", package = "cobalt")

earnings <- re78 ~ age + educ + race + married + nodegree + re74 + re75

set.seed(2026)
normal <- bartisan(earnings, data = lalonde, family = gaussian())

set.seed(2026)
compound <- bartisan(earnings, data = lalonde, family = tweedie())

loo_compare(list(gaussian = loo(normal), tweedie = loo(compound)))
#>     model elpd_diff se_diff p_worse diag_diff       diag_elpd
#>   tweedie       0.0     0.0      NA                          
#>  gaussian   -1247.4    95.1    1.00           1 k_psis > 0.66
```

The tweedie is ahead by more than a thousand points, and that number
should not be believed. At an observation where earnings are exactly
zero the tweedie reports a probability and the gaussian reports a
density, which are not the same kind of quantity and do not belong in
the same sum. Splitting the contributions at the zeros shows that the
gap is nothing but that mismatch:

``` r

zero <- lalonde$re78 == 0

rbind(gaussian = c(at_zero = sum(loo(normal)$pointwise[zero, "elpd_loo"]),
                   positive = sum(loo(normal)$pointwise[!zero, "elpd_loo"])),
      tweedie  = c(at_zero = sum(loo(compound)$pointwise[zero, "elpd_loo"]),
                   positive = sum(loo(compound)$pointwise[!zero, "elpd_loo"])))
#>          at_zero positive
#> gaussian -1459.2    -4843
#> tweedie   -213.5    -4841
```

At the positive outcomes, where both report a density, the two are
within a couple of points of each other. Everything else is the atom.

**The one case where the comparison is valid** is an outcome that is
really recorded on a grid rather than being continuous: whole dollars,
whole counts, tenths of a millimeter. Then neither model is reporting a
density at all, properly speaking, and both can be put on the
probability of landing in one cell of that grid, after which the totals
are comparable. Doing it means turning each density into a probability
by multiplying it by the cell width, which is a constant per observation
and so a correction of the same shape as the Jacobian above. The width
has to come from knowing how the outcome was recorded; no fitted model
can supply it, and choosing it to favor a model is choosing the answer.

When the zeros are exact rather than rounded, there is no width to use
and no correction to make. In that case do not compare the two families
by [`loo()`](https://mc-stan.org/loo/reference/loo.html) at all. Two
things answer the question instead, and between them they cover it.

The first is to compare them only where both describe the outcome the
same way, which is the positive observations, as above: on those, here,
neither family predicts better.

The second is to notice that the disagreement is entirely about whether
the outcome has a point mass, and that this is a question about the
shape of the predictive distribution rather than about a density at a
point. A posterior predictive check settles it in one line:

``` r

set.seed(2026)

zero_share <- function(y) mean(y == 0)

c(observed = zero_share(lalonde$re78),
  gaussian = mean(apply(rstantools::posterior_predict(normal), 1, zero_share)),
  tweedie  = mean(apply(rstantools::posterior_predict(compound), 1, zero_share)))
#> observed gaussian  tweedie 
#>   0.2329   0.0000   0.2185
```

A gaussian fit never produces an exact zero and about a quarter of these
outcomes are exactly zero, so
[`tweedie()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is the family to use here and no log score was needed to say so.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers the families that carry a point mass and
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
covers the checks.

The rule behind both halves of this section is that
[`loo()`](https://mc-stan.org/loo/reference/loo.html) compares models
whose densities are taken with respect to the same measure. A
transformation of the outcome is undone with its Jacobian, so those
comparisons are repairable; a point mass against a density is repairable
only when the outcome is recorded on a grid, and otherwise belongs to a
predictive check rather than to a log score.

## Comparing Against a Bayesian GLM

Often it is a good idea to compare a flexible model to a more easily
interpretable (generalized) linear model to assess whether the
flexibility buys us anything. There is no obstacle to this, and it is
worth doing. If a logistic regression predicts as well as the forest,
that is evidence the relationship is close to linear on the log-odds
scale, and the simpler model is the easier one to report.

The comparison has to be like for like, which means both models have to
produce a pointwise log density of the same outcome on the same scale.
Fitting the regression in a Bayesian framework is what arranges that:
*rstanarm* fits it with `stan_glm()` and gives it a
[`loo()`](https://mc-stan.org/loo/reference/loo.html) method, and the
resulting object goes into
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
beside ours.

``` r

logistic <- rstanarm::stan_glm(model, data = rhc, family = binomial(),
                               chains = 4, refresh = 0, seed = 2026)

loo(logistic)
#> 
#> Computed from 4000 by 1500 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -844.9 18.4
#> p_loo        16.6  0.6
#> looic      1689.7 36.8
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume independent draws (r_eff=1).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
```

``` r

loo_compare(list(bart = loo(full), logistic = loo(logistic)))
#>     model elpd_diff se_diff p_worse       diag_diff diag_elpd
#>  logistic       0.0     0.0      NA                          
#>      bart      -2.9     5.1    0.72 |elpd_diff| < 4
```

The forest is behind by roughly half a standard error of the difference,
which is to say the two predict this outcome equally well[^1]. Nothing
was lost by fitting the forest and nothing was gained, and the honest
report of that is the one above: the flexible model was tried and did
not find anything the linear one missed.

`p_loo` says where that came from. The regression has 16 coefficients
and an effective number of parameters to match; the forest’s is about
twice that, which is the flexibility it spent looking for curvature and
interaction that turned out not to be there. A forest that predicts no
better while using twice the parameters is a forest reporting that the
log-odds are close to linear here, which is a finding rather than a
disappointment.

Comparing `elpd_loo` from this package against
[`AIC()`](https://rdrr.io/r/stats/AIC.html) from
[`glm()`](https://rdrr.io/r/stats/glm.html) is not a comparison and
should not be reported as one. The two differ by a factor of \\-2\\
before anything else, and [`AIC()`](https://rdrr.io/r/stats/AIC.html)’s
penalty is a count of parameters, which a forest has no fixed number of;
the `p_loo` above is an estimate rather than a count. Fitting the
regression the Bayesian way is what lets us compare
[`loo()`](https://mc-stan.org/loo/reference/loo.html) against
[`loo()`](https://mc-stan.org/loo/reference/loo.html), as above.

## Tuning

The number of trees, `k`, and the gate should not be chosen by
cross-validation as a matter of routine. The priors are chosen so that
the defaults work across a wide range of problems, and tuning them
typically produces small gains together with an optimistically biased
estimate of performance when the same data chose the setting.

When there is a reason to tune, the shape of it is a grid fixed in
advance, chosen on one set of observations and assessed on another. The
split from the earlier section is already in hand, so the grid goes
through [`loo()`](https://mc-stan.org/loo/reference/loo.html) on the
training set:

``` r

trees <- c(20, 50, 200)

tuned <- lapply(trees, function(n) {
  set.seed(2026)
  bartisan(model, data = train, family = binomial(), num_trees = n)
})

names(tuned) <- paste0("trees_", trees)

loo_compare(lapply(tuned, loo))
#>      model elpd_diff se_diff p_worse       diag_diff diag_elpd
#>  trees_200       0.0     0.0      NA                          
#>   trees_20      -0.4     1.9    0.59 |elpd_diff| < 4          
#>   trees_50      -0.9     1.3    0.76 |elpd_diff| < 4
```

Nothing separates them: every difference is smaller than its own
standard error, and *loo* flags both comparisons as too small to read.
Scoring the same three fits on the held-out observations says it again,
and is the assessment the selection is not allowed to see:

``` r

sapply(tuned, function(f) sum(predict(f, newdata = held, type = "density",
                                      log = TRUE)))
#>  trees_20  trees_50 trees_200 
#>    -171.9    -172.7    -172.2
```

The two orderings disagree, which is the practical content of the
warning rather than a contradiction: when the spread across a grid is
smaller than the noise in the estimate, whichever setting comes out top
is the one the noise favored. Reporting its score as the model’s
performance is exactly the optimism this section opened with, and it is
why the assessment has to come from observations the grid never touched.
Here the honest summary is that `num_trees` does not matter on these
data, and the default is what to keep.

That is the usual outcome, and it is worth knowing before spending a
grid on it. When tuning does pay, it tends to be for cost rather than
accuracy, and
[`vignette("faq")`](https://ngreifer.github.io/bartisan/articles/faq.md)
lists the settings that buy speed;
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
documents what each one changes.

## Variable Selection

Variables should not be selected by fitting many models and keeping the
best. With a search over all subsets, the winner is chosen partly for
fitting the noise, for the same reason stepwise regression is not to be
trusted.

Two things work instead, and which one to use depends on the question.
When the question is whether a *set* of variables earns its place, the
comparison is the one from earlier in this vignette: name the
specifications in advance and put them through
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html),
as with the full and demographic models above. A handful of comparisons
decided beforehand is informative in a way a search is not.

When the question is which *individual* predictor matters, dropping one
at a time and reading
[`loo()`](https://mc-stan.org/loo/reference/loo.html) is a weak
instrument. It works for a dominant predictor and fails for the rest:

``` r

drop_one <- function(v) {
  set.seed(2026)
  reduced <- bartisan(update(model, paste(". ~ . -", v)), data = rhc,
                      family = binomial())

  d <- loo(full)$pointwise[, "elpd_loo"] - loo(reduced)$pointwise[, "elpd_loo"]

  c(elpd_diff = sum(d), se_diff = sqrt(length(d)) * sd(d))
}

rbind(surv2m = drop_one("surv2m"), rhc = drop_one("rhc"))
#>        elpd_diff se_diff
#> surv2m    40.454   9.357
#> rhc        1.787   2.544
```

Losing `surv2m`, the prognostic score, costs about forty points and is
plainly detectable. Losing `rhc` costs about two, which is smaller than
its own standard error, and yet `rhc` is real enough that
[`vignette("causal")`](https://ngreifer.github.io/bartisan/articles/causal.md)
puts its effect at around six percentage points on the probability of
death. The remaining predictors carry that information too, so the
forest reroutes around the one that was removed, and a comparison of
predictive density cannot see a variable whose job something else can
do. Running one such comparison per predictor would also be the search
this section began by warning against.

So that question belongs to
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
and to the sparsity prior, which does the selection inside the model
rather than across refits.
[`vignette("importance")`](https://ngreifer.github.io/bartisan/articles/importance.md)
covers reading it, including the noise-predictor calibration that says
which end of the table carries information, and
[`vignette("effects")`](https://ngreifer.github.io/bartisan/articles/effects.md)
covers the step after: asking what a predictor does to the outcome
rather than whether it is used.

## Further Reading

[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
covers checking that each model fits before comparing them, which is the
step most often skipped.
[`vignette("families")`](https://ngreifer.github.io/bartisan/articles/families.md)
covers what the families assume, which is what a comparison between them
is really about.

## References

Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2017. “Practical Bayesian
Model Evaluation Using Leave-One-Out Cross-Validation and WAIC.”
*Statistics and Computing* 27 (5): 1413–32.
<https://doi.org/10.1007/s11222-016-9696-4>.

[^1]: [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
    warns that the models do not have the same `y` variable, which is
    not what has happened here. It compares a hash of the response that
    *rstanarm* and *brms* attach to their
    [`loo()`](https://mc-stan.org/loo/reference/loo.html) output and
    this package does not, so it is holding a hash up against nothing.
    Attaching one would not help: the hash is taken of the response as
    each package stores it, and an outcome held as `integer` on one side
    and `double` on the other hashes differently, so the warning would
    start firing on identical data rather than stop firing. What we
    should check instead, and directly, is that both models were fitted
    to the same rows.
