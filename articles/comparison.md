# Choosing between models

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
the diagnostics that say when the approximation cannot be trusted. Next
we’ll compare two sets of variables and then three links, and we’ll work
through the one case in this package where two log densities are not on
the same scale. Finally we’ll cover the choices that should not be made
by cross-validation at all.

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

`p_loo` is the effective number of parameters, about 32 here; for a
forest with hundreds of leaves across its trees that number is small,
because the prior shrinks most of them toward zero. It is a useful
measure of how much of the data the model is actually using.

The Pareto \\k\\ diagnostics are the thing to check: leave-one-out by
importance sampling is trustworthy only when the importance weights are
well behaved, and a \\k\\ above .7 says that for that observation they
are not. Here they are all good.

### When the Approximation Fails

A forest is a flexible function, so a single observation can have a good
deal of influence on the leaves it falls into, and high \\k\\ values are
therefore more common than they are for a parametric model. If many
observations are flagged, the estimate is unreliable, and the remedy is
held-out data rather than a different diagnostic:

``` r

train_id <- sample.int(nrow(rhc), 1200)

train <- rhc[train_id, ]
held  <- rhc[-train_id, ]

fit <- bartisan(model, data = train, family = binomial())

sum(predict(fit, newdata = held, type = "density", log = TRUE))
```

This is the same quantity `elpd_loo` approximates, computed directly.

## Comparing Two Models (`loo_compare()`)

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

This is the right way to ask whether a set of variables earns its place,
and it is a better question than the one variable importance answers,
because it is about prediction rather than about how the forest happened
to spend its splits. Here the answer is not in doubt: how sick a patient
is on arrival predicts whether they die, and demographics alone do not.

## Comparing Links (`family`)

The other choice is the likelihood. For a binary outcome the family is
settled, and what remains is the link (i.e., the function that maps the
forest’s output onto a probability).

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
assumes that both densities are taken with respect to the same measure.
That is normally automatic, but there is one case in this package where
it is not.

The accelerated failure time families report the density of \\\log T\\,
while
[`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
reports the density of \\T\\. The two differ by a Jacobian, so a
comparison across that boundary is off by \\\sum \log t\\, which can
amount to hundreds of points and can reverse the ordering.
[`vignette("survival")`](https://ngreifer.github.io/bartisan/articles/survival.md)
gives the correction.

The general rule is that
[`loo()`](https://mc-stan.org/loo/reference/loo.html) compares models
fitted to the same outcome on the same scale. If we transform the
outcome, refit, and compare, the comparison is invalid unless the
Jacobian is accounted for; comparing a model of `y` to a model of
`log(y)` is the everyday version of this mistake.

## Comparing Against a Logistic Regression

There is no obstacle to this, and it is worth doing. If a logistic
regression predicts as well as the forest, that is evidence the
relationship is close to linear on the log-odds scale, and the simpler
model is the easier one to report.

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
which is to say the two predict this outcome equally well. Nothing was
lost by fitting the forest and nothing was gained, and the honest report
of that is the one above: the flexible model was tried and did not find
anything the linear one missed.

`p_loo` says where that came from. The regression has 16 coefficients
and an effective number of parameters to match; the forest’s is about
twice that, which is the flexibility it spent looking for curvature and
interaction that turned out not to be there. A forest that predicts no
better while using twice the parameters is a forest reporting that the
log-odds are close to linear here, which is a finding rather than a
disappointment.

Two cautions on the mechanics.
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
warns that the models do not have the same `y` variable, which is not
what has happened here. It compares a hash of the response that
*rstanarm* and *brms* attach to their
[`loo()`](https://mc-stan.org/loo/reference/loo.html) output and this
package does not, so it is holding a hash up against nothing. Attaching
one would not help: the hash is taken of the response as each package
stores it, and an outcome held as `integer` on one side and `double` on
the other hashes differently, so the warning would start firing on
identical data rather than stop firing. What we should check instead,
and directly, is that both models were fitted to the same rows.

And comparing `elpd_loo` from this package against
[`AIC()`](https://rdrr.io/r/stats/AIC.html) from
[`glm()`](https://rdrr.io/r/stats/glm.html) is not a comparison and
should not be reported as one. The two differ by a factor of \\-2\\
before anything else, and [`AIC()`](https://rdrr.io/r/stats/AIC.html)’s
penalty is a count of parameters, which a forest has no fixed number of;
the `p_loo` above is an estimate rather than a count. Fitting the
regression the Bayesian way is what lets us compare
[`loo()`](https://mc-stan.org/loo/reference/loo.html) against
[`loo()`](https://mc-stan.org/loo/reference/loo.html), as above.

## Tuning and Variable Selection

The number of trees, `k`, and the gate should not be chosen by
cross-validation as a matter of routine. The priors are chosen so that
the defaults work across a wide range of problems, and tuning them
typically produces small gains together with an optimistically biased
estimate of performance when the same data chose the setting. When we do
tune, data should be held out for the final assessment.

Nor should variables be selected by fitting many models and keeping the
best. With a handful of candidate specifications chosen in advance,
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html) is
informative; with a search over all subsets it is not, for the same
reason stepwise regression is not.

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
