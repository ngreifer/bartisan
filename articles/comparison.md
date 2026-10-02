# Model Comparison

## Introduction

Model selection means something different for a BART model than it does
for a linear model, and the difference is worth stating before any code.

With [`glm()`](https://rdrr.io/r/stats/glm.html), we choose which terms
enter, whether to add a squared term, and whether to include an
interaction; those choices are the model. With BART, the forest makes
them, so they are not ours to make. What remains is a shorter list:

1.  Which variables the model is allowed to see.
2.  Which likelihood is used, i.e., which family and link.
3.  Occasionally, a sampler setting (e.g., the number of trees).

In this guide we will compare those choices using leave-one-out
cross-validation. First we’ll fit a model to the `rhc` data and read its
[`loo()`](https://mc-stan.org/loo/reference/loo.html) output, including
the diagnostics that say when the approximation cannot be trusted, and
the two forms of cross-validation to fall back on when it cannot. Next
we’ll compare two sets of variables and then two links, cover the
families whose log densities are on different scales and what to do when
comparing them, and compare a logistic BART model to a Bayesian logistic
regression. Finally we’ll cover how to tune a setting and how to select
variables.

``` r

library(bartisan)
library(loo)

data("rhc")

model <- death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
  hema + pafi + paco2 + crea + surv2m + card

set.seed(2026)

full <- bartisan(model, data = rhc, family = binomial())
```

The fits here use the default single chain for brevity; leave-one-out
needs draws rather than chains. Convergence should still be checked
separately, and the fit controls possibly changed, as described in
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md),
before we trust any of these comparisons.

## Leave-One-Out Cross-Validation (`loo()`)

Leave-one-out cross-validation scores a model by how well it predicts
each observation from a fit that did not see it: the model is fit to
every observation but one, the log predictive density of the held-out
outcome is recorded, and the sum over observations, the expected log
pointwise predictive density (ELPD), measures how well the model
predicts new data. Done literally, that takes one refit per observation.
[`loo()`](https://mc-stan.org/loo/reference/loo.html) approximates it
from the single fit we already have by Pareto-smoothed importance
sampling ([Vehtari et al. 2017](#ref-vehtari2017)). Leaving out
observation \\i\\ changes the posterior only through that observation’s
likelihood, so reweighting each posterior draw by the inverse of its
likelihood for \\y_i\\ turns draws from the full posterior into
approximate draws from the posterior without \\i\\. Pareto smoothing
stabilizes the largest of those weights, and the Pareto \\k\\ diagnostic
below reports how heavy their tail is. All
[`loo()`](https://mc-stan.org/loo/reference/loo.html) needs is the log
likelihood of every observation at every draw, which is why it takes
seconds rather than a refit per observation.

``` r

loo(full)
#> 
#> Computed from 800 by 1500 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -847.6 17.3
#> p_loo        35.7  1.0
#> looic      1695.2 34.6
#> ------
#> MCSE of elpd_loo is 0.5.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.0, 0.5]).
#> 
#> All Pareto k estimates are good (k < 0.66).
#> See help('pareto-k-diagnostic') for details.
```

`elpd_loo` reports the leave-one-out ELPD, with higher values being
better. `looic` is the same measure multiplied by -2 to be on the
deviance scale.

`p_loo` is the effective number of parameters, a little over 35 here;
for a forest with hundreds of leaves across its trees, that number is
small; the prior shrinks most of them toward zero. It can be a useful
measure of how much of the data the model is actually using.

The Pareto \\k\\ diagnostics are worth checking: the leave-one-out
approximation by importance sampling is trustworthy only when the
importance weights are well behaved, and a \\k\\ above the threshold
*loo* prints with them, which is at most .7 and a little lower for a fit
with this many draws, says that for that observation they are not. Here
they are all good.

### Failures of the Approximation

When high Pareto \\k\\ values are observed, the posterior importance
sampling approximation to leave-one-out cross-validation can be
inaccurate. An alternative is to use k-fold cross-validation, which
avoids this approximation but requires refitting the model several
times, which can be computationally expensive. K-fold cross-validation
involves splitting the sample into \\K\\ parts, fitting the model \\K\\
times (leaving out one fold each time), and using that fold to compute
the log likelihood contribution of each unit. This can be done using
[`loo::kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html):

``` r

set.seed(2026)

# K = 5 here; more is better but slower
folds <- kfold_split_random(K = 5, N = nrow(rhc))

# Refit the models K times
kfold_full <- kfold(full, folds = folds)

kfold_full
#> 
#> Based on 5-fold cross-validation.
#> 
#>            Estimate   SE
#> elpd_kfold   -843.5 17.0
#> p_kfold        31.5  2.0
#> kfoldic      1686.9 33.9
```

[`kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html)
produces results that can be interpreted like those from
[`loo()`](https://mc-stan.org/loo/reference/loo.html), but without
relying on the importance sampling approximation. When the outcome is
rare enough that a random split could leave a fold with no events in it,
[`kfold_split_stratified()`](https://mc-stan.org/loo/reference/kfold-helpers.html)
can be used instead of
[`kfold_split_random()`](https://mc-stan.org/loo/reference/kfold-helpers.html)
to assign the folds instead.

## Comparing Two Models (`loo_compare()`)

We can compare the fit of two models by supplying their
[`loo()`](https://mc-stan.org/loo/reference/loo.html) or
[`kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html) output
to
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).
Below, we fit a model that only includes demographic variables as
predictors to compare to our full model.

``` r

set.seed(2026)
demographics <- bartisan(death ~ rhc + age + sex + race + edu,
                         data = rhc, family = binomial())

loo_compare(list(full         = loo(full),
                 demographics = loo(demographics)))
#>         model elpd_diff se_diff p_worse diag_diff diag_elpd
#>          full       0.0     0.0      NA                    
#>  demographics     -67.9    11.2    1.00
```

The full model predicts better by around six times the standard error of
the difference. The standard error is the important half: a difference
of many standard errors is clear, and a difference smaller than its own
standard error is not evidence of anything.

`p_worse` reports that reading as a probability. It is
`pnorm(abs(epld_diff / se_diff))`, the probability that the model on
that row is really the worse of the two given how far apart they came
out and how precisely the difference is known. The best-ranked model has
nothing to be compared against and gets `NA`. Some of its properties are
worth knowing before reading one. Because the models are sorted by
`elpd_loo` before `p_worse` is computed, every reported value is at
least .5 by construction, so .5 does not mean “even odds after weighing
the evidence” but “the ranking is arbitrary and another sample could
reverse it”; and because it comes from the normal approximation behind
`se_diff`, it inherits that approximation’s failures, which is why *loo*
flags them in `diag_diff` and `diag_elpd` when it detects them. Here it
is 1.00, and the reading is that a sample like this one would
essentially never put the demographic model ahead.

When using k-fold cross-validation for model comparison instead of the
importance sampling approximation, it’s important that the same folds
are used across fits so that each model is scored on the same split,
which makes them comparable. We can then supply the
[`kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html) output
to [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
to compare the models.

``` r

kfold_demographics <- kfold(demographics, folds = folds)

loo_compare(list(full         = kfold_full,
                 demographics = kfold_demographics))
#>         model elpd_diff se_diff p_worse diag_diff diag_elpd
#>          full       0.0     0.0      NA                    
#>  demographics     -72.3    10.9    1.00
```

Here, we find the same reading as before, consistent with the Pareto
\\k\\ diagnostics’s indication that the approximation is valid.

This is the right way to ask whether a set of variables earns its place,
and it is a better question than the one variable importance answers,
because it is about prediction rather than about how the forest happened
to spend its splits. Here the answer is not in doubt: how sick a patient
is on arrival predicts whether they die, and demographics alone do not.

## Comparing Links (`family`)

Another model choice is the likelihood. For a binary outcome, the family
is settled, and what remains is the link (i.e., the function that maps
the forest’s output onto a probability). Below we compare a probit BART
model to our original logistic BART model.

``` r

set.seed(2026)
probit <- bartisan(model, data = rhc, family = binomial("probit"))

loo_compare(list(logit  = loo(full),
                 probit = loo(probit)))
#>   model elpd_diff se_diff p_worse       diag_diff diag_elpd
#>   logit       0.0     0.0      NA                          
#>  probit      -0.2     1.0    0.60 |elpd_diff| < 4
```

The two are within a point of each other, and the differences are
smaller than their standard errors, which *loo* flags directly. The
reading is that the link does not matter here.

That is a useful negative result and worth reporting as one. It is also
the usual outcome: with a flexible function on the inside, the link has
little left to do, because the forest can absorb the difference between
one link and another. This is not true of a generalized linear model,
where the link carries the whole shape of the relationship.

## Comparing Against a Bayesian GLM

Often it is a good idea to compare a flexible model to a more easily
interpretable (generalized) linear model to assess whether the
flexibility buys us anything. There is no obstacle to this, and it is
worth doing. If a logistic regression predicts as well as the forest,
that is evidence the relationship is close to linear on the log-odds
scale, and the simpler model is the easier one to report.

The comparison has to be like for like, which means both models have to
produce a pointwise log density of the same outcome on the same scale.
Fitting the regression in a Bayesian framework arranges that: *rstanarm*
fits it with `stan_glm()` and gives it a
[`loo()`](https://mc-stan.org/loo/reference/loo.html) method, and the
resulting object goes into
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
beside ours.

``` r

# Fit a Bayesian logistic GLM
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

loo_compare(list(bart     = loo(full),
                 logistic = loo(logistic)))
#>     model elpd_diff se_diff p_worse       diag_diff diag_elpd
#>  logistic       0.0     0.0      NA                          
#>      bart      -2.8     5.0    0.71 |elpd_diff| < 4
```

The forest is behind by roughly half a standard error of the difference,
which is to say the two predict this outcome equally well[^1]. Nothing
was lost by fitting the forest and nothing was gained, and the honest
report of that is the one above: the flexible model was tried and did
not find anything the linear one missed.

## Tuning

The number of trees, `k`, and the gate should not be chosen by
cross-validation as a matter of routine. The priors are chosen so that
the defaults work across a wide range of problems, and tuning them
typically produces small gains together with an optimistically biased
estimate of performance when the same data chose the setting. When there
is a reason to tune, the shape of it is a grid fixed in advance.
Leave-one-out or K-fold cross-validation can be used to compare the
fits:

``` r

trees <- c(20, 50, 200)

tuned <- lapply(trees, function(n) {
  set.seed(2026)
  bartisan(model, data = rhc, family = binomial(),
           num_trees = n)
})

names(tuned) <- paste0("trees_", trees)

# LOOCV
loo_compare(lapply(tuned, loo))
#>      model elpd_diff se_diff p_worse       diag_diff diag_elpd
#>  trees_200       0.0     0.0      NA                          
#>   trees_50       0.0     1.2    0.51 |elpd_diff| < 4          
#>   trees_20      -2.2     2.0    0.87 |elpd_diff| < 4

# K-fold CV
# loo_compare(lapply(tuned, kfold, folds = folds))
```

In this case, nothing separates them: the differences are a couple of
points at most, about the size of their standard errors. Here the honest
summary is that `num_trees` does not matter on these data, and the
default can be kept, which is the usual outcome. When tuning does pay,
it tends to be for computational cost rather than accuracy, and
[`vignette("faq")`](https://ngreifer.github.io/bartisan/articles/faq.md)
lists the settings that buy speed;
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
documents what each one changes.

## Additional Topics

### A Constant Coefficient Against a Varying One (`vc()`)

`vc(z)` gives `z` a coefficient that is itself a forest, free to vary
with the other predictors, and `vc(z, ~ 1)` pins that coefficient to a
single number, so `z` enters as a linear term while everything else
stays nonparametric. Comparing the two by
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
asks a question neither variable importance nor a dropped predictor
answers: not whether `z` matters, but whether what it does depends on
anything else.
[`vignette("varying")`](https://ngreifer.github.io/bartisan/articles/varying.md)
works that comparison through on this data, and covers what to keep in
mind when reading it: that the constant coefficient is drawn under the
leaf prior and so is shrunk toward zero, and that the comparison is far
better powered on a Gaussian outcome than on a binary one.

### The Scale of the Log Density (`scale`)

[`loo()`](https://mc-stan.org/loo/reference/loo.html) compares models by
the log density each assigns to the observed outcomes, so the models
being compared have to describe the same quantity on the same scale.
That holds for every comparison above, and for any comparison between
two families of the same kind. It needs a step from us when a comparison
crosses between the two kinds of survival family.

The accelerated failure time families
([`weibull_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
[`loglogistic_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
[`lognormal_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
and
[`dpm_aft()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md))
model the log of the survival time, and
[`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
models the time itself. When a comparison includes both kinds,
`scale = "time"` should be supplied to every call to
[`loo()`](https://mc-stan.org/loo/reference/loo.html), which puts each
model on the scale of the time. A fit already on that scale is returned
unchanged, so the same argument can be given to every model in the
comparison. Below, we compare a log-normal accelerated failure time
model with a proportional hazards model for the survival outcome in
`rhc`:

``` r

library(survival)

surv_model <- Surv(days, death) ~ age + sex + race + edu +
                 aps + meanbp + surv2m

set.seed(2026)
aft <- bartisan(surv_model, data = rhc, family = lognormal_aft())

set.seed(2026)
prop_haz <- bartisan(surv_model, data = rhc, family = ph())

loo_compare(list(aft = loo(aft,      scale = "time"),
                 ph  = loo(prop_haz, scale = "time")))
#>  model elpd_diff se_diff p_worse diag_diff       diag_elpd
#>     ph       0.0     0.0      NA           1 k_psis > 0.66
#>    aft     -90.0    15.4    1.00
```

The proportional hazards model predicts these survival times better, by
about 90 points with a standard error of 15. The flag in the last column
says that one of the pointwise estimates for
[`ph()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
is unreliable (see the section on failures of the approximation above);
however, with a difference of nearly six standard errors, it does not
change the reading.

### Variable Selection

Variables should not be selected by fitting many models and keeping the
best. With a search over all subsets, the winner is chosen partly for
fitting the noise, for the same reason stepwise regression is not to be
trusted.

Other approaches work instead, and which one to use depends on the
question. When the question is whether a *set* of variables earns its
place, the comparison is the one from earlier in this vignette: name the
specifications in advance and put them through
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html),
as with the full and demographic models above. A handful of comparisons
decided beforehand is informative in a way a search is not.

When the question is which *individual* predictor matters, dropping one
at a time and reading
[`loo()`](https://mc-stan.org/loo/reference/loo.html) is a weak
instrument. That question belongs to
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
on a fit with `sparsity = TRUE`, whose prior does the selection inside
the model rather than across refits.
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
covers what the families assume, which is the real subject of a
comparison between them.

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
    What one should check instead, and directly, is that both models
    were fitted to the same rows.
