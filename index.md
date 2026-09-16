# *bartisan*: Generalized Bayesian Additive Regression Trees

## Overview

*bartisan* fits Bayesian additive regression trees (BART) for responses
that standard BART cannot reach, using the Laplace-approximation
reversible-jump sampler of Linero (2025), which lifts the requirement
that the leaf parameters be integrable in closed form. The primary
function is
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
and a model is specified the way it is in
[`glm()`](https://rdrr.io/r/stats/glm.html): a formula, a data frame,
and a family, with a forest in place of the linear predictor, so nothing
has to be said about which terms are curved or which interact. The
[`stats::family`](https://rdrr.io/r/stats/family.html) objects
[`glm()`](https://rdrr.io/r/stats/glm.html) takes are accepted
unchanged, links included. The families that have no
[`glm()`](https://rdrr.io/r/stats/glm.html) counterpart are supplied in
the same style: negative binomial, ordinal, multinomial, beta and
ordered beta, zero-inflated counts, a Tweedie compound Poisson,
accelerated failure time and proportional hazards models for
right-censored times, location-scale regression, a Dirichlet process
mixture for the error distribution, and a log density written as an R
function. Decision rules are soft by default, in the manner of Linero
and Yang (2018), which gives smoother fits than the step functions of
standard BART, and `gate = "hard"` gives those back. Treatment effects
come from
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md),
and a fit is otherwise read through the packages that already do that
work:
[*marginaleffects*](https://CRAN.R-project.org/package=marginaleffects)
for the wider range of counterfactual estimands,
[*loo*](https://CRAN.R-project.org/package=loo) for model comparison,
[*bayesplot*](https://CRAN.R-project.org/package=bayesplot) and
[*posterior*](https://CRAN.R-project.org/package=posterior) for
convergence diagnostics, and
[*performance*](https://CRAN.R-project.org/package=performance) for fit
statistics.

Check out the *bartisan*
[website](https://ngreifer.github.io/bartisan/)!

## Installation

*bartisan* is not yet on CRAN. You can install the development version
from [GitHub](https://github.com/ngreifer/bartisan) with:

``` r

# install.packages("pak")
pak::pak("ngreifer/bartisan")
```

Installation compiles C++, so it requires a C++17 toolchain along with
*Rcpp* and *RcppArmadillo*.

## Example

### Right heart catheterization in critically ill patients

`rhc` holds 1500 patients from the SUPPORT study, with an indicator for
whether each was given right heart catheterization within 24 hours of
admission to an intensive care unit, whether they died during follow-up,
and how long that took (Connors et al., 1996). Thirteen covariates
recorded before the catheter went in come with them, and
[`?rhc`](https://ngreifer.github.io/bartisan/reference/rhc.md) lists
them. Catheterization was not randomized, so sicker patients were
likelier to receive it, and those covariates are what adjustment has to
work with. Nothing in the model supplies the assumption that they
suffice; what it supplies is a flexible estimate of the outcome under
each treatment given the covariates, which is what g-computation then
averages.

``` r

library(bartisan)

data("rhc")
set.seed(123)

# Death as the outcome, with every variable but the timing of the same event
# as a candidate predictor. The family is the one `glm()` would be given
fit <- bartisan(death ~ . - days, data = rhc,
                family = binomial("probit"))

fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ . - days, data = rhc, family = binomial("probit"))
#> 
#> Family: "binomial" with the "probit" link
#> Observations: 1500
#> Structure: 1 forest of 50 trees, soft decision rules
#> Draws: 800 kept after 200 warmup
```

The forest finds its own interactions and nonlinearity, so what it
reports in place of a coefficient table is how it spent its splitting
rules:

``` r

head(variable_importance(fit), 5)
#> Variable importance
#> 
#>  variable prop_used prop_splits splits
#>    surv2m     1.000       0.318   23.9
#>       age     1.000       0.147   11.0
#>      pafi     1.000       0.099    7.5
#>     paco2     0.954       0.108    8.1
#>       rhc     0.945       0.062    4.7
#> 
#> ℹ splits_lower and splits_upper hold the 95% interval, not shown above.
```

Having no coefficients also means the effect of catheterization is a
contrast between what the model predicts with `rhc` set to one and with
it set to zero, averaged over the covariate distribution in the sample.
That is g-computation, and every posterior draw goes through it, so the
interval is a credible interval rather than a delta-method
approximation:

``` r

estimate_effect(fit, treat = "rhc")
#> Average treatment effect (difference)
#> 
#> Treatment: "rhc"
#> Averaged over 1500 units
#> 
#>     contrast estimate lower upper    n
#>  Y[1] - Y[0]   0.0562     0 0.105 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.634 0.603 0.665
#>      Y[1]    0.690 0.652 0.726
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with "rhc" set to a.
```

The same event is also a right-censored survival time, and asking for it
that way changes the left-hand side of the formula and the family and
nothing else. Here the contrast is a ratio, so the estimate is the
factor by which catheterization multiplies expected survival:

``` r

sfit <- bartisan(survival::Surv(days, death) ~ ., data = rhc,
                 family = lognormal_aft())

estimate_effect(sfit, treat = "rhc", comparison = "ratio")
#> Average treatment effect (ratio)
#> 
#> Treatment: "rhc"
#> Averaged over 1500 units
#> 
#>     contrast estimate lower upper    n
#>  Y[1] / Y[0]    0.595  0.43 0.806 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]      248   201   308
#>      Y[1]      146   113   189
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with "rhc" set to a.
```

[`vignette("bartisan")`](https://ngreifer.github.io/bartisan/articles/bartisan.md)
works through a complete analysis, from fitting to convergence checks to
the estimates, and the other vignettes take those pieces one at a time.
[`vignette("implementation")`](https://ngreifer.github.io/bartisan/articles/implementation.md)
sets out how the sampler works and compares the package’s features
against *dbarts*, *BART*, *flexBART*, *SoftBart*, *bartMachine*, and
*stochtree*.

## Citing *bartisan*

*bartisan* implements the sampler of Linero (2025), and its MCMC engine
is adapted from that paper’s `FlexBart` reference implementation. Please
cite both the method and the package, the latter with the version
number, which `citation("bartisan")` supplies:

Linero, A. R. (2025). Generalized Bayesian additive regression trees
models: beyond conditional conjugacy. *Journal of the American
Statistical Association*, 120(549), 356–369.
<https://doi.org/10.1080/01621459.2024.2337156>

Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles
that adapt to smoothness and sparsity. *Journal of the Royal Statistical
Society Series B*, 80(5), 1087–1110.
<https://doi.org/10.1111/rssb.12293>

## Questions and Bug Reports

Please file an issue at <https://github.com/ngreifer/bartisan/issues>,
with a reproducible example where the report concerns a fit.
