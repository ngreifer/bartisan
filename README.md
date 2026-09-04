
<!-- README.md is generated from README.Rmd. Please edit that file -->

# *bartisan*: Generalized Bayesian Additive Regression Trees

<!-- badges: start -->

<!-- badges: end -->

## Overview

*bartisan* fits Bayesian additive regression trees (BART) for
likelihoods outside the conditionally conjugate Gaussian case, which
standard BART implementations cannot reach. The primary function is
`bartisan()`, and its interface follows that of `glm()`: a model is
specified with a formula, a data frame, and a family, with the forest
taking the place of the linear predictor, so nothing has to be said
about which terms enter the model, which are curved, or which interact.
Ordinary `family` objects are used unchanged, including their links, so
moving a model from `glm()` to `bartisan()` is a one-word change.
Supported families include binomial, Poisson, negative binomial, gamma,
ordinal, multinomial, zero-inflated counts, beta and ordered beta, three
accelerated failure time models and a discrete proportional hazards
model for right-censored data, location-scale regression for a variance
that varies with the predictors, a Dirichlet process mixture for the
error distribution itself, and a likelihood written as an R function.
Decision rules may be hard, as in standard BART, or soft, as in the
SoftBart model of Linero and Yang (2018), which yields smoother fits. A
fitted model is read through the packages that already do that work:
[*marginaleffects*](https://CRAN.R-project.org/package=marginaleffects)
for counterfactual estimands with posterior intervals,
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

Installation compiles C++ code, so it requires a C++17 toolchain along
with *Rcpp* and *RcppArmadillo*.

## Examples

### Right heart catheterization in critically ill patients

The `rhc` dataset holds 1500 patients from the SUPPORT study, with an
indicator for whether each received right heart catheterization within
24 hours of admission to an intensive care unit, whether they died
during follow-up, how long that took, and thirteen physiological
covariates recorded before catheterization (Connors et al., 1996).
Catheterization was not randomized, so the comparison is confounded by
how sick each patient was on admission, and the covariates are what
adjustment has to work with. Two features make the data a useful
illustration: the covariates enter through a forest, so no functional
form has to be committed to, and the same event is available both as a
binary indicator and as a right-censored survival time, so it supports
two families without changing the call’s shape.

``` r
library(bartisan)

data("rhc")
set.seed(123)

# `days` is the timing of the same event as `death`, so it is excluded
# rather than conditioned on. The family is read off the response
fit <- bartisan(death ~ . - days, data = rhc,
                num_trees = 20, num_burn = 200, num_draws = 200)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ . - days, data = rhc, num_trees = 20, 
#>     num_burn = 200, num_draws = 200)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 1 forest of 20 trees, soft decision rules
#> Draws: 200 kept after 200 warmup
```

`summary()` reports how often each predictor was split on, which is the
forest’s account of what it used:

``` r
summary(fit)
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ . - days, data = rhc, num_trees = 20, 
#>     num_burn = 200, num_draws = 200)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 1 forest of 20 trees, soft decision rules
#> Draws: 200
#> 
#> Predictor usage
#> Splitting rules per draw, and how often used at all.
#>         mean    sd lower upper prop_used
#> age    4.230 1.689     2 8.000     1.000
#> pafi   2.690 1.461     1 6.000     1.000
#> surv2m 5.010 1.881     2 9.000     1.000
#> rhc    2.315 1.154     1 5.000     0.995
#> paco2  2.500 1.662     0 6.000     0.965
#> aps    3.165 1.972     0 7.025     0.935
#> card   1.695 1.401     0 5.000     0.795
#> meanbp 1.490 1.211     0 4.025     0.790
#> edu    1.315 1.110     0 4.000     0.750
#> hema   1.870 1.636     0 5.025     0.745
#> resp   1.655 1.406     0 5.000     0.740
#> crea   1.635 1.751     0 6.000     0.675
#> race   1.285 1.447     0 5.000     0.610
#> sex    0.670 0.967     0 3.000     0.425
```

A forest has no coefficients, so the effect of catheterization is a
contrast between what the model predicts under one value of the
treatment and under another, averaged over the covariate distribution in
the sample. That is g-computation, and *marginaleffects* performs it on
a fit directly, reporting a posterior interval with the estimate:

``` r
marginaleffects::avg_comparisons(fit, variables = "rhc")
#> 
#>  Estimate  2.5 % 97.5 %
#>    0.0583 0.0148  0.106
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
```

The same event as a right-censored survival outcome, given as a `Surv`
object on the left-hand side of the formula:

``` r
bartisan(survival::Surv(days, death) ~ . , data = rhc,
         family = lognormal_aft())
```

Prediction returns the full posterior rather than a point estimate,
since every draw of every tree is retained:

``` r
draws <- predict(fit, newdata = rhc[1:5, ], type = "response", draws = TRUE)

apply(draws, 2, quantile, c(.025, .975))
#>            [,1]      [,2]      [,3]      [,4]      [,5]
#> 2.5%  0.7173040 0.7054182 0.1091663 0.2504515 0.2803132
#> 97.5% 0.8687608 0.8884151 0.3109229 0.5071383 0.5910605
```

`vignette("bartisan")` walks through a complete analysis, and the
vignettes it links to take the pieces on their own:
`vignette("families")` for choosing a likelihood, `vignette("survival")`
for censored responses, `vignette("effects")` for reading the fitted
function, `vignette("diagnostics")` for convergence and fit,
`vignette("importance")` for which predictors matter,
`vignette("comparison")` for choosing between models, and
`vignette("causal")` for causal inference.

## Comparison With Other BART Packages

Checked against the installed versions of each package rather than from
memory: *dbarts* 0.9.34, *BART* 2.9.10, *flexBART* 2.0.3, *SoftBart*
1.0.3, *bartMachine* 1.4.2, and *stochtree* 0.4.5. A dash means the
package does not offer the feature, not that it fits it badly.

|  | bartisan | dbarts | BART | flexBART | SoftBart | bartMachine | stochtree |
|----|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **Responses** |  |  |  |  |  |  |  |
| Gaussian | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Binary, probit | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Binary, logit | ✓ | — | ✓ | ✓ | — | — | — |
| Count: Poisson, negative binomial | ✓ | — | — | — | — | — | — |
| Gamma, beta, ordered beta | ✓ | — | — | — | — | — | — |
| Ordinal | ✓ 3 links | — | — | — | — | — | ✓ cloglog |
| Multinomial | ✓ logit, probit | — | ✓ | — | — | — | — |
| Zero-inflated counts | ✓ | — | — | — | — | — | — |
| Survival | ✓ 3 AFT, PH | — | ✓ **+ recurrent, competing risks** | — | — | — | — |
| Heteroskedastic | ✓ | — | — | ✓ | — | — | ✓ |
| Error distribution itself modeled | ✓ | — | — | — | — | — | — |
| A likelihood written by the user | ✓ | — | — | — | — | — | — |
| **Rules and priors** |  |  |  |  |  |  |  |
| Hard decision rules | ✓ | ✓ | ✓ | ✓ | — | ✓ | ✓ |
| Soft decision rules | ✓ 4 gates | — | — | — | ✓ | — | — |
| Dirichlet sparsity prior (DART) | ✓ | — | ✓ | ✓ | ✓ | — | — |
| Splitting weights fixed by the user | ✓ | — | — | — | ✓ groups | ✓ | — |
| Categorical splits on level subsets | ✓ | — | — | ✓ **+ nested, network** | — | — | — |
| Missing predictors, no imputation | ✓ default | — | — | — | — | ✓ | — |
| **Structure** |  |  |  |  |  |  |  |
| Random intercepts | ✓ | ✓ | — | — | — | — | ✓ |
| Varying coefficients | ✓ | — | — | ✓ | ✓ | — | ✓ |
| Treatment-effect (BCF) structure | ✓ | — | — | — | — | — | ✓ |
| Formula interface | ✓ | ✓ | — | ✓ | — | — | — |
| **Running it** |  |  |  |  |  |  |  |
| Several chains | ✓ | ✓ | ✓ | ✓ | — | ✓ | ✓ |
| Threads inside one chain | — | ✓ | ✓ | — | — | ✓ | ✓ |
| Grow-from-root warm start | — | — | — | — | — | — | ✓ |
| Cross-validation over hyperparameters | — | ✓ | — | — | — | ✓ | — |
| Save and reload a fitted model | RDS | — | — | — | — | ✓ | ✓ JSON |
| **Reading the fit** |  |  |  |  |  |  |  |
| `predict()` on new data | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Convergence diagnostics built in | ✓ | — | ✓ | — | — | ✓ | — |
| Variable importance | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| Formal variable-selection test | ✎ | — | ✓ | — | ✓ | ✓ | — |
| Partial dependence | ✎ | ✓ | — | — | ✓ | ✓ | — |
| Interaction detection | ✎ | — | — | — | — | ✓ | — |
| Counterfactual estimands with intervals | ✓ | — | — | — | — | — | — |
| Cross-validated model comparison | ✓ | — | — | — | — | — | — |

**✎ means a helper package covers it, not this one.** All three of those
rows are things *bartMachine* does natively and *bartisan* does through
*marginaleffects*, because a fit works with it and every estimand there
is computed by pushing the draws through: `plot_predictions()` is a
partial dependence plot, and
`avg_comparisons(variables = "x", by = "z")` and `hypotheses(~pairwise)`
are interaction detection, neither of which needs code here. What is
genuinely missing is a *formal* variable-selection test: this package
reports split counts and the share of draws that used a predictor, which
is not a test, where *bartMachine* permutes the response and *SoftBart*
reports posterior inclusion probabilities.

The other columns worth reading as gaps rather than as differences are
these. There are **no threads inside a chain**, so a single chain is
single-core here where *dbarts*, *BART*, *bartMachine*, and *stochtree*
all use several; chains do run in parallel, which is the cheaper win,
but it does not help one chain. There is **no cross-validation** over
`k`, `num_trees`, and the tree prior, which `dbarts::xbart()` and
`bartMachine::bartMachineCV()` both automate. There is **no
grow-from-root warm start**, which is *stochtree*’s way of shortening
burn-in. And there is **no JSON serialization**, so a fit round-trips
through `saveRDS()` and not into another language. On survival, *BART*
is ahead: recurrent events and competing risks are there and here they
are not.

Three more packages are single-purpose rather than general, so they are
not columns above: [*bcf*](https://CRAN.R-project.org/package=bcf) fits
Bayesian causal forests only,
[*VCBART*](https://github.com/skdeshpande91/VCBART) varying-coefficient
models only, and *bartCause* wraps *dbarts* for causal estimands.

## Citing *bartisan*

*bartisan* implements the sampler of Linero (2025), and its MCMC engine
is adapted from that paper’s `FlexBart` reference implementation. Please
cite both the method and the package, the latter with its version
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
