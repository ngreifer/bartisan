# Interfaces to other packages

Methods that let a `<bartisan_fit>` object be used by the packages that
assess model fit. There is nothing to set up: load the other package and
call its function on the fit.

## Usage

``` r
# S3 method for class 'bartisan_fit'
posterior_predict(
  object,
  newdata = NULL,
  iterations = NULL,
  offset = NULL,
  weights = NULL,
  ...
)

# S3 method for class 'bartisan_fit'
posterior_epred(object, newdata = NULL, ...)

# S3 method for class 'bartisan_fit'
posterior_linpred(object, transform = FALSE, newdata = NULL, ...)

# S3 method for class 'bartisan_fit'
log_lik(object, newdata = NULL, ...)

# S3 method for class 'bartisan_fit'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'bartisan_fit'
fitted(object, type = "response", ...)

# S3 method for class 'bartisan_fit'
residuals(object, ...)

# S3 method for class 'bartisan_fit'
weights(object, ...)

# S3 method for class 'bartisan_fit'
sigma(object, ...)

# S3 method for class 'bartisan_fit'
loo(x, ...)

# S3 method for class 'bartisan_fit'
waic(x, ...)

# S3 method for class 'bartisan_fit'
pp_check(object, type = "dens_overlay", ndraws = 10, ...)

# S3 method for class 'bartisan_fit'
as_draws(x, eta = TRUE, ...)

# S3 method for class 'bartisan_fit'
r2_posterior(model, verbose = TRUE, ...)

# S3 method for class 'bartisan_fit'
r2(model, ...)

# S3 method for class 'bartisan_fit'
model_performance(model, metrics = "all", verbose = TRUE, ...)
```

## Arguments

- object, model, x:

  a `<bartisan_fit>` object; the output of a call to
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- newdata:

  optional; a data frame at which to predict. Default is `NULL` to use
  the data the model was fit to. Missing predictor values are allowed in
  the columns that had them when the model was fit, since only those
  columns' splitting rules carry an answer for one; a missing value
  anywhere else is an error. See section *Missing Predictor Values* in
  the Details of
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- iterations:

  `numeric`; optional indices of the stored draws to use, between 1 and
  the number of draws the fit retains. Default is `NULL` to use all of
  them.

- offset, weights:

  an offset and prior weights for `newdata`, as in
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md).
  Defaults are `NULL` to use those the model was fit with. For a
  binomial response the weights are the numbers of trials, and so are
  what a replicate outcome is a fraction of; they must be given
  alongside `newdata` when the model was fit with more than one trial,
  since the number of trials is not a function of the predictors and
  cannot be reconstructed.

- ...:

  further arguments, passed to whatever the method calls (the bayesplot
  check from `pp_check()`,
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) and
  [`loo::waic()`](https://mc-stan.org/loo/reference/waic.html) from
  `loo()` and `waic()`, and
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
  from the rest) and ignored where there is nowhere to pass them.

- transform:

  `logical`; for `posterior_linpred()`, whether to map the predictor
  through the inverse link, which is what `posterior_epred()` does.
  Default is `FALSE`.

- nsim, ndraws:

  `numeric`; the number of posterior draws to use, chosen at random from
  the retained ones. Defaults are 1 for
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) and 10 for
  `pp_check()`.

- seed:

  optional seed, set with
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before drawing and
  restored afterwards, following the
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html)
  convention. Default is `NULL` to leave the stream alone.

- type:

  string; for [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  the prediction scale, passed to
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md);
  default is `"response"`. For `pp_check()`, the name of the bayesplot
  check to run without its `ppc_` prefix, so that `"dens_overlay"` (the
  default) calls
  [`bayesplot::ppc_dens_overlay()`](https://mc-stan.org/bayesplot/reference/PPC-distributions.html)
  ;
  [`bayesplot::available_ppc()`](https://mc-stan.org/bayesplot/reference/available_ppc.html)
  lists them.

- eta:

  for `as_draws()`, which columns of the additive predictor to carry
  into the draws array alongside the scalar parameters, given as either
  a logical value or a numeric vector. Default is `TRUE`, which takes a
  representative ten spread across the range of the fitted function;
  `FALSE` takes none, and a numeric vector takes those observations. The
  default takes a handful rather than all of them because there is one
  column per observation, and an array with thousands of them is not
  something
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
  or a trace plot can be pointed at. The predictor is the quantity whose
  convergence usually matters, and the one
  [`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
  reports on, so it is included by default.

- verbose:

  `logical`; whether to report problems that do not stop the
  computation, such as a family that has no mean and so no Bayesian
  \\R^2\\. Default is `TRUE`.

- metrics:

  `character`; for `model_performance()`, which fit statistics to
  report. Allowable options include `"all"` (the default), `"ELPD"`,
  `"LOOIC"`, `"WAIC"`, `"R2"`, `"RMSE"`, and `"SIGMA"`, and a vector of
  them selects several.

## Value

`posterior_predict()`, `posterior_epred()`, `posterior_linpred()` and
`log_lik()` return a matrix of draws by observations.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) returns a data
frame of one column per replicate. `loo()` and `waic()` return the
`<loo>` and `<waic>` objects those functions produce, and
`model_performance()` a one-row data frame of class
`<performance_model>`. `as_draws()` returns a `<draws_array>` of
iterations by chains by parameters. The accessors return what their
names suggest.

## Details

### What Is Available

**Posterior predictions.**
[`rstantools::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
draws replicate outcomes from the fitted model,
[`rstantools::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
gives their mean and
[`rstantools::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html)
the additive predictor, following the rstantools conventions that brms
and rstanarm follow, and
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) is the same
thing in the shape base R expects. **Pointwise likelihood.**
[`rstantools::log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
returns the draws-by-observations matrix of log-likelihood
contributions, which is what
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) and
[`loo::waic()`](https://mc-stan.org/loo/reference/waic.html) need; both
have methods here.

**Graphical checks.** `pp_check()` runs any of the bayesplot
posterior-predictive checks on the fit. **Summaries.**
[`performance::model_performance()`](https://easystats.github.io/performance/reference/model_performance.html)
collects the fit statistics in one table,
[`performance::r2()`](https://easystats.github.io/performance/reference/r2.html)
gives the Bayesian \\R^2\\, and
[`posterior::as_draws()`](https://mc-stan.org/posterior/reference/draws.html)
hands the scalar parameters to
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
or to the bayesplot MCMC diagnostics. **Basic accessors.**
[`stats::fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`stats::residuals()`](https://rdrr.io/r/stats/residuals.html),
[`stats::weights()`](https://rdrr.io/r/stats/weights.html) and
[`stats::sigma()`](https://rdrr.io/r/stats/sigma.html) do what they do
for a `glm`, which is also most of what insight needs to make the fit
legible to the easystats packages.

### Leave-One-Out Is Approximate, and Strained Here

[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) estimates the
leave-one-out predictive density by importance sampling from the
full-data posterior, and the estimate is trustworthy only when the
importance weights have a finite variance, which is what the Pareto
\\k\\ diagnostic reports on. A forest is a very flexible function of the
predictors, so a single observation can have a lot of influence on the
leaves it lands in, and high \\k\\ values are common rather than
exceptional. Note that the warning loo prints in that case is not
boilerplate; it says the number is not reliable, and held-out data are
the alternative. A log score on data the model has not seen is available
directly:

    predict(fit, newdata = held_out, type = "density", log = TRUE)

### What a Posterior Predictive Draw Is On

The replicate outcomes are on the scale the likelihood was written on,
which is the scale
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
stored the response on. A **binomial** response is a proportion, so
binary data come back as 0 and 1, and data given as two columns or with
prior weights come back as a fraction of the trials. A response with
**categories** comes back as an integer category index, from 1 to the
number of categories, because a matrix cannot hold a factor;
`fit$levels` names them, and
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) returns
factors instead, since its result is a data frame and can. An
**accelerated failure time** response comes back as a time rather than a
log time, and it is an event time: the predictive distribution of the
outcome does not know about the censoring that may have hidden it, so
comparing replicates against censored observations is not like for like,
and `pp_check()` says so. A
[`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
fit has **no posterior predictive distribution at all**, because a log
density supplies no way to draw from it, so those methods error.

### What Is Deliberately Absent

There is no [`logLik()`](https://rdrr.io/r/stats/logLik.html) method,
and that is a choice rather than a gap. The generic exists so that
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) can be computed, and
both need a count of parameters, which a forest does not have, since the
number of leaves is itself drawn from the posterior.
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) and
[`loo::waic()`](https://mc-stan.org/loo/reference/waic.html) are the
corresponding quantities for a model like this one, and they are
computed from the posterior rather than from a parameter count.

For the same reason
[`performance::check_normality()`](https://easystats.github.io/performance/reference/check_normality.html)
and
[`performance::check_outliers()`](https://easystats.github.io/performance/reference/check_outliers.html)
do not work: they ask for a likelihood-ratio test and for Cook's
distance, neither of which is defined here.
[`performance::check_predictions()`](https://easystats.github.io/performance/reference/check_predictions.html)
does work, through
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html).

### The Bayesian R-Squared

[`performance::r2()`](https://easystats.github.io/performance/reference/r2.html)
returns the quantity of Gelman et al. (2019): per draw, the variance of
the fitted means across observations divided by that variance plus the
variance of the residuals. Being a per-draw quantity it has a posterior,
which is why it is reported with an interval and why it can fall as the
model is made more flexible. It needs a mean, so it is not available for
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
or
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md).

## References

Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019). R-squared for
Bayesian regression models. *The American Statistician*, 73(3), 307–309.

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27(5), 1413–1432.

## See also

[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
for the predictions these methods are built on;
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
for the convergence and mixing diagnostics;
[bartisan-marginaleffects](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)
for reading effects off a fit;
[`vignette("diagnostics")`](https://ngreifer.github.io/bartisan/articles/diagnostics.md)
for the fuller treatment

## Examples

``` r
data("rhc")
set.seed(123)

fit <- bartisan(death ~ . - days, data = rhc, num_trees = 10,
                num_burn = 50, num_draws = 50, chains = 2, verbose = FALSE)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

# Replicate outcomes, one per draw per observation, whose mean is what
# `fitted()` reports
yrep <- rstantools::posterior_predict(fit)
range(colMeans(rstantools::posterior_epred(fit)) - fitted(fit))
#> [1] 0 0

# Pointwise log likelihood, and the fit statistics built on it
loo::waic(rstantools::log_lik(fit))
#> 
#> Computed from 100 by 1500 log-likelihood matrix.
#> 
#>           Estimate   SE
#> elpd_waic   -849.2 17.1
#> p_waic        24.7  0.7
#> waic        1698.5 34.1

# Whether replicate outcomes look like the observed ones
if (rlang::is_installed("bayesplot")) {
  bayesplot::pp_check(fit, type = "bars")
}


# The scalar parameters and a spread of the predictor, as a draws array
if (rlang::is_installed("posterior")) {
  posterior::summarise_draws(posterior::as_draws(fit))
}
#> # A tibble: 12 × 10
#>    variable         mean   median     sd    mad       q5      q95  rhat ess_bulk
#>    <chr>           <dbl>    <dbl>  <dbl>  <dbl>    <dbl>    <dbl> <dbl>    <dbl>
#>  1 loglik       -8.37e+2 -8.36e+2 5.18   5.93   -845.    -828.    1.29      6.14
#>  2 sigma_mu.eta  4.67e-1  4.72e-1 0.0946 0.0901    0.343    0.659 1.13     13.7 
#>  3 eta[119]     -1.33e+0 -1.35e+0 0.286  0.231    -1.80    -0.844 1.04     51.1 
#>  4 eta[24]      -3.87e-1 -3.51e-1 0.434  0.469    -1.13     0.207 1.54      4.34
#>  5 eta[648]     -9.89e-4  1.71e-2 0.352  0.389    -0.591    0.488 1.16     14.1 
#>  6 eta[1432]     2.74e-1  2.58e-1 0.264  0.266    -0.195    0.686 1.09     18.1 
#>  7 eta[975]      5.48e-1  5.83e-1 0.339  0.252    -0.156    1.04  1.24      7.41
#>  8 eta[1083]     8.23e-1  8.27e-1 0.345  0.313     0.302    1.36  0.995    36.4 
#>  9 eta[144]      1.28e+0  1.29e+0 0.285  0.329     0.848    1.73  1.21      7.97
#> 10 eta[347]      1.65e+0  1.65e+0 0.355  0.392     1.14     2.19  1.40      5.30
#> 11 eta[513]      2.00e+0  2.02e+0 0.471  0.453     1.16     2.77  1.46      4.66
#> 12 eta[1135]     2.99e+0  2.97e+0 0.422  0.467     2.38     3.71  1.39      5.18
#> # ℹ 1 more variable: ess_tail <dbl>
```
