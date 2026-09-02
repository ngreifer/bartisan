# Interfaces to other packages

Methods that let a `bartisan` fit be used by the packages that assess
model fit, rather than requiring the posterior draws to be pulled out
and handled by hand. There is nothing to set up: load the other package
and call its function on the fit.

- **Posterior predictions.**
  [`rstantools::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  draws replicate outcomes from the fitted model,
  [`rstantools::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  gives their mean and
  [`rstantools::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html)
  the additive predictor, following the rstantools conventions that brms
  and rstanarm follow.
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) is the
  same thing in the shape base R expects.

- **Pointwise likelihood.**
  [`rstantools::log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
  returns the draws-by-observations matrix of log-likelihood
  contributions, which is what
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) and
  [`loo::waic()`](https://mc-stan.org/loo/reference/waic.html) need;
  both have methods here.

- **Graphical checks.** `pp_check()` runs any of the bayesplot
  posterior-predictive checks on the fit.

- **Summaries.**
  [`performance::model_performance()`](https://easystats.github.io/performance/reference/model_performance.html)
  collects the fit statistics in one table,
  [`performance::r2()`](https://easystats.github.io/performance/reference/r2.html)
  gives the Bayesian \\R^2\\, and
  [`posterior::as_draws()`](https://mc-stan.org/posterior/reference/draws.html)
  hands the scalar parameters to
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
  or to the bayesplot MCMC diagnostics.

- **Basic accessors.**
  [`stats::fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`stats::residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`stats::weights()`](https://rdrr.io/r/stats/weights.html) and
  [`stats::sigma()`](https://rdrr.io/r/stats/sigma.html) do what they do
  for a `glm`, which is also most of what insight needs to make the fit
  legible to the easystats packages.

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

  a fitted model from
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- newdata:

  optional data frame at which to evaluate, as in
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md).

- iterations:

  optional integer vector selecting which stored draws to use. Defaults
  to all of them.

- offset, weights:

  an offset and prior weights for `newdata`, as in
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md).
  For a binomial response the weights are the numbers of trials, and so
  are what a replicate outcome is a fraction of.

- ...:

  further arguments, passed on where the method has somewhere to pass
  them and ignored otherwise.

- transform:

  for `posterior_linpred()`, return the predictor mapped through the
  inverse link, which is what `posterior_epred()` does.

- nsim, ndraws:

  the number of posterior draws to use, chosen at random from the
  retained ones.

- seed:

  optional seed, set with
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before drawing and
  restored afterwards, following the
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html)
  convention.

- type:

  for [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), the
  prediction scale, passed to
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md).
  For `pp_check()`, the name of the bayesplot check to run without its
  `ppc_` prefix, so that `"dens_overlay"` calls
  [`bayesplot::ppc_dens_overlay()`](https://mc-stan.org/bayesplot/reference/PPC-distributions.html)
  .

- eta:

  for `as_draws()`, which columns of the additive predictor to carry
  into the draws array alongside the scalar parameters. `TRUE`, the
  default, takes a representative ten spread across the range of the
  fitted function – there is one column per observation, and an array
  with thousands of them is not something
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
  or a trace plot can be pointed at. `FALSE` takes none. A numeric
  vector takes those observations. The predictor is the quantity whose
  convergence usually matters, and the one `fit$rhat` reports on, so it
  is included by default.

- verbose:

  whether to report problems that do not stop the computation.

- metrics:

  for `model_performance()`, `"all"` or a character vector selecting
  from `"ELPD"`, `"LOOIC"`, `"WAIC"`, `"R2"`, `"RMSE"` and `"SIGMA"`.

## Value

`posterior_predict()`, `posterior_epred()`, `posterior_linpred()` and
`log_lik()` return a matrix of draws by observations.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) returns a data
frame of one column per replicate. `loo()` and `waic()` return the
objects those functions return. `model_performance()` returns a one-row
data frame. `as_draws()` returns a `draws_array` of iterations by chains
by parameters. The accessors return what their names suggest.

## Leave-one-out is approximate, and the approximation is strained here

[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) estimates the
leave-one-out predictive density by importance sampling from the
full-data posterior, and the estimate is trustworthy only when the
importance weights have a finite variance – which is what the Pareto
\\k\\ diagnostic reports on. A forest is a very flexible function of the
predictors, so a single observation can have a lot of influence on the
leaves it lands in, and high \\k\\ values are common rather than
exceptional. The warning loo prints in that case is not boilerplate;
treat it as saying that the number is not reliable, and reach for
held-out data instead. A log score on data the model has not seen is
available directly:

    predict(fit, newdata = held_out, type = "density", log = TRUE)

## What a posterior predictive draw is on

The replicate outcomes are on the scale the likelihood was written on,
which is the scale
[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
stored the response on:

- A binomial response is a **proportion**, so binary data come back as 0
  and 1, and data given as two columns or with prior weights come back
  as a fraction of the trials.

- A response with categories comes back as an **integer category
  index**, from 1 to the number of categories, because a matrix cannot
  hold a factor. `fit$levels` names them.
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) returns
  factors instead, since its result is a data frame and can.

- An accelerated failure time response comes back as a **time**, not a
  log time, and it is an event time: the predictive distribution of the
  outcome does not know about the censoring that may have hidden it.
  Comparing replicates against censored observations is therefore not
  like for like, and `pp_check()` says so.

- A
  [`custom_family()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
  fit has no posterior predictive distribution at all, because a log
  density supplies no way to draw from it. Those methods error.

## What is deliberately absent

There is no [`logLik()`](https://rdrr.io/r/stats/logLik.html) method,
and that is a choice rather than a gap. The generic exists so that
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) can be computed, and
both need a count of parameters – which a forest does not have, since
the number of leaves is itself drawn from the posterior.
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

## The Bayesian R-squared

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

[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md),
[bartisan-marginaleffects](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md)

## Examples

``` r
set.seed(1)
n <- 200
d <- data.frame(x1 = runif(n), x2 = runif(n))
d$y <- rpois(n, exp(0.5 + d$x1))

fit <- bartisan(y ~ x1 + x2, data = d, family = poisson(),
               control = bartisan_control(num_trees = 10, num_burn = 50,
                                         num_draws = 50, verbose = FALSE))

# Replicate outcomes, one per draw per observation
yrep <- rstantools::posterior_predict(fit)
dim(yrep)
#> [1]  50 200

# Their mean, which is what fitted() reports
range(colMeans(rstantools::posterior_epred(fit)) - fitted(fit))
#> [1] 0 0

# Pointwise log likelihood, and the fit statistics built on it
log_likelihood <- rstantools::log_lik(fit)
dim(log_likelihood)
#> [1]  50 200
loo::waic(log_likelihood)
#> 
#> Computed from 50 by 200 log-likelihood matrix.
#> 
#>           Estimate   SE
#> elpd_waic   -382.5  8.8
#> p_waic         5.6  0.6
#> waic         765.0 17.7
```
