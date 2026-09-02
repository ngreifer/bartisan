# Counterfactual estimands with marginaleffects

A `bartisan` fit works with the marginaleffects package, so that
predictions, comparisons and slopes – and hypothesis tests on any of
them – can be computed without extracting draws by hand. Load
marginaleffects and call its functions on the fit directly; there is
nothing to set up.

## Usage

``` r
# S3 method for class 'bartisan_fit'
formula(x, ...)

# S3 method for class 'bartisan_fit'
terms(x, ...)

# S3 method for class 'bartisan_fit'
model.frame(formula, ...)

# S3 method for class 'bartisan_fit'
nobs(object, ...)

# S3 method for class 'bartisan_fit'
family(object, ...)

# S3 method for class 'bartisan_fit'
get_predict(model, newdata = NULL, type = NULL, ...)

# S3 method for class 'bartisan_fit'
get_group_names(model, ...)

# S3 method for class 'bartisan_fit'
get_coef(model, ...)

# S3 method for class 'bartisan_fit'
set_coef(model, coefs, ...)

# S3 method for class 'bartisan_fit'
get_vcov(model, ...)

# S3 method for class 'bartisan_fit'
get_data(x, ...)
```

## Arguments

- ...:

  further arguments. `values`, `iterations`, `offset`, `weights` and
  `log` are passed on to
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md);
  anything else is ignored, since marginaleffects puts arguments of its
  own here too. marginaleffects warns that it does not recognize
  `values`, which is expected – it is this package's argument, not one
  of its own – and the value is used regardless.

- model, x, object, formula:

  a fitted `bartisan` object. The last is named that only because
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  names its first argument that way.

- newdata:

  data at which to evaluate the fit. Defaults to the data the model was
  fitted to, which is retained in the fit for this purpose.

- type:

  the scale to work on. `"response"` is the fitted mean, `"link"` the
  additive predictor, `"prob"` the per-category probabilities of a
  categorical family, `"mean"` the mean of a categorical response with
  its labels read as numbers, and `"stdlv"` the standardized latent
  variable of an ordinal fit; the last three are described under
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md).
  A categorical family has no mean, so `"response"` gives `"prob"`
  there. `"probs"`, `"lp"` and `"lv"` are accepted as aliases for
  `"prob"`, `"link"` and `"link"`, since those are the names the same
  quantities go by for other ordinal fits in marginaleffects.

- coefs:

  ignored; a forest has no coefficient vector.

## Value

[`get_predict()`](https://rdrr.io/pkg/marginaleffects/man/get_predict.html)
returns a data frame with columns `rowid`, `group` and `estimate`,
carrying the draws in a `"posterior_draws"` attribute of observations by
draws, which is the interface marginaleffects documents. The other
methods exist to satisfy the generic and return what their names
suggest.

## How the uncertainty is computed

A forest has no coefficient vector and no variance-covariance matrix, so
the delta method marginaleffects uses for a frequentist model has
nothing to work with. It has something better here: the posterior draws.
Every estimand is computed by pushing all of the draws through the same
transformation and summarizing at the end, so an interval is a posterior
quantile rather than a normal approximation, and a nonlinear estimand
needs no approximation at all. This is the same path marginaleffects
takes for brms and rstanarm fits.

One consequence worth knowing: marginaleffects centers a posterior at
its **median**, where
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
reports its **mean**. The two are summarizing the same draws, so a
difference between them is the skewness of the posterior and not a
disagreement.

## A contrast of exactly zero is usually real

[`avg_comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
reporting an estimate of exactly `0` is the most common surprise here,
and it is neither package computing anything wrong. Two facts meet to
produce it.

The posterior of a contrast has an **atom at exactly zero**. In any draw
where no tree in the forest splits on the variable being contrasted, the
fit does not depend on that variable at all, so the two counterfactual
predictions are identical to the last bit and their difference is
exactly zero. That is not a near-zero value that rounding flattered; it
is a point mass. The Dirichlet sparsity prior on the splitting
proportions, `update_s` in
[`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md),
is what makes those draws common: it is a variable selection prior, and
dropping a weak predictor from every tree is what it is for.

And marginaleffects centers a posterior at its **median**. So once the
atom holds more than half the mass, the reported estimate is exactly
zero however large the rest of the posterior is. On `MatchIt::lalonde`
with the default settings, `treat` was absent from all 50 trees in 64%
of draws and the contrast came out exactly zero in 65%, which put the
median at 0 while the posterior mean was 197 and the upper limit was
above 2000.

Four things to do about it, in the order worth trying:

1.  **Look at the inclusion probability**, which is what the zero is
    telling you. `summary(fit)` reports it as `prop_used`: the posterior
    probability that each predictor group appears anywhere in the
    forest. A contrast whose median is zero is a predictor the model is
    not sure belongs.

2.  **Ask for the mean instead**, with
    `options(marginaleffects_posterior_center = mean)`. The mean is the
    summary
    [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
    reports, and it is the one that behaves sensibly against an atom.

3.  **Reconsider the sparsity prior** if variable selection is not what
    you want from the fit. `sparsity = FALSE` in
    [`bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.md)
    removes the atom almost entirely: on the example above it fell from
    20% of draws to none by 50 trees. A larger `num_trees` does *not*
    remove it, which is worth knowing because it looks as though it
    should – with the prior on, the contrast was exactly zero in 20% of
    draws at 50 trees and 18% at 200. Turning the prior off is a
    modeling choice rather than a fix, so make it for a reason: it is
    the right one when a contrast on a particular predictor is the
    estimand, and the wrong one when there are many predictors and most
    are irrelevant.

4.  **Run several chains and compare them.** The variable selection
    state mixes slowly, because a predictor whose splitting proportion
    has gone small is rarely proposed and so is hard to get back in. On
    the example above, four chains disagreed by more than 100% of the
    estimate at every tree count from 20 to 200 with the prior on; with
    `sparsity = FALSE` and 50 trees they agreed to within 9%. A single
    chain can look much more settled than the posterior is.

## Slopes need a linear predictor transform

A slope is a numerical derivative, and taking one requires the fitted
function to be differentiable in the predictor *as the caller supplies
it*. The default `x_transform = "quantile"` maps each predictor through
its empirical distribution function before any rule sees it, and an
empirical distribution function is a step function – so the fit is a
step function of the original predictor whatever the decision rules are,
and its difference quotient grows without bound as the step shrinks.
Measured on a smooth surface where the average derivative is zero:

|      |                            |                         |
|------|----------------------------|-------------------------|
| step | `x_transform = "quantile"` | `x_transform = "range"` |
| 1e-4 | -4.79                      | -0.28                   |
| 1e-2 | -0.40                      | -0.30                   |
| 5e-2 | -0.29                      | -0.25                   |

So **use `x_transform = "range"` if slopes are the estimand**, which
maps each predictor linearly and leaves a soft-rule fit differentiable.
Hard rules give a piecewise-constant fit under either transform, and a
derivative of one is not a meaningful quantity however it is computed.

None of this affects
[`marginaleffects::predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html)
or
[`marginaleffects::comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html)
, which evaluate the fit at two points a substantive distance apart
rather than dividing by a vanishing one. Those are the estimands to
reach for with the default transform.

## The usual survival estimand

For a survival family, the estimand is usually a contrast in survival at
a horizon rather than in the predictor. `type = "survival"` with `times`
gives it:

    # The difference in one-year survival between treated and untreated.
    avg_comparisons(fit, variables = "trt", type = "survival", times = 1)

One time per call. marginaleffects checks the dots against a whitelist
of its own, hardcoded per model class, so it warns that it does not
recognize `times` – while passing it through, which is what the warning
says. There is no hook for registering an argument with it, so the
warning is expected and the result is correct.

## What is not covered

`type = "class"` and `type = "density"` are not available, because
neither is one number per observation that an average or a contrast
could be taken of: a class is a factor, and a density needs the outcome,
which a counterfactual grid does not have. Call
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
for those.

`type = "link"` is refused for a family with more than one additive
predictor –
[`location_scale()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md),
the zero-inflated families,
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
– because there is no single link there to be talking about. Those
families work on the response scale, which is one number per observation
whatever the family, and
[`multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
and
[`ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.md)
work on the probability scale, which gives one group per category. To
reach a *particular* predictor of a multi-predictor family, call
[`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md)
directly.

Extrapolation is worth keeping in mind for
[`comparisons()`](https://rdrr.io/pkg/marginaleffects/man/comparisons.html):
a forest is constant outside the range of the predictor it was fitted
on, so a contrast that steps a predictor beyond that range reports the
boundary value rather than an extrapolated one.

## Examples

``` r
set.seed(1)
n <- 200
d <- data.frame(x1 = runif(n), x2 = runif(n))
d$y <- 2 * sin(pi * d$x1) - d$x2 + rnorm(n)

fit <- bartisan(y ~ x1 + x2, data = d, family = gaussian(),
               control = bartisan_control(num_trees = 10, num_burn = 50,
                                         num_draws = 50))

marginaleffects::avg_predictions(fit)
#> 
#>  Estimate 2.5 % 97.5 %
#>     0.939   0.8   1.07
#> 
#> Type: response
#> 
marginaleffects::avg_comparisons(fit)
#> 
#>  Term Estimate  2.5 % 97.5 %
#>    x1   -0.979 -1.576 -0.505
#>    x2   -0.528 -0.976 -0.187
#> 
#> Type: response
#> Comparison: +1
#> 
```
