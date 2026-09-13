# Partial dependence on one or two predictors

Averages the fitted surface over the sample at each value of one or two
predictors, so that what is left is how the prediction moves with them.
`partial_dependence()` returns the values and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws them.

## Usage

``` r
partial_dependence(
  object,
  variables,
  newdata = NULL,
  grid = 25L,
  values = NULL,
  level = 0.95,
  type = "response",
  plot = FALSE
)

# S3 method for class 'bartisan_partial'
print(x, digits = 3L, ...)

# S3 method for class 'bartisan_partial'
plot(x, ...)

# S3 method for class 'bartisan_fit'
plot(x, y, ...)
```

## Arguments

- object:

  a `<bartisan_fit>` object; the output of a call to
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- variables:

  a one-sided formula naming one or two predictors, as in `~ age` or
  `~ age + sex`, or a character vector of their names.

- newdata:

  optional; a data frame to average over. Default is the data the model
  was fit to.

- grid:

  `integer`; how many values of a numeric predictor to evaluate. Default
  is 25. A factor is evaluated at each of its levels whatever this is.

- values:

  optional; a named list giving the values to evaluate a predictor at,
  which overrides `grid` for the predictors it names.

- level:

  `numeric`; the level of the credible interval. Default is `.95`.

- type:

  `string`; the prediction scale, passed to
  [`predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.md).
  Default is `"response"`.

- plot:

  `logical`; whether to draw the result rather than return it. Default
  is `FALSE`. `plot = TRUE` calls `plot.bartisan_partial()`, so the
  argument and the method cannot disagree.

- x:

  a `<bartisan_partial>` object; the output of a call to
  `partial_dependence()`. For `plot.bartisan_fit()`, a `<bartisan_fit>`.

- ...:

  for `plot.bartisan_fit()`, further arguments passed to
  `partial_dependence()`; otherwise ignored.

- y:

  for `plot.bartisan_fit()`, the predictors to plot, as `variables`
  above.

## Value

A `<bartisan_partial>` object, a data frame with one row per grid point,
columns for the predictors, and `estimate`, `lower` and `upper`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
ggplot2 object.

## Details

At each grid value every unit is assigned that value, the prediction is
taken for all of them, and the average over units is taken *within each
posterior draw*. The interval is then a quantile of those averages, so
it is an interval on the average prediction and not on any one unit's.

What this is and is not worth reading as a description of the fit is the
same caveat that applies to any partial dependence plot. Averaging over
the other predictors evaluates the model at covariate combinations that
may not occur, so a curve over a region where the predictor's values are
sparse says more about the prior than about the data. It is a summary of
the fitted function rather than a causal claim; for a contrast that is
meant causally, see
[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md).

## See also

[`estimate_effect()`](https://ngreifer.github.io/bartisan/reference/estimate_effect.md)
for treatment effects rather than fitted surfaces;
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
for which predictors the forest uses;
[bartisan-marginaleffects](https://ngreifer.github.io/bartisan/reference/bartisan-marginaleffects.md),
since
[`marginaleffects::plot_predictions()`](https://rdrr.io/pkg/marginaleffects/man/plot_predictions.html)
draws the same thing with more control over the grid

## Examples

``` r
data("rhc")
set.seed(123)

fit <- bartisan(death ~ age + sex + meanbp + aps, data = rhc,
                num_trees = 10, num_burn = 50, num_draws = 50,
                verbose = FALSE)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

# How the fitted risk moves with mean blood pressure
pd <- partial_dependence(fit, ~ meanbp)
pd
#> Partial dependence
#> 
#> Predictor: "meanbp"
#> Averaged over 1500 units, on the
#>            "response" scale
#> 
#>  meanbp estimate lower upper
#>    0.00    0.655 0.608 0.705
#>    9.25    0.655 0.608 0.705
#>   18.50    0.655 0.608 0.705
#>   27.75    0.655 0.608 0.705
#>   37.00    0.656 0.613 0.705
#>   46.25    0.660 0.627 0.705
#>   55.50    0.655 0.621 0.688
#>   64.75    0.656 0.633 0.674
#>   74.00    0.656 0.633 0.676
#>   83.25    0.655 0.632 0.677
#>   92.50    0.655 0.632 0.678
#>  101.75    0.655 0.631 0.678
#>  111.00    0.655 0.628 0.681
#>  120.25    0.654 0.625 0.685
#>  129.50    0.651 0.625 0.688
#>  138.75    0.649 0.618 0.688
#>  148.00    0.646 0.604 0.695
#>  157.25    0.645 0.600 0.696
#>  166.50    0.645 0.599 0.696
#>  175.75    0.645 0.598 0.696
#>  185.00    0.645 0.597 0.697
#>  194.25    0.645 0.596 0.697
#>  203.50    0.645 0.596 0.697
#>  212.75    0.645 0.596 0.697
#>  222.00    0.645 0.596 0.697
#> 
#> ℹ lower and upper bound the 95% credible interval on the average prediction,
#>   not on any one unit's.

plot(pd)


# The same thing from the fit, which is what the `plot()` method is for
plot(fit, ~ meanbp)


# Two predictors, one of them a factor, which gives a curve per level
plot(fit, ~ meanbp + sex)
```
