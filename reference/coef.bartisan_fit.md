# Varying coefficients

The coefficient functions of a model fitted with
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) terms,
evaluated at each observation. A forest has no coefficient vector, so a
fit with no
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) term has
no coefficient to report and this errors rather than returning anything.

## Usage

``` r
# S3 method for class 'bartisan_fit'
coef(object, newdata = NULL, draws = FALSE, ...)
```

## Arguments

- object:

  a `<bartisan_fit>` object; the output of a call to
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md),
  fitted with at least one
  [`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) term.

- newdata:

  optional; a data frame at which to evaluate the coefficients. Default
  is `NULL` to use the data the model was fitted to.

- draws:

  `logical`; whether to return every posterior draw of each coefficient
  rather than its posterior mean at each observation. Default is `FALSE`
  to return the posterior means.

- ...:

  not used.

## Value

With `draws = FALSE`, a matrix with one row per observation and one
column per coefficient. With `draws = TRUE`, a named list of
draws-by-observations matrices, one per coefficient.

## Details

The coefficients are the varying ones alone. The control function is the
surface at the value each covariate was centered on, which is a
prediction rather than a coefficient, and `predict(object)` reports it.

For a factor the coefficients are recentered to sum to zero across its
levels, which makes them the deviations they are reported as. The
symmetric coding carries one spare function-valued dimension, so this is
exact rather than an approximation, and it is the reason a factor's
reference level is a choice made here rather than at fitting time.

## See also

[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) for
declaring a varying coefficient;
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
for which predictors a forest uses at all

## Examples

``` r
data("rhc")
set.seed(123)

# The effect of catheterization is allowed to vary with
# the other predictors, so its coefficient is a function
# rather than a number
fit <- bartisan(death ~ age + aps + surv2m + vc(rhc),
                data = rhc, num_trees = 10, num_burn = 50,
                num_draws = 50, verbose = FALSE)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` explicitly to silence this message.

# One coefficient per patient, on the link scale
head(coef(fit))
#>               rhc
#> [1,]  0.370664770
#> [2,] -0.005113301
#> [3,]  0.129570231
#> [4,]  0.307429405
#> [5,]  0.269363681
#> [6,]  0.186869057

# How much it varies across patients
quantile(coef(fit)[, "rhc"])
#>           0%          25%          50%          75%         100% 
#> -0.005113301  0.171367444  0.271716353  0.395163994  0.671462755 
```
