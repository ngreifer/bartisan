# Varying coefficients

The coefficient functions of a model fitted with
[`vc()`](https://ngreifer.github.io/bartisan/reference/vc.md) terms,
evaluated at each observation. A forest has no coefficient vector, so
for any other model this returns nothing; for a varying-coefficient
model the coefficients are functions and this is what they come to.

## Usage

``` r
# S3 method for class 'bartisan_fit'
coef(object, newdata = NULL, draws = FALSE, ...)
```

## Arguments

- object:

  a fitted
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)
  model.

- newdata:

  optional data to evaluate the coefficients at. The default, `NULL`,
  uses the data the model was fitted to.

- draws:

  `FALSE`, the default, returns the posterior mean of each coefficient
  at each observation. `TRUE` returns every draw, as a list of
  draws-by-observations matrices, one per coefficient.

- ...:

  ignored.

## Value

With `draws = FALSE`, a matrix with one row per observation and one
column per coefficient. With `draws = TRUE`, a named list of matrices.

## Details

The control function is not among them. It is the surface at the value
each covariate was centered on, which is a prediction rather than a
coefficient; `predict(object)` is what reports predictions.

For a factor the coefficients are recentered to sum to zero across its
levels, which is what makes them the deviations they are reported as.
The symmetric coding carries one spare function-valued dimension, so
this is exact rather than an approximation, and it is the reason a
factor's reference level is a choice made here rather than at fitting
time.

## Examples

``` r
set.seed(1)
n <- 200
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z = rbinom(n, 1, 0.5))
d$y <- d$x1 + d$z * (1 + d$x2) + rnorm(n)

fit <- bartisan(y ~ x1 + x2 + vc(z), data = d, family = gaussian(),
                control = bartisan_control(num_trees = 10, num_burn = 50,
                                           num_draws = 50, verbose = FALSE))

head(coef(fit))
#>               z
#> [1,]  1.7620856
#> [2,]  2.4075532
#> [3,]  2.3611641
#> [4,]  0.4148748
#> [5,] -0.2074336
#> [6,]  2.5594020
```
