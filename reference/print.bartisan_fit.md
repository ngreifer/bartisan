# Summarize a generalized BART model

[`print()`](https://rdrr.io/r/base/print.html) reports what was fit and
how long the chain is.
[`summary()`](https://rdrr.io/r/base/summary.html) adds posterior
summaries of the nuisance parameters and of how often each predictor was
used in a splitting rule, which is the model's variable selection
output.

## Usage

``` r
# S3 method for class 'bartisan_fit'
print(x, digits = 3L, ...)

# S3 method for class 'bartisan_fit'
summary(object, level = 0.95, ...)

# S3 method for class 'summary.bartisan_fit'
print(x, digits = 3, ...)
```

## Arguments

- x, object:

  a fitted model from
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).

- digits:

  number of significant digits to print.

- ...:

  ignored, present for compatibility with the generics.

- level:

  width of the reported posterior intervals.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument
invisibly. [`summary()`](https://rdrr.io/r/base/summary.html) returns a
list of class `summary.bartisan_fit`, with a `usage` element giving the
posterior summary of the splitting counts for each predictor group, and
an `aux` element for the nuisance parameters when the family has any.

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md)

## Examples

``` r
set.seed(1)

n <- 150
d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
d$y <- 2 * d$x1 + rnorm(n, sd = 0.3)

fit <- bartisan(y ~ x1 + x2 + x3, data = d, family = gaussian(),
               control = bartisan_control(num_trees = 10, num_burn = 50,
                                         num_draws = 50, verbose = FALSE))
fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = y ~ x1 + x2 + x3, data = d, family = gaussian(), 
#>     control = bartisan_control(num_trees = 10, num_burn = 50, 
#>         num_draws = 50, verbose = FALSE))
#> 
#> Family: "gaussian" with the "identity" link
#> Observations: 150
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 50 kept after 50 warmup
#> 
#> Posterior means: sigma = 0.293
summary(fit)
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = y ~ x1 + x2 + x3, data = d, family = gaussian(), 
#>     control = bartisan_control(num_trees = 10, num_burn = 50, 
#>         num_draws = 50, verbose = FALSE))
#> 
#> Family: "gaussian" with the "identity" link
#> Observations: 150
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 50
#> 
#> Nuisance parameters
#>        mean    sd lower upper
#> sigma 0.293 0.018 0.269 0.336
#> 
#> Predictor usage
#> Splitting rules per draw, and how often used at all.
#>    mean    sd lower upper prop_used
#> x1 13.9 2.613    10    18         1
#> x2  0.0 0.000     0     0         0
#> x3  0.0 0.000     0     0         0
```
