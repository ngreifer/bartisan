# Summarize a generalized BART model

[`print()`](https://rdrr.io/r/base/print.html) reports what was fit and
how long the chain is.
[`summary()`](https://rdrr.io/r/base/summary.html) adds posterior
summaries of the nuisance parameters, of the random-effect scales, and
of how often each predictor was used in a splitting rule, which is the
model's variable-selection output.

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

  a `<bartisan_fit>` object; the output of a call to
  [`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md).
  For `print.summary.bartisan_fit()`, `x` is instead the
  `<summary.bartisan_fit>` object
  [`summary()`](https://rdrr.io/r/base/summary.html) returned.

- digits:

  `numeric`; the number of significant digits to print. Default is 3.

- ...:

  not used.

- level:

  `numeric`; the width of the posterior intervals
  [`summary()`](https://rdrr.io/r/base/summary.html) reports. Default is
  .95 for 95% intervals.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument
invisibly. [`summary()`](https://rdrr.io/r/base/summary.html) returns a
`<summary.bartisan_fit>` object, a list of the posterior summaries its
own [`print()`](https://rdrr.io/r/base/print.html) method displays: the
nuisance parameters when the family has any, the scale of each
random-effect term, and the splitting counts of each predictor group.

## See also

[`bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.md);
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
for the splitting counts as a data frame rather than printed;
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
for whether the chains the summaries are computed from have converged

## Examples

``` r
data("rhc")
set.seed(123)

# Whether a patient died, with catheterization among the predictors
fit <- bartisan(death ~ . - days, data = rhc, num_trees = 10,
                num_burn = 50, num_draws = 50, chains = 2, verbose = FALSE)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` to choose another, which also silences this message.

# What was fit, and how many draws it rests on
fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ . - days, data = rhc, num_trees = 10, 
#>     num_burn = 50, num_draws = 50, chains = 2, verbose = FALSE)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 100 kept across 2 chains after 50 warmup

# The splitting counts, which say which predictors the forest reaches for
summary(fit)
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ . - days, data = rhc, num_trees = 10, 
#>     num_burn = 50, num_draws = 50, chains = 2, verbose = FALSE)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 100
#> 
#> Predictor usage
#> Splitting rules per draw, and how often used at all.
#>        mean    sd lower upper prop_used
#> age    1.94 0.897     1 4.000      1.00
#> paco2  1.40 0.725     1 3.000      1.00
#> surv2m 2.35 0.626     2 4.000      1.00
#> rhc    1.82 0.770     1 3.000      0.99
#> pafi   1.98 1.189     1 5.000      0.99
#> hema   1.42 1.017     0 3.000      0.78
#> aps    1.11 0.942     0 3.000      0.71
#> resp   0.65 0.757     0 2.000      0.49
#> edu    0.54 0.744     0 2.000      0.40
#> card   0.49 0.674     0 2.000      0.39
#> race   0.41 0.588     0 2.000      0.36
#> meanbp 0.45 0.770     0 3.000      0.32
#> crea   0.23 0.489     0 1.525      0.20
#> sex    0.19 0.465     0 1.525      0.16

# A family with a nuisance parameter reports its posterior too, here the
# residual standard deviation of the log survival time
fit2 <- bartisan(log(days) ~ . - death, data = rhc, family = gaussian(),
                 num_trees = 10, num_burn = 50, num_draws = 50,
                 verbose = FALSE)

summary(fit2, level = .8)
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = log(days) ~ . - death, data = rhc, family = gaussian(), 
#>     num_trees = 10, num_burn = 50, num_draws = 50, verbose = FALSE)
#> 
#> Family: "gaussian" with the "identity" link
#> Observations: 1500
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 50
#> 
#> Nuisance parameters
#>        mean    sd lower upper
#> sigma 1.541 0.028 1.505 1.573
#> 
#> Predictor usage
#> Splitting rules per draw, and how often used at all.
#>        mean    sd lower upper prop_used
#> rhc    3.42 0.992     2   5.0      1.00
#> aps    3.82 1.976     2   7.0      1.00
#> surv2m 5.04 1.603     3   7.0      1.00
#> meanbp 0.78 1.148     0   3.0      0.44
#> paco2  0.56 0.972     0   2.0      0.32
#> race   0.56 0.907     0   2.0      0.30
#> pafi   0.32 0.621     0   1.0      0.24
#> resp   0.22 0.418     0   1.0      0.22
#> hema   0.28 0.607     0   1.0      0.22
#> crea   0.10 0.303     0   0.1      0.10
#> card   0.10 0.303     0   0.1      0.10
#> sex    0.08 0.274     0   0.0      0.08
#> edu    0.04 0.198     0   0.0      0.04
#> age    0.00 0.000     0   0.0      0.00
```
