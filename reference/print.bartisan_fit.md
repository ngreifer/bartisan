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
#> age    2.50 1.202     1 5.000      1.00
#> surv2m 3.36 1.133     2 6.000      1.00
#> pafi   1.22 0.504     1 2.000      0.98
#> rhc    1.29 0.556     0 2.000      0.96
#> edu    2.63 1.709     0 5.525      0.93
#> paco2  1.43 0.902     0 3.000      0.87
#> card   0.96 1.063     0 3.000      0.56
#> crea   0.75 0.845     0 2.525      0.52
#> resp   0.71 0.946     0 3.000      0.50
#> aps    0.70 0.847     0 2.525      0.48
#> hema   0.45 0.716     0 2.000      0.33
#> meanbp 0.28 0.514     0 1.525      0.25
#> race   0.29 0.574     0 2.000      0.23
#> sex    0.17 0.378     0 1.000      0.17

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
#>       mean    sd lower upper
#> sigma 1.53 0.025   1.5 1.558
#> 
#> Predictor usage
#> Splitting rules per draw, and how often used at all.
#>        mean    sd lower upper prop_used
#> rhc    2.08 1.066   1.0   4.0      1.00
#> race   4.52 0.953   4.0   6.0      1.00
#> aps    3.44 1.296   1.9   5.0      1.00
#> surv2m 4.44 1.128   3.0   6.0      1.00
#> age    0.84 1.037   0.0   2.1      0.52
#> meanbp 0.38 0.667   0.0   1.0      0.30
#> paco2  0.46 0.930   0.0   1.1      0.28
#> hema   0.34 0.772   0.0   2.0      0.18
#> pafi   0.20 0.452   0.0   1.0      0.18
#> resp   0.20 0.571   0.0   1.0      0.14
#> crea   0.10 0.303   0.0   0.1      0.10
#> edu    0.02 0.141   0.0   0.0      0.02
#> card   0.02 0.141   0.0   0.0      0.02
#> sex    0.00 0.000   0.0   0.0      0.00
```
