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

# Whether a patient died, with catheterization among the
# predictors
fit <- bartisan(death ~ rhc + age + sex + race + edu +
                  aps + meanbp + resp + hema + pafi +
                  paco2 + crea + surv2m + card,
                data = rhc, num_trees = 10, num_burn = 50,
                num_draws = 50, chains = 2, verbose = FALSE)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` explicitly to silence this message.

# What was fit, and how many draws it rests on
fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ rhc + age + sex + race + edu + aps + 
#>     meanbp + resp + hema + pafi + paco2 + crea + surv2m + card, 
#>     data = rhc, num_trees = 10, num_burn = 50, num_draws = 50, 
#>     chains = 2, verbose = FALSE)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 100 kept across 2 chains after 50 warmup

# The splitting counts, which say which predictors the
# forest reaches for
summary(fit)
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ rhc + age + sex + race + edu + aps + 
#>     meanbp + resp + hema + pafi + paco2 + crea + surv2m + card, 
#>     data = rhc, num_trees = 10, num_burn = 50, num_draws = 50, 
#>     chains = 2, verbose = FALSE)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 100
#> 
#> Predictor usage
#> Splitting rules per draw, and how often used at all.
#>        mean    sd lower upper prop_used
#> age    1.84 0.929     1 4.000      1.00
#> aps    1.97 0.834     1 4.000      1.00
#> pafi   1.78 0.824     1 3.000      1.00
#> paco2  2.44 0.857     1 4.525      1.00
#> surv2m 2.94 0.908     2 5.000      1.00
#> rhc    1.20 0.725     0 3.000      0.89
#> meanbp 1.13 0.734     0 2.525      0.82
#> edu    1.16 0.992     0 3.525      0.75
#> crea   0.78 0.786     0 2.000      0.58
#> card   0.83 0.900     0 3.000      0.57
#> resp   0.70 0.798     0 2.525      0.54
#> hema   0.52 0.643     0 2.000      0.44
#> sex    0.39 0.665     0 2.000      0.30
#> race   0.20 0.471     0 1.525      0.17

# A family with a nuisance parameter reports its
# posterior too, here the residual standard deviation of
# the log survival time
fit2 <- bartisan(log(days) ~ rhc + age + sex + race +
                   edu + aps + meanbp + resp + hema +
                   pafi + paco2 + crea + surv2m + card,
                 data = rhc, family = gaussian(),
                 num_trees = 10, num_burn = 50,
                 num_draws = 50, verbose = FALSE)

summary(fit2, level = .8)
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = log(days) ~ rhc + age + sex + race + edu + 
#>     aps + meanbp + resp + hema + pafi + paco2 + crea + surv2m + 
#>     card, data = rhc, family = gaussian(), num_trees = 10, num_burn = 50, 
#>     num_draws = 50, verbose = FALSE)
#> 
#> Family: "gaussian" with the "identity" link
#> Observations: 1500
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 50
#> 
#> Nuisance parameters
#>       mean    sd lower upper
#> sigma 1.54 0.028 1.505 1.573
#> 
#> Predictor usage
#> Splitting rules per draw, and how often used at all.
#>        mean    sd lower upper prop_used
#> rhc    1.40 0.670   1.0     2      1.00
#> aps    1.24 0.431   1.0     2      1.00
#> surv2m 2.44 1.033   1.0     4      1.00
#> paco2  1.38 0.780   0.9     2      0.90
#> meanbp 1.04 0.755   0.0     2      0.78
#> resp   1.64 1.102   0.0     3      0.78
#> edu    0.92 0.944   0.0     2      0.58
#> race   0.46 0.613   0.0     1      0.40
#> crea   0.54 0.788   0.0     2      0.38
#> sex    0.42 0.642   0.0     1      0.34
#> card   0.24 0.476   0.0     1      0.22
#> hema   0.18 0.388   0.0     1      0.18
#> age    0.08 0.274   0.0     0      0.08
#> pafi   0.08 0.274   0.0     0      0.08
```
