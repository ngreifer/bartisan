# Summarize a generalized BART model

[`print()`](https://rdrr.io/r/base/print.html) reports what was fit and
how long the chain is.
[`summary()`](https://rdrr.io/r/base/summary.html) reports what
[`print()`](https://rdrr.io/r/base/print.html) does not: posterior
summaries of the nuisance parameters and the random-effect scales, a
quick check of convergence, and the most used predictors, with pointers
to the functions that examine each of these in full.

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
random-effect term, the splitting counts of each predictor group, the
output of
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md),
and R-hat and the effective sample sizes of the log likelihood.

## Details

The printed summary is a starting point for the functions that examine a
fit in full, and it computes only what is cheap. Its convergence line
gives R-hat and the bulk and tail effective sample sizes of the log
likelihood, which take milliseconds because the log likelihood is a
single series. It can look fine while individual fitted values mix
badly, which is what
[`diagnose()`](https://ngreifer.github.io/bartisan/reference/diagnose.md)
checks, at a cost that grows with the number of observations. The
variable importance table shows at most the five most used predictors of
the first forest;
[`variable_importance()`](https://ngreifer.github.io/bartisan/reference/variable_importance.md)
gives every predictor in every forest. The summary also points to
[loo()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
and
[kfold()](https://ngreifer.github.io/bartisan/reference/bartisan-interop.md)
for comparing fits and to
[`partial_dependence()`](https://ngreifer.github.io/bartisan/reference/partial_dependence.md)
for how the predictions depend on a predictor, none of which it runs,
since each can take seconds or more on a large fit.

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
#> Convergence and mixing
#> 
#> Log likelihood: R-hat 1.058, bulk ESS 38, tail ESS 81, over 2 chains
#> 
#> ℹ Use diagnose() (`?bartisan::diagnose`) to examine convergence and mixing
#>   diagnostics.
#> 
#> Variable importance
#> 
#> Predictors ranked by use; 5 of 14 shown.
#>  variable prop_used prop_splits splits
#>    surv2m         1       0.167    2.9
#>     paco2         1       0.138    2.4
#>       aps         1       0.111    2.0
#>       age         1       0.103    1.8
#>      pafi         1       0.098    1.8
#> 
#> ℹ Use variable_importance() (`?bartisan::variable_importance`) to examine
#>   variable importance.
#> 
#> Further tools
#> 
#> ℹ Use loo() (`?bartisan::loo.bartisan_fit`) to compare this fit with others, or
#>   kfold() (`?bartisan::kfold.bartisan_fit`) if `loo()` reports many Pareto k
#>   values above 0.7.
#> ℹ Use partial_dependence() (`?bartisan::partial_dependence`) and plot()
#>   (`?bartisan::plot.bartisan_fit`) to view the partial dependence of the
#>   predictions on a predictor.
#> 

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
#> Nuisance parameters
#>       mean    sd lower upper
#> sigma 1.54 0.028 1.505 1.573
#> 
#> Convergence and mixing
#> 
#> Log likelihood: R-hat 1.427, bulk ESS 6, tail ESS 14, over 1 chain
#> 
#> ℹ Use diagnose() (`?bartisan::diagnose`) to examine convergence and mixing
#>   diagnostics.
#> 
#> Variable importance
#> 
#> Predictors ranked by use; 5 of 14 shown.
#>  variable prop_used prop_splits splits
#>    surv2m      1.00       0.205    2.4
#>       rhc      1.00       0.115    1.4
#>       aps      1.00       0.106    1.2
#>     paco2      0.90       0.115    1.4
#>      resp      0.78       0.138    1.6
#> 
#> ℹ Use variable_importance() (`?bartisan::variable_importance`) to examine
#>   variable importance.
#> 
#> Further tools
#> 
#> ℹ Use loo() (`?bartisan::loo.bartisan_fit`) to compare this fit with others, or
#>   kfold() (`?bartisan::kfold.bartisan_fit`) if `loo()` reports many Pareto k
#>   values above 0.7.
#> ℹ Use partial_dependence() (`?bartisan::partial_dependence`) and plot()
#>   (`?bartisan::plot.bartisan_fit`) to view the partial dependence of the
#>   predictions on a predictor.
#> 
```
