
``` r
variable_importance(fit)
#> Variable importance
#> 
#>  variable prop_used prop_splits splits
#>    surv2m     1.000       0.116    8.5
#>       age     1.000       0.089    6.5
#>     paco2     1.000       0.073    5.4
#>       rhc     1.000       0.070    5.1
#>      pafi     1.000       0.069    5.1
#>       aps     0.999       0.082    6.1
#>       edu     0.998       0.062    4.6
#>      hema     0.996       0.065    4.8
#>      card     0.996       0.060    4.4
#>      resp     0.990       0.064    4.7
#>       sex     0.988       0.068    5.0
#>    meanbp     0.988       0.063    4.7
#>      race     0.986       0.059    4.4
#>      crea     0.980       0.060    4.5
#> 
#> ℹ Fitted with `sparsity = FALSE` (the default), so every predictor keeps a
#>   share of the rules and prop_used is near 1 throughout. Refit with `sparsity =
#>   TRUE` to read it as a selection rule.
#> ℹ splits_lower and splits_upper hold the 95% interval, not shown above.
```