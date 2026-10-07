
``` r
variable_importance(fit_aps)
#> Variable importance
#> 
#>    predictor variable prop_used prop_splits splits
#>          rhc      aps     1.000       0.500   37.7
#>          rhc      age     1.000       0.500   37.6
#>  (Intercept)   surv2m     1.000       0.122    9.2
#>  (Intercept)      age     1.000       0.096    7.2
#>  (Intercept)    paco2     1.000       0.081    6.1
#>  (Intercept)      aps     0.999       0.082    6.2
#>  (Intercept)     hema     0.999       0.073    5.5
#>  (Intercept)      edu     0.999       0.069    5.2
#>  (Intercept)     pafi     0.998       0.074    5.6
#>  (Intercept)   meanbp     0.996       0.073    5.5
#>  (Intercept)     crea     0.996       0.070    5.3
#>  (Intercept)      sex     0.996       0.065    4.9
#>  (Intercept)     card     0.994       0.066    5.0
#>  (Intercept)     resp     0.993       0.064    4.8
#>  (Intercept)     race     0.990       0.064    4.8
#>          rhc      sex     0.000       0.000    0.0
#>          rhc     race     0.000       0.000    0.0
#>          rhc      edu     0.000       0.000    0.0
#>          rhc   meanbp     0.000       0.000    0.0
#>          rhc     resp     0.000       0.000    0.0
#>          rhc     hema     0.000       0.000    0.0
#>          rhc     pafi     0.000       0.000    0.0
#>          rhc    paco2     0.000       0.000    0.0
#>          rhc     crea     0.000       0.000    0.0
#>          rhc   surv2m     0.000       0.000    0.0
#>          rhc     card     0.000       0.000    0.0
#> 
#> ℹ Fitted with `sparsity = FALSE` (the default), so every predictor keeps a
#>   share of the rules and prop_used is near 1 throughout. Refit with `sparsity =
#>   TRUE` to read it as a selection rule.
#> ℹ splits_lower and splits_upper hold the 95% interval, not shown above.
```