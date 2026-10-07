
``` r
variable_importance(fit_vc)
#> Variable importance
#> 
#>    predictor variable prop_used prop_splits splits
#>  (Intercept)   surv2m     1.000       0.117    8.9
#>  (Intercept)      age     1.000       0.097    7.4
#>  (Intercept)    paco2     1.000       0.084    6.4
#>          rhc      aps     1.000       0.082    6.2
#>  (Intercept)     pafi     1.000       0.081    6.1
#>          rhc     race     1.000       0.079    6.0
#>          rhc   meanbp     1.000       0.079    6.0
#>          rhc     resp     1.000       0.078    5.9
#>          rhc     card     1.000       0.077    5.8
#>  (Intercept)     race     1.000       0.067    5.1
#>          rhc     hema     0.999       0.083    6.3
#>  (Intercept)      aps     0.999       0.082    6.3
#>          rhc   surv2m     0.999       0.075    5.7
#>  (Intercept)      edu     0.999       0.073    5.6
#>  (Intercept)     card     0.999       0.069    5.2
#>          rhc      age     0.998       0.077    5.9
#>  (Intercept)     hema     0.998       0.071    5.4
#>  (Intercept)   meanbp     0.998       0.071    5.4
#>  (Intercept)      sex     0.998       0.064    4.9
#>          rhc     pafi     0.996       0.076    5.7
#>          rhc      sex     0.996       0.074    5.6
#>          rhc    paco2     0.996       0.073    5.5
#>          rhc      edu     0.995       0.075    5.7
#>          rhc     crea     0.995       0.073    5.5
#>  (Intercept)     crea     0.989       0.065    4.9
#>  (Intercept)     resp     0.984       0.060    4.5
#> 
#> ℹ Fitted with `sparsity = FALSE` (the default), so every predictor keeps a
#>   share of the rules and prop_used is near 1 throughout. Refit with `sparsity =
#>   TRUE` to read it as a selection rule.
#> ℹ splits_lower and splits_upper hold the 95% interval, not shown above.
```