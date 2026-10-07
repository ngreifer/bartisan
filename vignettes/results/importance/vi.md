
``` r
set.seed(2026)

fit <- bartisan(model, data = rhc, family = binomial(), sparsity = TRUE)

imp <- variable_importance(fit)

imp
#> Variable importance
#> 
#>  variable prop_used prop_splits splits
#>    surv2m     1.000       0.211   15.6
#>       age     1.000       0.182   13.6
#>       rhc     0.995       0.098    7.3
#>      pafi     0.995       0.080    6.0
#>     paco2     0.991       0.105    7.9
#>       aps     0.764       0.062    4.7
#>       edu     0.684       0.053    4.0
#>      resp     0.667       0.059    4.4
#>      race     0.637       0.031    2.3
#>      card     0.621       0.025    1.8
#>      hema     0.574       0.033    2.5
#>    meanbp     0.555       0.022    1.7
#>      crea     0.514       0.018    1.4
#>       sex     0.479       0.020    1.5
#> 
#> ℹ splits_lower and splits_upper hold the 95% interval, not shown above.
```