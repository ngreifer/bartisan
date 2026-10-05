
``` r
loo_compare(list(bart     = loo(full),
                 logistic = loo(logistic)))
#>     model elpd_diff se_diff p_worse       diag_diff diag_elpd
#>  logistic       0.0     0.0      NA                          
#>      bart      -2.8     5.0    0.71 |elpd_diff| < 4
```