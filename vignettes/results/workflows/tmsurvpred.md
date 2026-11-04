
``` r
predict(dpm_fit, rhc_test[1:3, ], type = "time")
#> # A tibble: 3 × 1
#>   .pred_time
#>        <dbl>
#> 1       121.
#> 2       293.
#> 3       228.

surv_pred <- predict(dpm_fit, rhc_test[1:3, ], type = "survival",
                     eval_time = c(30, 180),
                     interval = "confidence")

surv_pred$.pred[[1]]
#> # A tibble: 2 × 4
#>   .eval_time .pred_survival .pred_lower .pred_upper
#>        <dbl>          <dbl>       <dbl>       <dbl>
#> 1         30          0.644       0.602       0.694
#> 2        180          0.422       0.328       0.496
```