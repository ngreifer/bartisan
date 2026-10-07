
``` r
po_180 <- avg_predictions(fit_ph, variables = "rhc", type = "survival",
                          times = 180)

po_180
#> 
#>  rhc Estimate 2.5 % 97.5 %
#>    0    0.507 0.478  0.537
#>    1    0.445 0.410  0.481
#> 
#> Type: survival
```