
``` r
eff_180 <- avg_comparisons(fit_ph, variables = "rhc", type = "survival",
                           times = 180)

eff_180
#> 
#>  Estimate  2.5 %  97.5 %
#>   -0.0617 -0.107 -0.0155
#> 
#> Term: rhc
#> Type: survival
#> Comparison: 1 - 0
```