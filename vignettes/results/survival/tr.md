
``` r
tr <- avg_comparisons(fit_dpm, variables = "rhc", type = "link",
                      transform = exp)

tr
#> 
#>  Estimate 2.5 % 97.5 %
#>      0.82 0.676  0.999
#> 
#> Term: rhc
#> Type: link
#> Comparison: 1 - 0
```