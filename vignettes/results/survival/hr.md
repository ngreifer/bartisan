
``` r
hr <- avg_comparisons(fit_ph, variables = "rhc", type = "link",
                      transform = exp)

hr
#> 
#>  Estimate 2.5 % 97.5 %
#>      1.22  1.05   1.41
#> 
#> Term: rhc
#> Type: link
#> Comparison: 1 - 0
```