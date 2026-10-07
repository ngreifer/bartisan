
``` r
options(marginaleffects_posterior_center = mean)

avg_comparisons(fit, variables = "rhc")
#> 
#>  Estimate  2.5 % 97.5 %
#>     0.063 0.0145  0.113
#> 
#> Term: rhc
#> Type: response
#> Comparison: 1 - 0
```