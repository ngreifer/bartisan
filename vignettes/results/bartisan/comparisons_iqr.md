
``` r
marginaleffects::avg_comparisons(fit, variables = list(surv2m = "iqr"))
#> 
#>  Estimate  2.5 % 97.5 %
#>    -0.282 -0.355 -0.207
#> 
#> Term: surv2m
#> Type: response
#> Comparison: Q3 - Q1
```