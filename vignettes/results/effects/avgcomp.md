
``` r
comp <- avg_comparisons(fit, variables = c("rhc", "age", "card"))

comp
#> 
#>  Term Contrast Estimate    2.5 %  97.5 %
#>  age  +1        0.00309  0.00153 0.00467
#>  card yes - no  0.04164 -0.00868 0.08973
#>  rhc  1 - 0     0.05968  0.01182 0.11118
#> 
#> Type: response
```