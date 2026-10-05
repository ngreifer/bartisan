
``` r
avg_comparisons(fit, variables = "rhc", by = "card",
                hypothesis = ~pairwise)
#> 
#>    Hypothesis Estimate   2.5 % 97.5 %
#>  (yes) - (no) -0.00259 -0.0611 0.0248
#> 
#> Type: response
```