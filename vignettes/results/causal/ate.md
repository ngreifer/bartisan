
``` r
ate <- estimate_effect(fit, treat = "rhc")

ate
#> Average treatment effect (difference)
#> 
#> Treatment: `rhc`
#> Averaged over 1500 units
#> 
#>     contrast estimate  lower upper    n
#>  Y[1] - Y[0]    0.063 0.0145 0.113 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.631 0.600 0.659
#>      Y[1]    0.694 0.655 0.730
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `rhc` set to "a".
```