
``` r
estimate_effect(fit_bcf)
#> Average treatment effect (difference)
#> 
#> Treatment: `rhc`
#> Averaged over 1500 units
#> 
#>     contrast estimate   lower upper    n
#>  Y[1] - Y[0]   0.0563 0.00241 0.108 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.633 0.605 0.662
#>      Y[1]    0.689 0.650 0.728
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `rhc` set to "a".
```