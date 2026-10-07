
``` r
estimate_effect(fit, treat = "rhc", comparison = "lnor")
#> Average treatment effect (log odds ratio)
#> 
#> Treatment: `rhc`
#> Averaged over 1500 units
#> 
#>                contrast estimate  lower upper    n
#>  log(O(Y[1]) / O(Y[0]))    0.284 0.0633 0.513 1500
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]    0.631 0.600 0.659
#>      Y[1]    0.694 0.655 0.730
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `rhc` set to "a", and O(y) is the odds
#>   `y/(1-y)`.
```