
``` r
at <- quantile(rhc$aps, c(.1, .5, .9))

# Uses posterior mean
partial_dependence(fit, ~ aps, values = list(aps = at))
#> Partial dependence
#> 
#> Predictor: "aps"
#> Averaged over 1500 units, on the "response" scale
#> 
#>   aps estimate lower upper
#>  29.9    0.610 0.552 0.657
#>  54.0    0.656 0.623 0.695
#>  83.0    0.703 0.649 0.765
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.

# Uses posterior median by default
avg_predictions(fit, variables = list(aps = at))
#> 
#>   aps Estimate 2.5 % 97.5 %
#>  29.9    0.611 0.552  0.657
#>  54.0    0.656 0.623  0.695
#>  83.0    0.701 0.649  0.765
#> 
#> Type: response
```