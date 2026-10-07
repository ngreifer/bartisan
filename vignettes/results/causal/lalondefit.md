
``` r
fit_earn_bcf <- bcf(
  re78 ~ age + educ + race + married + nodegree + re74 + re75,
  treat = ~ treat,
  data = lalonde, family = tweedie()
)

estimate_effect(fit_earn_bcf, estimand = "ATT")
#> Average treatment effect on the treated (difference)
#> 
#> Treatment: `treat`
#> Averaged over the 185 units in group "1"
#> 
#>     contrast estimate lower upper   n
#>  Y[1] - Y[0]     1120  -363  2540 185
#> 
#> Average potential outcomes
#> 
#>  quantity estimate lower upper
#>      Y[0]     5230  4240  6330
#>      Y[1]     6350  5410  7410
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ Y[a] is the average response with `treat` set to "a".
```