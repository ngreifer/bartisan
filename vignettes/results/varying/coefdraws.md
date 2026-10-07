
``` r
draws <- coef(fit_vc, draws = TRUE)[["rhc"]]

# A 95% interval on the coefficient for the first three patients
t(apply(draws[, 1:3], 2L, quantile, c(.025, .975)))
#>           2.5%  97.5%
#> [1,] -0.002961 0.9387
#> [2,] -0.240273 0.6338
#> [3,] -0.214174 0.6488
```