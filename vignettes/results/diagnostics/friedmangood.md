
``` r
long_enough <- bartisan(y ~ ., fr, family = gaussian(), chains = 4)

diagnose(long_enough)$table
#>                              quantity rhat rhat_late ess_bulk ess_tail ess_frac
#> 1                              loglik 1.08      1.05     45.6    160.7  0.01424
#> 2                           aux.sigma 1.02      1.01    319.3   1324.8  0.09977
#> 3                          splits.eta 1.10      1.17     31.5    188.4  0.00983
#> 4 eta.eta (average over observations) 1.00      1.00   3357.2   3015.1  1.04912
#> 5  eta.eta (worst 5% of observations) 1.17      1.22     16.3     70.0  0.00511
#> 6      bandwidth (average over trees) 1.07      1.17     76.7    236.9  0.02398
#> 7       bandwidth (worst 5% of trees) 1.34      1.45     11.1     30.4  0.00348
#>   rhat_bad late_bad
#> 1    1.000    1.000
#> 2    1.000    0.000
#> 3    1.000    1.000
#> 4    0.000    0.000
#> 5    0.993    0.995
#> 6    1.000    1.000
#> 7    1.000    1.000
```