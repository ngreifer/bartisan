
``` r
long_enough <- bartisan(y ~ ., fr, family = gaussian(), chains = 4)

diagnose(long_enough)$table
#>                              quantity rhat rhat_late ess_bulk ess_tail ess_frac
#> 1                              loglik 1.09      1.14     29.4    106.2  0.00919
#> 2                           aux.sigma 1.02      1.03    140.3   1099.2  0.04384
#> 3                          splits.eta 1.14      1.16     21.1    175.4  0.00661
#> 4 eta.eta (average over observations) 1.00      1.00   2766.9   3137.1  0.86465
#> 5  eta.eta (worst 5% of observations) 1.16      1.19     17.2     77.5  0.00536
#> 6      bandwidth (average over trees) 1.14      1.30     21.9     67.0  0.00684
#> 7       bandwidth (worst 5% of trees) 1.41      1.64     10.4     25.0  0.00324
#>   rhat_bad late_bad
#> 1    1.000    1.000
#> 2    1.000    1.000
#> 3    1.000    1.000
#> 4    0.000    0.000
#> 5    0.993    0.998
#> 6    1.000    1.000
#> 7    1.000    1.000
```