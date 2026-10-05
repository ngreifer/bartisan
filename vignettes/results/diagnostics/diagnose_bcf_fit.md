
``` r
diagnose(bcf_fit)$table
#>                                       quantity rhat rhat_late ess_bulk ess_tail
#> 1                                       loglik 1.06      1.06     51.6    360.7
#> 2                                  aux.b.rhc.0 1.07      1.31     48.0     82.5
#> 3                                  aux.b.rhc.1 1.14      1.31     21.4     73.6
#> 4                           splits.(Intercept) 1.03      1.02    223.0    537.7
#> 5                                   splits.rhc 1.01      1.01    315.4    655.8
#> 6  eta.(Intercept) (average over observations) 1.23      1.61     12.8     31.8
#> 7   eta.(Intercept) (worst 5% of observations) 1.16      1.41     18.0     36.6
#> 8          eta.rhc (average over observations) 1.20      1.40     14.1     40.3
#> 9           eta.rhc (worst 5% of observations) 1.20      1.37     15.3     30.1
#> 10              bandwidth (average over trees) 1.01      1.01    656.7   1161.7
#> 11               bandwidth (worst 5% of trees) 1.02      1.04    274.3    312.9
#>    ess_frac rhat_bad late_bad
#> 1   0.01611    1.000    1.000
#> 2   0.01499    1.000    1.000
#> 3   0.00669    1.000    1.000
#> 4   0.06969    1.000    1.000
#> 5   0.09857    0.000    1.000
#> 6   0.00401    1.000    1.000
#> 7   0.00563    0.997    1.000
#> 8   0.00440    1.000    1.000
#> 9   0.00477    1.000    1.000
#> 10  0.20521    0.000    0.000
#> 11  0.08573    0.307    0.747
```