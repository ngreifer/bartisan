
``` r
diagnose(bcf_fit)$table
#>                                       quantity rhat rhat_late ess_bulk ess_tail
#> 1                                       loglik 1.07      1.15    38.21    105.0
#> 2                                  aux.b.rhc.0 1.15      1.16    25.94    131.8
#> 3                                  aux.b.rhc.1 1.07      1.09    44.10    169.6
#> 4                           splits.(Intercept) 1.01      1.02   241.54    389.4
#> 5                                   splits.rhc 1.02      1.03   243.32    554.6
#> 6  eta.(Intercept) (average over observations) 1.35      1.25     9.26     14.8
#> 7   eta.(Intercept) (worst 5% of observations) 1.32      1.18     9.93     17.6
#> 8          eta.rhc (average over observations) 1.38      1.51    10.65     34.1
#> 9           eta.rhc (worst 5% of observations) 1.43      1.45     9.96     21.3
#> 10              bandwidth (average over trees) 1.00      1.00   528.38   1104.5
#> 11               bandwidth (worst 5% of trees) 1.02      1.05   252.97    310.3
#>    ess_frac rhat_bad late_bad
#> 1   0.01194    1.000      1.0
#> 2   0.00811    1.000      1.0
#> 3   0.01378    1.000      1.0
#> 4   0.07548    1.000      1.0
#> 5   0.07604    1.000      1.0
#> 6   0.00289    1.000      1.0
#> 7   0.00310    1.000      1.0
#> 8   0.00333    1.000      1.0
#> 9   0.00311    1.000      1.0
#> 10  0.16512    0.000      0.0
#> 11  0.07905    0.467      0.8
```