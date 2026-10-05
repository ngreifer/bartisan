
``` r
loo(full)
#> 
#> Computed from 800 by 1500 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -847.6 17.3
#> p_loo        35.7  1.0
#> looic      1695.2 34.6
#> ------
#> MCSE of elpd_loo is 0.5.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.0, 0.5]).
#> 
#> All Pareto k estimates are good (k < 0.66).
#> See help('pareto-k-diagnostic') for details.
```