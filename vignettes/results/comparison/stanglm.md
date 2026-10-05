
``` r
# Fit a Bayesian logistic GLM
logistic <- rstanarm::stan_glm(model, data = rhc, family = binomial(),
                               chains = 4, refresh = 0, seed = 2026)

loo(logistic)
#> 
#> Computed from 4000 by 1500 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo   -844.9 18.4
#> p_loo        16.6  0.6
#> looic      1689.7 36.8
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume independent draws (r_eff=1).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
```