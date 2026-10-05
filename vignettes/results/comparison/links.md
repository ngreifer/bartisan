
``` r
set.seed(2026)
probit <- bartisan(model, data = rhc, family = binomial("probit"))

loo_compare(list(logit  = loo(full),
                 probit = loo(probit)))
#>   model elpd_diff se_diff p_worse       diag_diff diag_elpd
#>   logit       0.0     0.0      NA                          
#>  probit      -0.2     1.0    0.60 |elpd_diff| < 4
```