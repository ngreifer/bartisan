
``` r
loo::loo_compare(
  loo::loo(bartisan(count ~ ., d, family = poisson(), control = ctrl)),
  loo::loo(bartisan(count ~ ., d, family = negbin(), control = ctrl))
)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#>   model elpd_diff se_diff p_worse       diag_diff       diag_elpd
#>  model1       0.0     0.0      NA                 8 k_psis > 0.54
#>  model2      -1.0     1.6    0.73 |elpd_diff| < 4 1 k_psis > 0.54
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
```