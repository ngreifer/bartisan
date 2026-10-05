
``` r
bcf_fit <- bcf(death ~ age + sex + race + edu + aps + meanbp + resp +
                 hema + pafi + paco2 + crea + surv2m + card,
               treat = ~rhc, data = rhc,
               family = binomial("probit"), chains = 4)

estimate_effect(bcf_fit, estimand = "ATE") |>
  diagnose()
#> Convergence and mixing
#> 
#>     quantity rhat rhat_late ess_bulk ess_tail
#>  Y[1] - Y[0] 1.01      1.01      410     1041
#>         Y[0] 1.00      1.00      754     1661
#>         Y[1] 1.00      1.02      367      976
#> 
#> ✔ 4 chains, 3200 draws kept in total
#> ✔ R-hat is below 1.01 for every reported quantity
#> ✔ Warmup was long enough, since R-hat is already fine
#> ✖ Bulk ESS is 367 for Y[1], below 400
#> ✔ Tail ESS is at least 976 for every reported quantity, above 400
#> 
#> What to do
#> 
#> • Raise `num_draws`, which was 800. The chains agree and are stationary, so
#>   they simply have not run long enough. Do not reach for `num_thin`: thinning
#>   discards draws already paid for and lowers the effective sample size per unit
#>   of time.
```