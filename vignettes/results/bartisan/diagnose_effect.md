
``` r
diagnose(estimate_effect(fit, treat = "rhc"))
#> Convergence and mixing
#> 
#>     quantity  rhat rhat_late ess_bulk ess_tail
#>  Y[1] - Y[0] 1.008     0.998      492      393
#>         Y[0] 1.002     1.002      496      736
#>         Y[1] 1.007     1.007      564      615
#> 
#> ✖ Only one chain, so R-hat can only compare it with itself; set `chains = 4`
#> ✔ R-hat is below 1.01 for every reported quantity
#> ✔ Warmup was long enough, since R-hat is already fine
#> ✔ Bulk ESS is at least 492 for every reported quantity, above 400
#> ✖ Tail ESS is 393 for Y[1] - Y[0], below 400
#> 
#> What to do
#> 
#> • Refit with `chains = 4`. R-hat compares chains against each other, and one
#>   chain can only be compared with itself, so nothing below is reliable until
#>   there are several. With future installed the chains run in parallel.
#> • Raise `num_draws`, which was 800. The chains agree and are stationary, so
#>   they simply have not run long enough. Do not reach for `num_thin`: thinning
#>   discards draws already paid for and lowers the effective sample size per unit
#>   of time.
#> • The tail is the binding constraint, so a posterior mean is already fine and
#>   an interval endpoint is not. Raise `num_draws`, which was 800 if intervals
#>   are what gets reported.
```