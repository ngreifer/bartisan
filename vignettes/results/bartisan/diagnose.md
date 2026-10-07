
``` r
diagnose(fit)
#> Convergence and mixing
#> 
#>                             quantity  rhat rhat_late ess_bulk ess_tail
#>                               loglik 1.039     1.049       95      187
#>                           splits.eta 1.073     1.152       23      134
#>  eta.eta (average over observations) 1.000     1.003      435      723
#>   eta.eta (worst 5% of observations) 1.038     1.069       77      197
#>       bandwidth (average over trees) 1.022     1.037      172      299
#>        bandwidth (worst 5% of trees) 1.064     1.111       50       59
#> 
#> ✖ Only one chain, so R-hat can only compare it with itself; set `chains = 4`
#> ✖ R-hat is above 1.01 for loglik
#> ✖ That R-hat rests on only 95 effective draws, where a single chain averages
#>   1.011 even when its two halves agree
#> ℹ A longer warmup is not the fix: R-hat stays high on the second half of the
#>   draws alone as well
#> ✖ The chains disagree about how many splitting rules the forest has (R-hat
#>   1.07)
#> ✖ The chains disagree about how wide the decision rules are (R-hat 1.06)
#> ✖ Bulk ESS is 77 for eta.eta (worst 5% of observations), below 400
#> ✖ Tail ESS is 187 for loglik, below 400
#> 
#> What to do
#> 
#> • Refit with `chains = 4`. R-hat compares chains against each other, and one
#>   chain can only be compared with itself, so nothing below is reliable until
#>   there are several. With future installed the chains run in parallel.
#> • Raise `num_draws`, which was 800. R-hat is above the threshold for a quantity
#>   that carries too few effective draws for the threshold to mean anything: with
#>   this many chains it would sit about where it does even if the chains agreed
#>   exactly, as the check above reports. A larger effective sample size makes it
#>   readable, and that grows with the total number of draws; using fewer chains
#>   lowers the bar as well, since R-hat's null rises with the number of chains
#>   being compared.
#> • If that does not settle it, reduce `num_trees`, which was 50. A smaller
#>   forest has fewer ways to represent the same fit, so the sampler has less room
#>   to move between them.
#> • Then check the family. A likelihood that fits the data badly can give a
#>   posterior with no single place to be; `bayesplot::pp_check()` is the
#>   diagnostic.
#> • The forest's own size is a different kind of failure from the others and does
#>   not take the same advice. A sum of trees represents one function through many
#>   different partitions, so two chains can agree about every fitted value while
#>   disagreeing about how many rules they used to get there, and the quantities a
#>   fit reports are integrals over that structure. Raising `num_draws` moves this
#>   row slowly and may not clear the threshold at any affordable length. Act on
#>   it when split counts are themselves what gets reported --
#>   `variable_importance()` and `vignette("importance")` -- and not when fitted
#>   values, predictions or effects are.
```