
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
#>  Y[1] - Y[0] 1.02      1.04      245      919
#>         Y[0] 1.01      1.01      534     1526
#>         Y[1] 1.02      1.03      244      875
#> 
#> ✔ 4 chains, 3200 draws kept in total
#> ✖ R-hat is above 1.01 for Y[1] - Y[0]
#> ✖ That R-hat rests on only 245 effective draws, where 4 chains average 1.016
#>   even when they agree
#> ℹ A longer warmup is not the fix: R-hat stays high on the second half of the
#>   draws alone as well
#> ✖ Bulk ESS is 244 for Y[1], below 400
#> ✔ Tail ESS is at least 875 for every reported quantity, above 400
#> 
#> What to do
#> 
#> • Raise `num_draws`, which was 800. R-hat is above the threshold for a quantity
#>   that carries too few effective draws for the threshold to mean anything: with
#>   this many chains it would sit about where it does even if the chains agreed
#>   exactly, as the check above reports. A larger effective sample size makes it
#>   readable, and that grows with the total number of draws; using fewer chains
#>   lowers the bar as well, since R-hat's null rises with the number of chains
#>   being compared.
#> • If that does not settle it, reduce `num_trees`. A smaller forest has fewer
#>   ways to represent the same fit, so the sampler has less room to move between
#>   them.
#> • Then check the family. A likelihood that fits the data badly can give a
#>   posterior with no single place to be; `bayesplot::pp_check()` is the
#>   diagnostic.
```