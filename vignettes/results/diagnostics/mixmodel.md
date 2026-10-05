
``` r
model <- death ~ rhc + age + aps + meanbp + surv2m

set.seed(2026)
short <- bartisan(model, data = rhc, family = binomial(), chains = 4)

diagnose(short)
#> Convergence and mixing
#> 
#>                             quantity rhat rhat_late ess_bulk ess_tail
#>                               loglik 1.04      1.11       85      277
#>                           splits.eta 1.00      1.02      349      444
#>  eta.eta (average over observations) 1.00      1.00     1809     2614
#>   eta.eta (worst 5% of observations) 1.02      1.03      343      676
#>       bandwidth (average over trees) 1.01      1.01      481     1126
#>        bandwidth (worst 5% of trees) 1.02      1.03      338      349
#> 
#> ✔ 4 chains, 3200 draws kept in total
#> ✖ R-hat is above 1.01 for loglik
#> ✖ That R-hat rests on only 85 effective draws, where 4 chains average 1.047
#>   even when they agree
#> ℹ A longer warmup is not the fix: R-hat stays high on the second half of the
#>   draws alone as well
#> ✔ The chains agree about the size of the forest
#> ✖ The chains disagree about how wide the decision rules are (R-hat 1.02)
#> ✖ Bulk ESS is 85 for loglik, below 400
#> ✖ Tail ESS is 277 for loglik, below 400
#> ℹ The chains disagree about individual observations and agree about their
#>   average (R-hat 1.01, 481 effective draws)
#> ℹ Per-draw efficiency is lowest for loglik, which carries 2.6 effective draws
#>   per hundred kept
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
#> • If that does not settle it, reduce `num_trees`, which was 50. A smaller
#>   forest has fewer ways to represent the same fit, so the sampler has less room
#>   to move between them.
#> • Then check the family. A likelihood that fits the data badly can give a
#>   posterior with no single place to be; `bayesplot::pp_check()` is the
#>   diagnostic.
#> • Note that the chains disagree about the fitted values of individual
#>   observations and not about their average, which is the usual shape of this in
#>   a forest. What that means for an estimand cannot be read off this table
#>   either way, since an estimand is a contrast and a contrast can mix badly
#>   where the function it contrasts mixes well. Compute it: `diagnose()` takes
#>   the output of `estimate_effect()`, and `posterior::as_draws()` hands the
#>   draws to `posterior::summarise_draws()` for anything else.
```