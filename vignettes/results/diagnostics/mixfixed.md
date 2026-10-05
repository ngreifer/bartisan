
``` r
set.seed(2026)
longer <- update(short, num_burn = 2000, num_draws = 8000)

diagnose(longer)
#> Convergence and mixing
#> 
#>                             quantity rhat rhat_late ess_bulk ess_tail
#>                               loglik 1.02      1.03      208      977
#>                           splits.eta 1.00      1.00     2710     5914
#>  eta.eta (average over observations) 1.00      1.00    17037    23899
#>   eta.eta (worst 5% of observations) 1.01      1.01     1759     3630
#>       bandwidth (average over trees) 1.00      1.00     5219    10097
#>        bandwidth (worst 5% of trees) 1.00      1.00     3972     4697
#> 
#> ✔ 4 chains, 32000 draws kept in total
#> ✖ R-hat is above 1.01 for loglik
#> ✖ That R-hat rests on only 208 effective draws, where 4 chains average 1.019
#>   even when they agree
#> ℹ A longer warmup is not the fix: R-hat stays high on the second half of the
#>   draws alone as well
#> ✔ The chains agree about the size of the forest
#> ✔ The chains agree about how wide the decision rules are
#> ✖ Bulk ESS is 208 for loglik, below 400
#> ✔ Tail ESS is at least 977 for every reported quantity, above 400
#> ℹ Per-draw efficiency is lowest for loglik, which carries 0.6 effective draws
#>   per hundred kept
#> 
#> What to do
#> 
#> • Raise `num_draws`, which was 8000. R-hat is above the threshold for a
#>   quantity that carries too few effective draws for the threshold to mean
#>   anything: with this many chains it would sit about where it does even if the
#>   chains agreed exactly, as the check above reports. A larger effective sample
#>   size makes it readable, and that grows with the total number of draws; using
#>   fewer chains lowers the bar as well, since R-hat's null rises with the number
#>   of chains being compared.
#> • If that does not settle it, reduce `num_trees`, which was 50. A smaller
#>   forest has fewer ways to represent the same fit, so the sampler has less room
#>   to move between them.
#> • Then check the family. A likelihood that fits the data badly can give a
#>   posterior with no single place to be; `bayesplot::pp_check()` is the
#>   diagnostic.
```