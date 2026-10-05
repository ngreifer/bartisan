
``` r
library(ggplot2)

pd <- partial_dependence(fit, ~aps)

pd
#> Partial dependence
#> 
#> Predictor: "aps"
#> Averaged over 1500 units, on the "response" scale
#> 
#>     aps estimate lower upper
#>    4.00    0.592 0.491 0.657
#>    9.72    0.593 0.495 0.657
#>   15.44    0.593 0.499 0.657
#>   21.16    0.597 0.513 0.656
#>   26.88    0.606 0.543 0.658
#>   --- 16 rows omitted ---
#>  124.12    0.724 0.651 0.811
#>  129.84    0.724 0.651 0.811
#>  135.56    0.724 0.651 0.811
#>  141.28    0.724 0.651 0.812
#>  147.00    0.724 0.651 0.812
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ `n_print` in `print()` (`?bartisan::print.bartisan_partial()`) sets how many
#>   rows are shown, half from each end; `print(., n_print = Inf)` shows all of
#>   them.

plot(pd) +
  labs(x = "APACHE III score on day 1", y = "Fitted probability of death")
```

<img src="results/effects/pdp-1.png" alt="" style="display: block; margin: auto;" />

``` r

# Same thing:
## plot(fit, ~aps)
```