
``` r

pd <- partial_dependence(fit, ~ surv2m)

pd
#> Partial dependence
#> 
#> Predictor: "surv2m"
#> Averaged over 1500 units, on the "response" scale
#> 
#>  surv2m estimate lower upper
#>  0.0000    0.859 0.806 0.907
#>  0.0376    0.858 0.805 0.905
#>  0.0752    0.858 0.805 0.904
#>  0.1128    0.857 0.806 0.903
#>  0.1504    0.856 0.807 0.903
#>   --- 16 rows omitted ---
#>  0.7896    0.495 0.434 0.554
#>  0.8272    0.483 0.419 0.548
#>  0.8648    0.477 0.404 0.547
#>  0.9024    0.475 0.395 0.546
#>  0.9400    0.474 0.393 0.546
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ `n_print` in `print()` (`?bartisan::print.bartisan_partial()`) sets how many
#>   rows are shown, half from each end; `print(., n_print = Inf)` shows all of
#>   them.

plot(pd) +
  ggplot2::labs(x = "Estimated probability of surviving two months",
                y = "Fitted probability of death")
```

<img src="results/bartisan/pdp-1.png" alt="" style="display: block; margin: auto;" />