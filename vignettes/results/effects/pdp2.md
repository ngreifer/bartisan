
``` r
pd2 <- partial_dependence(fit, ~ aps + rhc)

pd2
#> Partial dependence
#> 
#> Predictors: "aps" and "rhc"
#> Averaged over 1500 units, on the "response" scale
#> 
#>     aps rhc estimate lower upper
#>    4.00   0    0.568 0.460 0.638
#>    9.72   0    0.568 0.462 0.637
#>   15.44   0    0.569 0.472 0.637
#>   21.16   0    0.573 0.493 0.637
#>   26.88   0    0.581 0.518 0.640
#>     --- 42 rows omitted ---
#>  124.12   1    0.758 0.680 0.840
#>  129.84   1    0.758 0.680 0.840
#>  135.56   1    0.759 0.680 0.841
#>  141.28   1    0.759 0.680 0.841
#>  147.00   1    0.759 0.680 0.841
#> 
#> ℹ estimate is the posterior mean; lower and upper bound the 95% equal-tailed
#>   credible interval.
#> ℹ `n_print` in `print()` (`?bartisan::print.bartisan_partial()`) sets how many
#>   rows are shown, half from each end; `print(., n_print = Inf)` shows all of
#>   them.

plot(pd2) +
  labs(x = "APACHE III score on day 1", y = "Fitted probability of death",
       color = "Catheterized", fill = "Catheterized")
```

<img src="results/effects/pdp2-1.png" alt="" style="display: block; margin: auto;" />