
``` r
d$het <- rnorm(n, 2 * d$x1, exp(-1 + 1.5 * d$x2))

fit_sub <- bartisan(list(mean   = het ~ x1 + x2,
                         log_sd =     ~ x2),
                    data = d, control = ctrl,
                    family = gaussian_ls(),
                    num_trees = c(mean = 10, log_sd = 5),
                    sparsity = c(mean = TRUE, log_sd = FALSE))

variable_importance(fit_sub)
#> Variable importance
#> 
#>  predictor variable prop_used prop_splits splits
#>       mean       x1     1.000       0.997   15.1
#>     log_sd       x2     1.000       1.000    7.2
#>       mean       x2     0.053       0.003    0.1
#>     log_sd       x1     0.000       0.000    0.0
#> 
#> ℹ Fitted with `sparsity = FALSE` (the default), so every predictor keeps a
#>   share of the rules and prop_used is near 1 throughout. Refit with `sparsity =
#>   TRUE` to read it as a selection rule.
#> ℹ splits_lower and splits_upper hold the 95% interval, not shown above.
```