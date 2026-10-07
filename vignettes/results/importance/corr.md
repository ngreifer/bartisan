
``` r
set.seed(5)
n <- 400
cd <- data.frame(x1 = runif(n))
cd$x1_copy <- cd$x1 + rnorm(n, sd = 0.01)   # almost the same variable
cd$x2 <- runif(n)
cd$x3 <- runif(n)
cd$y <- 3 * cd$x1 + rnorm(n, sd = 0.3)      # only x1 is in the truth

fit_corr <- bartisan(y ~ ., data = cd, family = gaussian(), sparsity = TRUE)

variable_importance(fit_corr)
#> Variable importance
#> 
#>  variable prop_used prop_splits splits
#>   x1_copy     1.000       0.660   50.3
#>        x1     0.797       0.317   24.3
#>        x3     0.445       0.017    1.3
#>        x2     0.152       0.006    0.4
#> 
#> ℹ splits_lower and splits_upper hold the 95% interval, not shown above.
```