
``` r
friedman <- function(n) {
  x <- as.data.frame(matrix(runif(n * 10), n, 10))
  names(x) <- paste0("x", 1:10)
  x$y <- 10 * sin(pi * x$x1 * x$x2) + 20 * (x$x3 - 0.5)^2 +
    10 * x$x4 + 5 * x$x5 + rnorm(n)
  x
}

set.seed(7)
fit_fr <- bartisan(y ~ ., data = friedman(500), family = gaussian(),
                   sparsity = TRUE)

variable_importance(fit_fr)
#> Variable importance
#> 
#>  variable prop_used prop_splits splits
#>        x2     1.000       0.418   34.1
#>        x1     1.000       0.219   17.8
#>        x3     1.000       0.186   15.2
#>        x4     1.000       0.123   10.1
#>        x5     1.000       0.051    4.2
#>        x8     0.056       0.001    0.1
#>        x7     0.055       0.001    0.1
#>       x10     0.045       0.001    0.1
#>        x6     0.040       0.001    0.0
#>        x9     0.021       0.000    0.0
#> 
#> ℹ splits_lower and splits_upper hold the 95% interval, not shown above.
```