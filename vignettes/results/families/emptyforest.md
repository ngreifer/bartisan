
``` r
set.seed(8)
dz2 <- data.frame(x1 = runif(500), x2 = runif(500))
dz2$y <- rnorm(500, 2 * dz2$x1, 0.7)

flat <- bartisan(list(y ~ x1 + x2, ~ 1), data = dz2,
                 family = gaussian_ls(), control = ctrl)

# no splitting rules in the scale forest, and one value of sigma
c(scale_splits = sum(flat$counts$log_sd),
  sigma = mean(exp(flat$eta[[2L]])))
#> scale_splits        sigma 
#>        0.000        0.718
```