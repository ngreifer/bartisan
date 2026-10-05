
``` r
set.seed(3)
fr <- as.data.frame(matrix(runif(400 * 10), 400, 10))
names(fr) <- paste0("x", 1:10)
fr$y <- 10 * sin(pi * fr$x1 * fr$x2) + 20 * (fr$x3 - 0.5)^2 +
  10 * fr$x4 + 5 * fr$x5 + rnorm(400)

too_short <- bartisan(y ~ ., fr, family = gaussian(), chains = 4,
                      control = bartisan_control(num_trees = 20, num_burn = 50,
                                                 num_draws = 50))

diagnose(too_short)$table
#>                              quantity rhat rhat_late ess_bulk ess_tail ess_frac
#> 1                              loglik 3.11      3.32     5.59     14.4   0.0279
#> 2                           aux.sigma 2.56      2.30     5.92     25.4   0.0296
#> 3                          splits.eta 1.64      1.72     8.12     32.2   0.0406
#> 4 eta.eta (average over observations) 1.00      1.04   223.53    236.8   1.1177
#> 5  eta.eta (worst 5% of observations) 1.96      2.11     6.83     17.4   0.0342
#> 6      bandwidth (average over trees) 1.15      1.48    19.84    108.4   0.0992
#> 7       bandwidth (worst 5% of trees) 3.19      3.79     5.69     11.4   0.0284
#>   rhat_bad late_bad
#> 1        1        1
#> 2        1        1
#> 3        1        1
#> 4        0        1
#> 5        1        1
#> 6        1        1
#> 7        1        1
```