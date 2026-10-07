
``` r
set.seed(2026)

n_groups <- 40L
d <- data.frame(g = factor(sample(n_groups, 600L, replace = TRUE)),
                x = runif(600L))
u <- rnorm(n_groups, 0, 0.5)
d$y <- sin(2 * pi * d$x) + u[d$g] + rnorm(600L, 0, 0.5)

fit_re <- bartisan(y ~ x + (1 | g), data = d)

fit_re
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = y ~ x + (1 | g), data = d)
#> 
#> Family: "dpm" with the "identity" link
#> Observations: 600
#> Structure: 1 forest of 50 trees, soft decision rules
#> Draws: 800 kept after 200 warmup
#> Random intercepts: g (40 levels)
#> 
#> Posterior means: alpha = 3.37, clusters = 15.7, center = -0.0182, error_sd = 0.521
```