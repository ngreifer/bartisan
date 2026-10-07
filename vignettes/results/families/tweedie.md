
``` r
mu <- exp(1.5 + sin(pi * d$x1))
lambda <- mu^(2 - 1.5) / (4 * (2 - 1.5))
claims <- rpois(n, lambda)
d$spend <- ifelse(claims > 0,
                  rgamma(n, shape = claims * (2 - 1.5) / (1.5 - 1),
                         scale = 4 * (1.5 - 1) * mu^(1.5 - 1)),
                  0)

fit_tw <- bartisan(spend ~ x1 + x2, data = d,
                   family = tweedie(), control = ctrl)

c(zeros = mean(d$spend == 0), phi = mean(fit_tw$aux[, "phi"]))
#> zeros   phi 
#> 0.217 3.863
```