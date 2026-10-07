
``` r
d$rate <- rbeta(n, plogis(1.5 * sin(pi * d$x1)) * 12,
                12 - plogis(1.5 * sin(pi * d$x1)) * 12)

fit_beta <- bartisan(rate ~ x1 + x2, data = d,
                     family = Beta(), control = ctrl)

mean(fit_beta$aux[, "phi"])
#> [1] 13.8
```