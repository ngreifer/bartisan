
``` r
d$time <- rexp(n, exp(-(1 + d$x1)))
d$event <- rbinom(n, 1, 0.7)

fit_aft <- bartisan(survival::Surv(time, event) ~ x1 + x2,
                    data = d, control = ctrl,
                    family = weibull_aft())

colMeans(fit_aft$aux)
#> sigma 
#> 0.979
```