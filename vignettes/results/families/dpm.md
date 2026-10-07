
``` r
d$heavy <- 2 * sin(pi * d$x1) + d$x2 + rt(n, df = 3)

fit_dpm <- bartisan(heavy ~ x1 + x2, data = d,
                    family = dpm(), control = ctrl)
```