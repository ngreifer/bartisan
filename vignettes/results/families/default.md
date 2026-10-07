
``` r
d$binary <- rbinom(n, 1, 0.4)

fit <- bartisan(binary ~ x1 + x2, data = d, control = ctrl)
#> ℹ Using `family = binomial()`.
#> ℹ Set `family` explicitly to silence this message.

fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = binary ~ x1 + x2, data = d, control = ctrl)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 300
#> Structure: 1 forest of 10 trees, soft decision rules
#> Draws: 150 kept after 150 warmup
```