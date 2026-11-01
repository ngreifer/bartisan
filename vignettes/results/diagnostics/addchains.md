
``` r
set.seed(2026)
fit1 <- update(fit, chains = 1)

fit4 <- add_chains(fit1, chains = 3)

all.equal(fit4$eta, fit$eta)
#> [1] TRUE
```