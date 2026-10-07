
``` r
set.seed(2026)

fit <- bartisan(
  death ~ rhc + age + sex + race + edu + aps + meanbp + resp + hema +
    pafi + paco2 + crea + surv2m + card,
  data = rhc, family = binomial()
)

fit
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ rhc + age + sex + race + edu + aps + 
#>     meanbp + resp + hema + pafi + paco2 + crea + surv2m + card, 
#>     data = rhc, family = binomial())
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 1 forest of 50 trees, soft decision rules
#> Draws: 800 kept after 200 warmup
```