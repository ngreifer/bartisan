
``` r
set.seed(2026)

fit_vc <- bartisan(death ~ age + sex + race + edu + aps +
                     meanbp + resp + hema + pafi +
                     paco2 + crea + surv2m + card +
                     vc(rhc),
                   data = rhc, family = binomial())

fit_vc
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ age + sex + race + edu + aps + meanbp + 
#>     resp + hema + pafi + paco2 + crea + surv2m + card + vc(rhc), 
#>     data = rhc, family = binomial())
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 2 forests of 50 trees, soft decision rules
#> Draws: 800 kept after 200 warmup
```