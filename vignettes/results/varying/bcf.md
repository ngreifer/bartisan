
``` r
set.seed(2026)

fit_bcf <- bcf(death ~ age + sex + race + edu + aps +
                 meanbp + resp + hema + pafi +
                 paco2 + crea + surv2m + card,
               treat = ~ rhc, data = rhc, family = binomial())

fit_bcf
#> Generalized BART
#> 
#> Call:
#> bcf(formula = death ~ age + sex + race + edu + aps + meanbp + 
#>     resp + hema + pafi + paco2 + crea + surv2m + card, treat = ~rhc, 
#>     data = rhc, family = binomial())
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 2 forests of 50 and 25 trees, soft decision rules
#> Draws: 800 kept after 200 warmup
#> 
#> Posterior means: b.rhc.0 = -0.0312, b.rhc.1 = 0.0635
#> 
#> Treatment: "rhc"
#> Effect moderators: "age", "sex", "race", "edu", "aps", "meanbp", "resp", "hema", "pafi", "paco2", "crea", "surv2m", and "card"
#> ℹ `estimate_effect()` reports the treatment effect, with the average potential
#>   outcomes beside it; `plot()` draws the conditional ones.
```