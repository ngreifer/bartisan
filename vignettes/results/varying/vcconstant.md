
``` r
set.seed(2026)

constant <- bartisan(death ~ age + sex + race + edu + aps +
                       meanbp + resp + hema + pafi +
                       paco2 + crea + surv2m + card +
                       vc(rhc, ~ 1),
                     data = rhc, family = binomial())

constant
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = death ~ age + sex + race + edu + aps + meanbp + 
#>     resp + hema + pafi + paco2 + crea + surv2m + card + vc(rhc, 
#>     ~1), data = rhc, family = binomial())
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1500
#> Structure: 2 forests of 50 trees, soft decision rules
#> Draws: 800 kept after 200 warmup

library(loo)
loo_compare(list(constant = loo(constant), varying = loo(fit_vc)))
#>     model elpd_diff se_diff p_worse       diag_diff diag_elpd
#>  constant       0.0     0.0      NA                          
#>   varying      -0.5     1.0    0.68 |elpd_diff| < 4
```