
``` r
set.seed(2026)

fit_aps <- bartisan(death ~ age + sex + race + edu + aps +
                      meanbp + resp + hema + pafi +
                      paco2 + crea + surv2m + card +
                      vc(rhc, ~ aps + age),
                    data = rhc, family = binomial())
```