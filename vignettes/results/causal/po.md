
``` r
rhc$prop_score <- prop_score

fit <- bartisan(death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
                  hema + pafi + paco2 + crea + surv2m + card + prop_score,
                data = rhc, family = binomial())
```