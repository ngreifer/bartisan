
``` r
set.seed(2026)

prior <- bartisan(death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
                    hema + pafi + paco2 + crea + surv2m + card,
                  data = rhc, family = binomial(), prior_only = TRUE)
```