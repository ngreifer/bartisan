
``` r
library(bartisan)
library(loo)

data("rhc")

model <- death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
  hema + pafi + paco2 + crea + surv2m + card

set.seed(2026)

full <- bartisan(model, data = rhc, family = binomial())
```