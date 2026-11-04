
``` r
surv_model <- surv ~ rhc + age + sex + race + edu + aps + meanbp +
  resp + hema + pafi + paco2 + crea + surv2m + card

set.seed(1234)
dpm_fit <- bart(mode = "censored regression") |>
  set_engine("bartisan", family = "dpm_aft") |>
  fit(surv_model, data = rhc_train)

weibull_fit <- bart(mode = "censored regression") |>
  set_engine("bartisan", family = "weibull_aft") |>
  fit(surv_model, data = rhc_train)
```