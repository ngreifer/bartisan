
``` r
set.seed(1234)
bart_fit <- bart_spec |>
  fit(died ~ rhc + age + sex + race + edu + aps + meanbp + resp +
        hema + pafi + paco2 + crea + surv2m + card,
      data = rhc_train)
```