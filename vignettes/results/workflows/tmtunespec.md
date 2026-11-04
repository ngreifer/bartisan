
``` r
tune_spec <- bart(trees = tune(), prior_outcome_range = tune()) |>
  set_engine("bartisan", family = "binomial", num_draws = 400) |>
  set_mode("classification")

bart_wf <- workflow() |>
  add_formula(died ~ rhc + age + sex + race + edu + aps + meanbp +
                resp + hema + pafi + paco2 + crea + surv2m + card) |>
  add_model(tune_spec)
```