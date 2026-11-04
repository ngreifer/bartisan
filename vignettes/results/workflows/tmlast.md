
``` r
bart_final <- bart_wf |>
  finalize_workflow(select_best(bart_res, metric = "roc_auc")) |>
  last_fit(rhc_split, metrics = metric_set(roc_auc, brier_class))

glm_final <- workflow() |>
  add_formula(died ~ rhc + age + sex + race + edu + aps + meanbp +
                resp + hema + pafi + paco2 + crea + surv2m + card) |>
  add_model(logistic_reg()) |>
  last_fit(rhc_split, metrics = metric_set(roc_auc, brier_class))

bind_rows(BART = collect_metrics(bart_final),
          GLM = collect_metrics(glm_final),
          .id = "model")
#> # A tibble: 4 × 5
#>   model .metric     .estimator .estimate .config        
#>   <chr> <chr>       <chr>          <dbl> <chr>          
#> 1 BART  roc_auc     binary         0.736 pre0_mod0_post0
#> 2 BART  brier_class binary         0.190 pre0_mod0_post0
#> 3 GLM   roc_auc     binary         0.745 pre0_mod0_post0
#> 4 GLM   brier_class binary         0.187 pre0_mod0_post0
```