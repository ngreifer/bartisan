
``` r
surv_metrics <- function(fit) {
  aug <- augment(fit, rhc_test, eval_time = c(30, 90, 180))

  bind_rows(brier_survival(aug, truth = surv, .pred),
            concordance_survival(aug, truth = surv,
                                 estimate = .pred_time))
}

surv_eval <- bind_rows(dpm_aft = surv_metrics(dpm_fit),
                       weibull_aft = surv_metrics(weibull_fit),
                       .id = "family")

surv_eval
#> # A tibble: 8 × 5
#>   family      .metric              .estimator .eval_time .estimate
#>   <chr>       <chr>                <chr>           <dbl>     <dbl>
#> 1 dpm_aft     brier_survival       standard           30     0.204
#> 2 dpm_aft     brier_survival       standard           90     0.228
#> 3 dpm_aft     brier_survival       standard          180     0.228
#> 4 dpm_aft     concordance_survival standard           NA     0.667
#> 5 weibull_aft brier_survival       standard           30     0.198
#> 6 weibull_aft brier_survival       standard           90     0.205
#> 7 weibull_aft brier_survival       standard          180     0.217
#> 8 weibull_aft concordance_survival standard           NA     0.667
```