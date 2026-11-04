
``` r
set.seed(1234)
rhc_folds <- vfold_cv(rhc_train, v = 5, strata = died)

bart_grid <- expand.grid(trees = c(20, 50),
                         prior_outcome_range = c(1, 2, 4))

bart_res <- tune_grid(bart_wf, resamples = rhc_folds,
                      grid = bart_grid,
                      metrics = metric_set(roc_auc, brier_class))

show_best(bart_res, metric = "roc_auc")
#> # A tibble: 5 × 8
#>   trees prior_outcome_range .metric .estimator  mean     n std_err .config      
#>   <dbl>               <dbl> <chr>   <chr>      <dbl> <int>   <dbl> <chr>        
#> 1    50                   1 roc_auc binary     0.737     5  0.0211 pre0_mod4_po…
#> 2    50                   4 roc_auc binary     0.737     5  0.0206 pre0_mod6_po…
#> 3    20                   4 roc_auc binary     0.737     5  0.0216 pre0_mod3_po…
#> 4    50                   2 roc_auc binary     0.737     5  0.0218 pre0_mod5_po…
#> 5    20                   2 roc_auc binary     0.736     5  0.0219 pre0_mod2_po…
```