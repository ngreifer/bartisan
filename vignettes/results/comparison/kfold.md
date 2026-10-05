
``` r
set.seed(2026)

# K = 5 here; more is better but slower
folds <- kfold_split_random(K = 5, N = nrow(rhc))

# Refit the models K times
kfold_full <- kfold(full, folds = folds)

kfold_full
#> 
#> Based on 5-fold cross-validation.
#> 
#>            Estimate   SE
#> elpd_kfold   -843.5 17.0
#> p_kfold        31.5  2.0
#> kfoldic      1686.9 33.9
```