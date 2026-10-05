
``` r
trees <- c(20, 50, 200)

tuned <- lapply(trees, function(n) {
  set.seed(2026)
  bartisan(model, data = rhc, family = binomial(),
           num_trees = n)
})

names(tuned) <- paste0("trees_", trees)

# LOOCV
loo_compare(lapply(tuned, loo))
#>      model elpd_diff se_diff p_worse       diag_diff diag_elpd
#>  trees_200       0.0     0.0      NA                          
#>   trees_50       0.0     1.2    0.51 |elpd_diff| < 4          
#>   trees_20      -2.2     2.0    0.87 |elpd_diff| < 4

# K-fold CV
# loo_compare(lapply(tuned, kfold, folds = folds))
```