
``` r
kfold_demographics <- kfold(demographics, folds = folds)

loo_compare(list(full         = kfold_full,
                 demographics = kfold_demographics))
#>         model elpd_diff se_diff p_worse diag_diff diag_elpd
#>          full       0.0     0.0      NA                    
#>  demographics     -72.3    10.9    1.00
```