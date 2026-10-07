
``` r
summary(fit)
#> Convergence and mixing
#> 
#> Log likelihood: R-hat 1.039, bulk ESS 95, tail ESS 187, over 1 chain
#> 
#> ℹ Use diagnose() (`?bartisan::diagnose`) to examine convergence and mixing
#>   diagnostics.
#> 
#> Variable importance
#> 
#> Predictors ranked by use; 5 of 14 shown.
#>  variable prop_used prop_splits splits
#>    surv2m         1       0.116    8.5
#>       age         1       0.089    6.5
#>     paco2         1       0.073    5.4
#>       rhc         1       0.070    5.1
#>      pafi         1       0.069    5.1
#> 
#> ℹ Use variable_importance() (`?bartisan::variable_importance`) to examine
#>   variable importance.
#> 
#> Further tools
#> 
#> ℹ Use loo() (`?bartisan::loo.bartisan_fit`) to compare this fit with others, or
#>   kfold() (`?bartisan::kfold.bartisan_fit`) if `loo()` reports many Pareto k
#>   values above 0.7.
#> ℹ Use partial_dependence() (`?bartisan::partial_dependence`) and plot()
#>   (`?bartisan::plot.bartisan_fit`) to view the partial dependence of the
#>   predictions on a predictor.
```