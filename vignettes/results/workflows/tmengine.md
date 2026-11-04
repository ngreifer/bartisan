
``` r
extract_fit_engine(bart_fit)
#> Generalized BART
#> 
#> Call:
#> bartisan::bartisan(formula = died ~ rhc + age + sex + race + 
#>     edu + aps + meanbp + resp + hema + pafi + paco2 + crea + 
#>     surv2m + card, data = data, family = "binomial", num_trees = 50)
#> 
#> Family: "binomial" with the "logit" link
#> Observations: 1124
#> Structure: 1 forest of 50 trees, soft decision rules
#> Draws: 800 kept after 200 warmup
```