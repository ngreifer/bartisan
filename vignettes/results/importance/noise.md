
``` r
set.seed(11)
rhc_noise <- rhc
for (j in 1:3) rhc_noise[[paste0("noise", j)]] <- rnorm(nrow(rhc_noise))

set.seed(2026)
fit_noise <- bartisan(update(model, . ~ . + noise1 + noise2 + noise3),
                      data = rhc_noise, family = binomial(), sparsity = TRUE)

variable_importance(fit_noise)
#> Variable importance
#> 
#>  variable prop_used prop_splits splits
#>    surv2m     1.000       0.172   13.1
#>     paco2     1.000       0.151   11.6
#>       age     1.000       0.121    9.2
#>      pafi     1.000       0.072    5.5
#>       rhc     0.980       0.045    3.4
#>       edu     0.975       0.053    4.0
#>       aps     0.894       0.069    5.3
#>    noise1     0.774       0.039    3.0
#>      crea     0.686       0.039    3.0
#>    meanbp     0.674       0.047    3.6
#>    noise3     0.670       0.032    2.5
#>      card     0.662       0.030    2.3
#>      resp     0.639       0.026    2.0
#>    noise2     0.601       0.026    2.0
#>      hema     0.598       0.028    2.1
#>       sex     0.562       0.028    2.2
#>      race     0.498       0.022    1.6
#> 
#> ℹ splits_lower and splits_upper hold the 95% interval, not shown above.
```