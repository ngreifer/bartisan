
``` r
# Custom Gaussian family
gauss_by_hand <- custom_family(
  logdens   = function(y, eta, aux) dnorm(y, eta[, 1], exp(aux[1]), log = TRUE),
  aux_start = c(log_sigma = 0),
  start     = mean(d$heavy))

fit_aux <- bartisan(heavy ~ x1 + x2,
                    data = d, control = ctrl,
                    family = gauss_by_hand)

# Internal Gaussian family
fit_gauss <- bartisan(heavy ~ x1 + x2,
                      data = d, control = ctrl,
                      family = gaussian())

cor(predict(fit_aux, type = "link"),
    predict(fit_gauss, type = "link"))
#> [1] 0.991

# Estimate of the auxiliary parameter
summary(exp(fit_aux$aux[, "log_sigma"]))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    1.55    1.69    1.75    1.75    1.80    1.93
summary(fit_gauss$aux[, "sigma"])
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    1.59    1.70    1.75    1.75    1.80    2.00
```