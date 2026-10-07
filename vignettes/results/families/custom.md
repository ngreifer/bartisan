
``` r
# Custom Poisson family
pois_by_hand <- custom_family(
  logdens = function(y, eta) dpois(y, exp(eta[, 1]), log = TRUE),
  start   = log(mean(d$count)),
  name    = "hand-rolled Poisson")

fit_custom <- bartisan(count ~ x1 + x2,
                       data = d, control = ctrl,
                       family = pois_by_hand)

# Internal Poisson family
fit_pois <- bartisan(count ~ x1 + x2,
                     data = d, control = ctrl,
                     family = poisson())

cor(predict(fit_custom, type = "link"),
    predict(fit_pois, type = "link"))
#> [1] 0.997
```