
``` r
# Binning d$heavy
nbins <- 25

d$binned <- d$heavy |>
  quantile(seq(0, 1, length.out = nbins + 1)) |>
  unique() |>
  cut(x = d$heavy, include.lowest = TRUE, labels = FALSE) |>
  ave(x = d$heavy)

fit_oc <- bartisan(binned ~ x1 + x2,
                   data = d, control = ctrl,
                   family = ordinal("probit"))

head(predict(fit_oc, type = "mean"))
#> [1] 2.147 2.470 1.351 2.252 2.746 0.884
```