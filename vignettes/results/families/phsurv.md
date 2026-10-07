
``` r
fit_ph <- bartisan(survival::Surv(time, event) ~ x1 + x2,
                   data = d, control = ctrl,
                   family = ph())

predict(fit_ph, type = "survival", times = c(1, 2, 5)) |>
  head(3)
#>          1     2     5
#> [1,] 0.904 0.765 0.524
#> [2,] 0.896 0.745 0.493
#> [3,] 0.830 0.609 0.303
```