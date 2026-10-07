
``` r
b <- coef(constant, draws = TRUE)[["rhc"]][, 1]

c(log_OR = mean(b),
  OR = exp(mean(b)),
  lower = quantile(exp(b), .025, names = FALSE),
  upper = quantile(exp(b), .975, names = FALSE)) |>
  round(3)
#> log_OR     OR  lower  upper 
#>  0.338  1.402  1.059  1.799
```