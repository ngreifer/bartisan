
``` r
counts <- variable_importance(fit, draws = TRUE)

mean(counts[, "aps"] > counts[, "meanbp"])
#> [1] 0.649
```