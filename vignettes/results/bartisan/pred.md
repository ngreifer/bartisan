
``` r
new_patient <- rhc[1, ]
new_patient$rhc <- 1

predict(fit, newdata = new_patient)
#> [1] 0.8366
```