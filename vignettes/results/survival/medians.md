
``` r
med <- predict(fit_ph, type = "response")

summary(med)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     9.2    55.9   206.9   240.7   390.0   860.0
```