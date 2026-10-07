
``` r
hr_unit <- comparisons(fit_ph, variables = "rhc", type = "link",
                       transform = exp)

summary(hr_unit$estimate)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    1.10    1.20    1.22    1.22    1.25    1.31
```