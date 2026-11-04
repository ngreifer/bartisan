
``` r
predict(sl, newdata = rhc[1:5, covs], onlySL = TRUE)$pred
#>       [,1]
#> [1,] 0.767
#> [2,] 0.752
#> [3,] 0.199
#> [4,] 0.390
#> [5,] 0.363
```