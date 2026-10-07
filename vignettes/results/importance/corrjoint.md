
``` r
lo <- transform(cd, x1 = 0.25, x1_copy = 0.25)
hi <- transform(cd, x1 = 0.75, x1_copy = 0.75)

drawn <- rowMeans(predict(fit_corr, newdata = hi, draws = TRUE) -
                    predict(fit_corr, newdata = lo, draws = TRUE))

round(c(estimate = mean(drawn), quantile(drawn, c(0.025, 0.975))), 3)
#> estimate     2.5%    97.5% 
#>     1.59     1.46     1.73
```