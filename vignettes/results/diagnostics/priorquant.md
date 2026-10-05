
``` r
p <- c(.01, .1, .5, .9, .99)

quantile(predict(prior, type = "response", draws = TRUE), p)
#>      1%     10%     50%     90%     99% 
#> 0.00534 0.20691 0.65381 0.93274 0.99801

quantile(predict(fit, type = "response", draws = TRUE), p)
#>    1%   10%   50%   90%   99% 
#> 0.183 0.365 0.679 0.895 0.952
```