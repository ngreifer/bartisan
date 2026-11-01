
``` r
p <- c(.01, .1, .5, .9, .99)

quantile(predict(prior, type = "response", draws = TRUE), p)
#>       1%      10%      50%      90%      99% 
#> 3.02e-10 2.17e-02 6.58e-01 9.95e-01 1.00e+00

quantile(predict(fit, type = "response", draws = TRUE), p)
#>    1%   10%   50%   90%   99% 
#> 0.183 0.365 0.679 0.895 0.952
```