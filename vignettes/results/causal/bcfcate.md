
``` r
cate <- estimate_effect(fit_bcf, estimand = "CATE", comparison = "or")

quantile(cate$estimate, probs = c(0, .25, .5, .75, 1))
#>    0%   25%   50%   75%  100% 
#> 1.141 1.324 1.420 1.523 1.839
```