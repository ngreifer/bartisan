
``` r
library(marginaleffects)

avg_comparisons(fit_corr, variables = c("x1", "x1_copy"))
#> 
#>     Term Estimate   2.5 % 97.5 %
#>  x1         0.382 -0.0409   1.28
#>  x1_copy    0.998  0.1900   1.53
#> 
#> Type: response
#> Comparison: +1
```