
``` r
att <- estimate_effect(fit_earn_bcf, estimand = "ATT")

cate_att <- estimate_effect(fit_earn_bcf, estimand = "CATE",
                            newdata = subset(lalonde, treat == 1))

c(ATT = att$estimate, mean_CATE = mean(cate_att$estimate))
#>       ATT mean_CATE 
#>      1125      1125
```