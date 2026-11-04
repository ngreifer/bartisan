
``` r
rhc_mis <- rhc[c("death", "age", "sex", "aps", "meanbp", "pafi",
                 "crea", "card")]

set.seed(1234)
p_mis <- plogis(-2 + .03 * (rhc$age - 60) +
                  .02 * (rhc$aps - 55) + rhc$death)

rhc_mis$pafi[runif(nrow(rhc)) < p_mis] <- NA
rhc_mis$crea[runif(nrow(rhc)) < p_mis] <- NA

colMeans(is.na(rhc_mis))
#>  death    age    sex    aps meanbp   pafi   crea   card 
#>  0.000  0.000  0.000  0.000  0.000  0.256  0.258  0.000
```