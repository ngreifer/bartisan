
``` r
library(SuperLearner)

covs <- c("age", "sex", "race", "edu", "aps", "meanbp", "resp",
          "hema", "pafi", "paco2", "crea", "surv2m", "card")

set.seed(1234)
sl <- SuperLearner(Y = rhc$death, X = rhc[covs],
                   family = binomial(),
                   SL.library = c("SL.mean", "SL.glm", "SL.bartisan"),
                   cvControl = list(V = 5))

sl
#> 
#> Call:  
#> SuperLearner(Y = rhc$death, X = rhc[covs], family = binomial(), SL.library = c("SL.mean",  
#>     "SL.glm", "SL.bartisan"), cvControl = list(V = 5)) 
#> 
#> 
#>                  Risk  Coef
#> SL.mean_All     0.226 0.000
#> SL.glm_All      0.191 0.551
#> SL.bartisan_All 0.192 0.449
```