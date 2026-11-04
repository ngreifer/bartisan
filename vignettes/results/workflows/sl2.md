
``` r
set.seed(1234)
sl2 <- SuperLearner(Y = rhc$death, X = rhc[covs],
                    family = binomial(),
                    SL.library = c("SL.mean", "SL.glm",
                                   bart_learners$names),
                    cvControl = list(V = 5))

sl2
#> 
#> Call:  
#> SuperLearner(Y = rhc$death, X = rhc[covs], family = binomial(), SL.library = c("SL.mean",  
#>     "SL.glm", bart_learners$names), cvControl = list(V = 5)) 
#> 
#> 
#>                    Risk   Coef
#> SL.mean_All       0.226 0.0000
#> SL.glm_All        0.191 0.5239
#> SL.bartisan_1_All 0.192 0.0000
#> SL.bartisan_2_All 0.192 0.0000
#> SL.bartisan_3_All 0.191 0.4539
#> SL.bartisan_4_All 0.192 0.0223
```