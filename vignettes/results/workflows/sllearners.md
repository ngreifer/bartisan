
``` r
bart_learners <- create.Learner("SL.bartisan",
                                tune = list(num_trees = c(20, 50),
                                            k = c(2, 4)))

bart_learners$names
#> [1] "SL.bartisan_1" "SL.bartisan_2" "SL.bartisan_3" "SL.bartisan_4"
```