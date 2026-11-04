
``` r
rhc_train$surv <- survival::Surv(rhc_train$days, rhc_train$death)
rhc_test$surv <- survival::Surv(rhc_test$days, rhc_test$death)
```