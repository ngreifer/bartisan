
``` r
bartisan(death ~ age + sex + race + edu + aps +
           meanbp + resp + hema + pafi +
           paco2 + crea + surv2m + card + .propensity +
           vc(rhc, center = "estimate"),
         data = rhc, family = binomial(), num_trees = c(50L, 25L))
```