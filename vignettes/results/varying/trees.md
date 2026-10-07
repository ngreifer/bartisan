
``` r
bartisan(death ~ age + sex + race + edu + aps +
           meanbp + resp + hema + pafi +
           paco2 + crea + surv2m + card +
           vc(rhc),
         data = rhc, family = binomial(), num_trees = c(50L, 25L))

bartisan(death ~ age + sex + race + edu + aps +
           meanbp + resp + hema + pafi +
           paco2 + crea + surv2m + card +
           vc(rhc),
         data = rhc, family = binomial(),
         num_trees = c("(Intercept)" = 50L, rhc = 25L))
```