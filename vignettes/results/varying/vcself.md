
``` r
# The effect of `aps` may itself change across `aps`, so the dose response is
# a curve rather than a line through the origin
bartisan(death ~ age + sex + vc(aps, ~ . + aps),
         data = rhc, family = binomial())
```