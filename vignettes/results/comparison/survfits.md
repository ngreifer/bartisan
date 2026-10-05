
``` r
library(survival)

surv_model <- Surv(days, death) ~ age + sex + race + edu +
                 aps + meanbp + surv2m

set.seed(2026)
aft <- bartisan(surv_model, data = rhc, family = lognormal_aft())

set.seed(2026)
prop_haz <- bartisan(surv_model, data = rhc, family = ph())

loo_compare(list(aft = loo(aft,      scale = "time"),
                 ph  = loo(prop_haz, scale = "time")))
#>  model elpd_diff se_diff p_worse diag_diff       diag_elpd
#>     ph       0.0     0.0      NA           1 k_psis > 0.66
#>    aft     -90.0    15.4    1.00
```