
``` r
library(loo)

loo_fits <- loo_compare(list(dpm_aft   = loo(fit_dpm, scale = "time"),
                             lognormal = loo(fit_ln,  scale = "time"),
                             ph        = loo(fit_ph,  scale = "time")))

loo_fits
#>      model elpd_diff se_diff p_worse diag_diff       diag_elpd
#>         ph       0.0     0.0      NA           1 k_psis > 0.66
#>    dpm_aft     -36.9    13.1    1.00                          
#>  lognormal     -82.5    15.0    1.00
```