
``` r
set.seed(2026)
fit_ph <- update(fit_dpm, family = ph())

set.seed(2026)
fit_ln <- update(fit_dpm, family = lognormal_aft())
```