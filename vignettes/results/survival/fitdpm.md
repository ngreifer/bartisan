
``` r
set.seed(2026)

# Fit a BART accelerated failure time model whose error
# distribution is a Dirichlet process mixture of normals
fit_dpm <- bartisan(Surv(days, death) ~ rhc + age + sex + race + edu + aps +
                      meanbp + resp + hema + pafi + paco2 + crea + surv2m +
                      card,
                    data = rhc, family = dpm_aft())

fit_dpm
#> Generalized BART
#> 
#> Call:
#> bartisan(formula = Surv(days, death) ~ rhc + age + sex + race + 
#>     edu + aps + meanbp + resp + hema + pafi + paco2 + crea + 
#>     surv2m + card, data = rhc, family = dpm_aft())
#> 
#> Family: "dpm_aft" with the "identity" link
#> Observations: 1500
#> Structure: 1 forest of 50 trees, soft decision rules
#> Draws: 800 kept after 200 warmup
#> 
#> Posterior means: alpha = 0.764, clusters = 5.2, center = 0.464, error_sd = 1.91
```