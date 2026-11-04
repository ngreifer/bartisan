
``` r
mi_fits <- with(imp, glm(death ~ age + sex + aps + meanbp + pafi +
                           crea + card, family = binomial))

pooled <- summary(pool(mi_fits))

f <- death ~ age + sex + aps + meanbp + pafi + crea + card

cc_fit <- glm(f, data = rhc_mis, family = binomial)
full_fit <- glm(f, data = rhc, family = binomial)

vars <- c("pafi", "crea")

data.frame(
  term = vars,
  full = coef(full_fit)[vars],
  full_se = sqrt(diag(vcov(full_fit)))[vars],
  bart_mi = pooled$estimate[match(vars, pooled$term)],
  bart_mi_se = pooled$std.error[match(vars, pooled$term)],
  cc = coef(cc_fit)[vars],
  cc_se = sqrt(diag(vcov(cc_fit)))[vars],
  row.names = NULL
)
#>   term     full  full_se  bart_mi bart_mi_se       cc    cc_se
#> 1 pafi  0.00169 0.000554  0.00187   0.000643  0.00169 0.000689
#> 2 crea -0.07386 0.031159 -0.09623   0.034960 -0.05479 0.041331
```