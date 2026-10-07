
``` r
# Forests: mean, mean:z, log_sd
bartisan(list(mean = y ~ x1 + x2 + vc(z),
              log_sd = ~ x1 + x2),
         data = d,
         family = gaussian_ls())

# Forests: mean, mean:z, log_sd, log_sd:z. One formula reaches every parameter,
# which is the rule every per-forest argument follows too
bartisan(y ~ x1 + x2 + vc(z), data = d, family = gaussian_ls())

# Each coefficient with modifiers of its own
bartisan(list(mean = y ~ x1 + x2 + vc(z, ~ x2),
              log_sd = ~ x1 + x2 + vc(z, ~ x1)),
         data = d, family = gaussian_ls())
```