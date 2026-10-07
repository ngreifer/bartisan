
``` r
# The mean's control function and coefficient, then the log standard
# deviation's, in the order the print lists them
bartisan(y ~ x1 + x2 + vc(z), data = d, family = gaussian_ls(),
         num_trees = c(50L, 20L, 15L, 10L))

# The same, by name; the order then does not matter
bartisan(y ~ x1 + x2 + vc(z), data = d, family = gaussian_ls(),
         num_trees = c("mean"   = 50L, "mean:z"   = 20L,
                       "log_sd" = 15L, "log_sd:z" = 10L))
```