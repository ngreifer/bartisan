
``` r
# How much of the variation is between groups: the posterior of the group scale
quantile(fit_re$tau[[1L]][, "g"], c(.025, .5, .975))
#>   2.5%    50%  97.5% 
#> 0.4220 0.5221 0.6656

# The intercepts, as posterior means, against the values they were drawn from
re <- nlme::ranef(fit_re)$g
cor(re[, "(Intercept)"], u)
#> [1] 0.9613

# With intervals, for the first four groups
apply(nlme::ranef(fit_re, draws = TRUE)$g[["(Intercept)"]][, 1:4], 2L,
      quantile, c(.025, .975))
#>            1       2       3       4
#> 2.5%  0.3092 0.09548 -0.5970 -0.7757
#> 97.5% 0.8688 0.72409  0.0372 -0.2062
```