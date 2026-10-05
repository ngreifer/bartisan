
``` r
library(marginaleffects)

drawn <- avg_comparisons(fit, variables = "rhc") |>
  posterior_draws()

chains <- fit$chains
per_chain <- nrow(drawn) / chains

folded <- array(drawn$draw, dim = c(per_chain, chains, 1L),
                dimnames = list(NULL, NULL, "ATE"))

folded |>
  as_draws_array() |>
  summarise_draws("mean", "rhat", "ess_bulk", "ess_tail")
#> # A tibble: 1 × 5
#>   variable   mean  rhat ess_bulk ess_tail
#>   <chr>     <dbl> <dbl>    <dbl>    <dbl>
#> 1 ATE      0.0621  1.00    1872.    2645.
```