
``` r
bart_spec <- bart(trees = 50) |>
  set_engine("bartisan", family = "binomial") |>
  set_mode("classification")

bart_spec
#> BART Model Specification (classification)
#> 
#> Main Arguments:
#>   trees = 50
#> 
#> Engine-Specific Arguments:
#>   family = binomial
#> 
#> Computational engine: bartisan
```