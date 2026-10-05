
``` r
rstantools::prior_summary(fit)
#> Priors
#> 
#> Trees
#> • 50 trees per additive predictor, summed. A node at depth d branches with
#>   probability 0.95 * (1 + d)^-2, so the root splits with probability 0.95 and a
#>   node at depth 3 with 0.059.
#> 
#> Leaves
#> • Each leaf value is Normal(0, 0.212^2), that scale being 3 * s / (2 *
#>   sqrt(50)) with s the response's scale on the link scale. The scale is itself
#>   given a half-Cauchy prior centred there and is estimated.
#> 
#> Splitting variables
#> • Each of the 14 predictors is equally likely to be split on, and that is not
#>   drawn (`sparsity = FALSE`).
#> 
#> Decision rules
#> • Soft, with "smoothstep" gates. Each tree's bandwidth is drawn from an
#>   exponential with mean 0.1, on predictors mapped to [0, 1], and is estimated.
#> 
#> Family: binomial, logit link
#> • No parameters of its own beyond the additive predictors above.
#> 
#> ℹ The leaf scale, and any number above read off the response, are calibrated
#>   rather than fitted; that is how a BART prior is specified.
#> ℹ `prior_only = TRUE` in `bartisan()` draws from all of this, so that what it
#>   implies can be read on the outcome's own scale.
```