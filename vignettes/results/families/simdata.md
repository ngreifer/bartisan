
``` r
# The count response used here and through the rest of the vignette
n <- 300
d <- data.frame(x1 = runif(n), x2 = runif(n))
d$count <- rpois(n, exp(1.2 * sin(pi * d$x1) + 0.4))
```