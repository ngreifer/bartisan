
``` r
set.seed(11)

nz <- 1000
dz <- data.frame(x1 = runif(nz), x2 = runif(nz), x3 = runif(nz))

# No zero inflation anywhere in this: one Poisson draw per observation
dz$count <- rpois(nz, exp(-3 + 5 * dz$x1 + 0.5 * dz$x2))

c(mean = mean(dz$count), zeros = mean(dz$count == 0),
  var_over_mean = var(dz$count) / mean(dz$count))
#>          mean         zeros var_over_mean 
#>          1.94          0.46          4.56
```