
``` r
library(bartisan)

set.seed(2026)

# Small chains throughout, so that this vignette builds quickly. The defaults
# are 50 trees and 800 draws after 200 warmup iterations.
ctrl <- bartisan_control(num_trees = 10, num_burn = 150, num_draws = 150)
```