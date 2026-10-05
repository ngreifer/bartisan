
``` r
library(bayesplot)

as_draws(fit, eta = 1) |>
  mcmc_trace(pars = c("loglik", "eta[1]"))
```

<img src="results/diagnostics/trace-1.png" alt="" style="display: block; margin: auto;" />