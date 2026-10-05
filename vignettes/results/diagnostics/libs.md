
``` r
library(bartisan)

# For parallelization; optional
if (rlang::is_installed("future")) {
  future::plan(future::multisession)
}

data("rhc")

set.seed(2026)

fit <- bartisan(death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
                  hema + pafi + paco2 + crea + surv2m + card,
                data = rhc, family = binomial(), chains = 4)
```