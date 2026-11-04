
``` r
library(tidymodels)

rhc_tm <- rhc
rhc_tm$died <- factor(rhc$death, levels = c(1, 0),
                      labels = c("yes", "no"))

set.seed(1234)
rhc_split <- initial_split(rhc_tm, prop = .75, strata = died)
rhc_train <- training(rhc_split)
rhc_test <- testing(rhc_split)
```