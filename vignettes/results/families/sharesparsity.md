
``` r
set.seed(9)
dw <- data.frame(matrix(runif(400 * 25), 400, 25))
names(dw) <- paste0("x", 1:25)
dw$y <- rnorm(400, 2 + 1.4 * sin(pi * dw$x1 * dw$x2) + 0.9 * dw$x4, 0.5)

apart <- bartisan(y ~ ., dw, family = gaussian_ls(), control = ctrl,
                  sparsity = TRUE)
shared <- bartisan(y ~ ., dw, family = gaussian_ls(), control = ctrl,
                   sparsity = TRUE, share_sparsity = TRUE)

top3 <- function(fit, forest) {
  imp <- variable_importance(fit)
  part <- imp[imp$predictor == forest, ]
  paste(head(part$variable[order(-part$prop_splits)], 3), collapse = " ")
}

# the spread is constant here, so the scale forest has nothing of its own to
# find and whatever it concentrates on came from the mean forest
c(apart = top3(apart, "log_sd"), shared = top3(shared, "log_sd"))
#>        apart       shared 
#> "x8 x14 x17"   "x1 x4 x2"
```