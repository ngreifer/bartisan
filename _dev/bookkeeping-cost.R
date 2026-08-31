# flexBART was faster than bartisan at n = 2000 in the categorical benchmark,
# even with hard rules on both sides. Their paper attributes their speed partly
# to caching which observations reach which leaf instead of looping over the
# whole dataset on every tree update. bartisan already does that -- `Node::idx`
# is that cache and `split_support()` divides a parent's cache between its
# children -- so the gap has to be somewhere else. This looks for it.
#
# The ablation is against dbarts, which is a conjugate-Gaussian implementation
# with the same hard rules and no generalized machinery, so the distance from it
# is the price of that machinery rather than of bookkeeping.

library(bartisan)

TREES <- 50L
BURN <- 500L
DRAWS <- 500L
REPS <- 3L

friedman <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5]
}

time_it <- function(f) {
  best <- Inf
  for (r in seq_len(REPS)) {
    t0 <- proc.time()[["elapsed"]]
    f()
    best <- min(best, proc.time()[["elapsed"]] - t0)
  }
  best
}

for (n in c(500L, 2000L, 8000L)) {
  set.seed(1)
  p <- 10L
  x <- matrix(stats::runif(n * p), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  d <- data.frame(y = friedman(x) + stats::rnorm(n), x)

  ctrl <- function(...) {
    bartisan_control(num_trees = TREES, num_burn = BURN, num_draws = DRAWS,
                     gate = "hard", verbose = FALSE, ...)
  }

  rows <- c(
    `bartisan hard` = time_it(function() {
      set.seed(2); bartisan(y ~ ., data = d, family = gaussian(), control = ctrl())
    }),
    `bartisan, no sparsity draw` = time_it(function() {
      set.seed(2); bartisan(y ~ ., data = d, family = gaussian(),
                            control = ctrl(sparsity = FALSE))
    }),
    `bartisan, scale fixed` = time_it(function() {
      set.seed(2); bartisan(y ~ ., data = d, family = gaussian(),
                            control = ctrl(sparsity = FALSE,
                                           update_sigma_mu = FALSE))
    }),
    dbarts = time_it(function() {
      set.seed(2)
      dbarts::bart(x.train = x, y.train = d$y, ntree = TREES, nskip = BURN,
                   ndpost = DRAWS, verbose = FALSE, keeptrees = FALSE)
    }),
    flexBART = time_it(function() {
      set.seed(2)
      flexBART::flexBART(y ~ bart(.), train_data = d, M = TREES,
                         n.chains = 1L, nd = DRAWS, burn = BURN, verbose = FALSE)
    })
  )

  cat(sprintf("\nn = %d, p = %d, %d trees, %d + %d iterations, best of %d\n",
              n, p, TREES, BURN, DRAWS, REPS))
  for (nm in names(rows)) {
    cat(sprintf("  %-28s %6.2f s  %5.2fx dbarts\n", nm, rows[[nm]],
                rows[[nm]] / rows[["dbarts"]]))
  }
}
