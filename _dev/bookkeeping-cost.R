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

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

TREES <- 50L
BURN <- 500L
DRAWS <- 500L
REPS <- 3L
SIZES <- c(500L, 2000L, 8000L)
OUT <- "_dev/bookkeeping-cost.rds"

# Five arms per size, each timed as the best of REPS.
N_FITS <- length(SIZES) * 5L * REPS

pr <- prog_init(total = N_FITS, title = "What the bookkeeping costs",
                unit = "fit", kind = "benchmark")
on.exit(prog_end(pr, "failed", "aborted before the last size"), add = TRUE)

timings <- list()
section <- ""
done <- 0L

# Written after every size rather than at the end, so a run that is killed
# leaves its finished sizes readable. `complete` is what tells a reader which of
# the two it is looking at.
checkpoint <- function(complete) {
  saveRDS(list(res = timings, complete = complete, done = done,
               total = N_FITS, reps = REPS, sizes = SIZES,
               num_trees = TREES), OUT)
}

friedman <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5]
}

time_it <- function(label, f) {
  best <- Inf

  for (r in seq_len(REPS)) {
    t0 <- proc.time()[["elapsed"]]
    f()
    secs <- proc.time()[["elapsed"]] - t0
    best <- min(best, secs)
    done <<- done + 1L
    prog_tick(pr, i = done, secs = secs,
              label = sprintf("%s / %s rep %d", section, label, r))
  }

  best
}

for (n in SIZES) {
  set.seed(1)
  p <- 10L
  x <- matrix(stats::runif(n * p), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  d <- data.frame(y = friedman(x) + stats::rnorm(n), x)

  ctrl <- function(...) {
    bartisan_control(num_trees = TREES, num_burn = BURN, num_draws = DRAWS,
                     gate = "hard", verbose = FALSE, ...)
  }

  section <- sprintf("n = %d", n)

  rows <- c(
    `bartisan hard` = time_it("bartisan hard", function() {
      set.seed(2); bartisan(y ~ ., data = d, family = gaussian(), control = ctrl())
    }),
    `bartisan, no sparsity draw` = time_it("no sparsity draw", function() {
      set.seed(2); bartisan(y ~ ., data = d, family = gaussian(),
                            control = ctrl(sparsity = FALSE))
    }),
    `bartisan, scale fixed` = time_it("scale fixed", function() {
      set.seed(2); bartisan(y ~ ., data = d, family = gaussian(),
                            control = ctrl(sparsity = FALSE,
                                           update_sigma_mu = FALSE))
    }),
    dbarts = time_it("dbarts", function() {
      set.seed(2)
      dbarts::bart(x.train = x, y.train = d$y, ntree = TREES, nskip = BURN,
                   ndpost = DRAWS, verbose = FALSE, keeptrees = FALSE)
    }),
    flexBART = time_it("flexBART", function() {
      set.seed(2)
      flexBART::flexBART(y ~ bart(.), train_data = d, M = TREES,
                         n.chains = 1L, nd = DRAWS, burn = BURN, verbose = FALSE)
    })
  )

  timings[[section]] <- rows
  checkpoint(FALSE)

  cat(sprintf("\nn = %d, p = %d, %d trees, %d + %d iterations, best of %d\n",
              n, p, TREES, BURN, DRAWS, REPS))
  for (nm in names(rows)) {
    cat(sprintf("  %-28s %6.2f s  %5.2fx dbarts\n", nm, rows[[nm]],
                rows[[nm]] / rows[["dbarts"]]))
  }
}

checkpoint(TRUE)
on.exit()
prog_end(pr, "done", sprintf("%d fits", done))
cat("\nwrote", OUT, "\n")
