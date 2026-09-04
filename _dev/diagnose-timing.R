# How long the convergence pass takes, and what a `future` plan buys it.
#
# `diagnose()` walks one column per observation, computing a rank-normalized
# R-hat and two effective sample sizes for each. The columns are independent and
# no random numbers are drawn, so the pass splits cleanly over workers; what is
# not obvious is whether the split pays, because the draws have to reach the
# workers first (about 25 MB at n = 2000 with four chains of 400).
#
# This was never measured on real hardware: it was written on a machine where
# `parallelly::availableCores()` reports 1, so every worker contended for one
# core and the observed speedup of 1.16x said nothing. Run it somewhere with
# cores.
#
# Run with: Rscript _dev/diagnose-timing.R
# Writes:   _dev/diagnose-timing.rds

library(bartisan)

if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
  stop("this needs future and future.apply")
}

cat(sprintf("availableCores(): %d\n\n", parallelly::availableCores()))

sizes <- c(500L, 2000L, 8000L)
worker_counts <- c(1L, 2L, 4L, 8L)
reps <- 3L

# Only the counts the machine can actually run, so a reported speedup is not
# workers contending for one core.
worker_counts <- worker_counts[worker_counts <= parallelly::availableCores()]

sim <- function(n) {
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  x3 = stats::runif(n))
  d$y <- 2 * d$x1 + sin(pi * d$x2) + stats::rnorm(n)
  d
}

out <- list()

for (n in sizes) {
  set.seed(1)
  d <- sim(n)

  future::plan(future::sequential)

  fit_seconds <- system.time(
    fit <- bartisan(y ~ x1 + x2 + x3, data = d, chains = 4L,
                    family = gaussian(),
                    control = bartisan_control(num_trees = 50L,
                                               num_burn = 400L,
                                               num_draws = 400L,
                                               gate = "hard"))
  )[["elapsed"]]

  draws_mb <- as.numeric(utils::object.size(fit[["eta"]])) / 1e6

  for (workers in worker_counts) {
    if (workers == 1L) {
      future::plan(future::sequential)
    }
    else {
      future::plan(future::multisession, workers = workers)
    }

    # The first call under a new plan pays to start the workers, which is a cost
    # of the plan rather than of the pass, so it is not one of the replicates.
    invisible(diagnose(fit))

    seconds <- vapply(seq_len(reps), function(i) {
      system.time(diagnose(fit))[["elapsed"]]
    }, numeric(1L))

    out[[length(out) + 1L]] <- data.frame(
      n = n, workers = workers, fit_seconds = fit_seconds,
      draws_mb = draws_mb, seconds = stats::median(seconds),
      fastest = min(seconds), slowest = max(seconds)
    )

    cat(sprintf("n = %5d  workers = %d  diagnose() %5.2fs  (fit was %5.2fs, %5.1f MB of draws)\n",
                n, workers, stats::median(seconds), fit_seconds, draws_mb))
    utils::flush.console()
  }
}

future::plan(future::sequential)

results <- do.call(rbind, out)

saveRDS(results, "_dev/diagnose-timing.rds")

cat("\nspeedup against one worker, by size:\n")

for (n in sizes) {
  rows <- results[results$n == n, ]
  base <- rows$seconds[rows$workers == 1L]

  cat(sprintf("  n = %5d : %s\n", n,
              paste(sprintf("%dw %.2fx", rows$workers, base / rows$seconds),
                    collapse = "  ")))
}

cat("\nthe pass as a share of a four-chain fit, sequentially:\n")

for (n in sizes) {
  rows <- results[results$n == n & results$workers == 1L, ]
  cat(sprintf("  n = %5d : %.0f%% of the fit's own time\n", n,
              100 * rows$seconds / rows$fit_seconds))
}
