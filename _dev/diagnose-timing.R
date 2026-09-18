# How long the convergence pass takes, and what a `future` plan buys it.
#
# `diagnose()` walks one column per observation, computing a rank-normalized
# R-hat and two effective sample sizes for each. The columns are independent and
# no random numbers are drawn, so the pass splits cleanly over workers; what is
# not obvious is whether the split pays, because the draws have to reach the
# workers first (about 25 MB at n = 2000 with four chains of 400).
#
# Run on ten cores it gives about 2.6x, 2.8x and 3.5x on eight workers at
# n = 500, 2000 and 8000, and the pass is 53% to 77% of a four-chain fit's own
# time when run sequentially. The ceiling is the number of fast cores rather
# than a serial section: that machine is an M4, four performance cores and six
# efficiency ones, and per-worker time is nearly flat to four workers and climbs
# after. A compute-bound loop over 8 KB of data gives the same curve, so it is
# not the size of the draws. Nothing in `diagnose()` is left to remove, and the
# same script should scale further on more than four equal cores.
#
# Needs a machine with cores. On one where `parallelly::availableCores()` reports
# 1 every worker contends for the same core and the numbers say nothing.
#
# Run with: Rscript _dev/diagnose-timing.R
# Writes:   _dev/diagnose-timing.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

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

OUT <- "_dev/diagnose-timing.rds"

# One fit per size, then a pass of `reps` diagnose() calls per worker count.
n_total <- length(sizes) * (1L + length(worker_counts))

pr <- prog_init(total = n_total, title = "diagnose() across worker counts",
                unit = "step", kind = "benchmark")
on.exit(prog_end(pr, "failed", "aborted before the last size"), add = TRUE)

out <- list()
done <- 0L

# Written after every worker count rather than at the end, so a run that is
# killed leaves its finished cells readable. `complete` is what tells a reader
# which of the two it is looking at.
checkpoint <- function(complete) {
  saveRDS(list(res = do.call(rbind, out), complete = complete, done = done,
               total = n_total, reps = reps, sizes = sizes,
               worker_counts = worker_counts), OUT)
}

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

  done <- done + 1L
  prog_tick(pr, i = done, secs = fit_seconds, label = sprintf("fit n=%d", n))

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

    done <- done + 1L
    prog_tick(pr, i = done, secs = sum(seconds),
              label = sprintf("diagnose() n=%d workers=%d", n, workers))
    checkpoint(FALSE)
  }
}

future::plan(future::sequential)

checkpoint(TRUE)
results <- do.call(rbind, out)

on.exit()
prog_end(pr, "done", sprintf("%d cells", nrow(results)))

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
