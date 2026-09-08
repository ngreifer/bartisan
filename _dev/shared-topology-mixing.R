# Is the speedup real, or does the chain just need more of it?
#
# Sharing a topology means a sweep performs one tree move per shared group
# instead of one per component, so a sweep is cheaper -- but it is also a
# smaller step, and a fit that runs twice as fast and needs twice as many draws
# has not gained anything. The number that settles it is effective sample size
# per second, not seconds per fit, so this script measures both at a matched
# number of draws.
#
# Reported on the worst-mixing quantities rather than the average, which is the
# convention the rest of this package's mixing work uses: an average over
# observations mixes about as well as the draw count whatever the sampler is
# doing, and the quantities that actually hold a fit back are the per-unit ones.
# Run on a quiet machine, after `_dev/shared-topology-timing.R`.

library(bartisan)

# Sharing the topology is not an argument: it is engaged by the package's
# internal flag, for the reasons recorded beside it in `R/control.R`. These
# scripts are the record of how it was measured, so they reach it the same way
# the tests do.
share_forests <- function(on) {
  flags <- bartisan:::the
  flags$share_forests <- isTRUE(on)
}

# Live progress, per the `/live-progress` skill. One tick per fit rather than
# per row: a row here is six four-chain fits and takes twenty minutes at
# n = 2000, so a per-row tick would leave the estimate blind for most of the run.
source(path.expand(Sys.getenv(
  "LIVE_PROGRESS", "~/.claude/skills/live-progress/assets/progress.R")))

REPS <- 3L
CHAINS <- 4L
NUM_BURN <- 1000L
NUM_DRAWS <- 1000L
OUT <- file.path("_dev", "shared-topology-mixing.rds")

sim_ls <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  d <- as.data.frame(X)
  d$y <- 3 * sin(pi * d$x1) + 2 * d$x2 +
    rnorm(n, 0, exp(-0.9 + 0.9 * d$x1 + 0.7 * d$x2))
  d
}

sim_zi <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  d <- as.data.frame(X)
  d$y <- ifelse(runif(n) < plogis(-0.6 + 1.6 * (0.9 * d$x1 + 0.7 * d$x2)),
                0L, rpois(n, exp(0.4 + 1.2 * d$x1)))
  d
}

# The quantities `diagnose()` reports, reduced to the ones that hold a fit back:
# the smallest bulk effective sample size, and the largest R-hat. `diagnose()`
# reports R-hat against a null that depends on the chain count and the
# effective sample size, so `ess_frac` is carried along to make the R-hat
# readable rather than compared with 1.01.
worst <- function(fit) {
  tab <- diagnose(fit)[["table"]]

  # The averages over observations are dropped: an average mixes about as well
  # as the draw count whatever the sampler is doing, so keeping it in a minimum
  # would only dilute the comparison.
  tab <- tab[!grepl("average over observations", tab[["quantity"]],
                    fixed = TRUE), ]

  c(ess_min = min(tab[["ess_bulk"]], na.rm = TRUE),
    ess_median = stats::median(tab[["ess_bulk"]], na.rm = TRUE),
    rhat_max = max(tab[["rhat"]], na.rm = TRUE))
}

one <- function(model, n, p, shared, seed) {
  d <- if (identical(model, "gaussian_ls")) sim_ls(n, p, seed) else
    sim_zi(n, p, seed)
  fam <- if (identical(model, "gaussian_ls")) gaussian_ls() else zi_poisson()
  share_forests(shared)
  ctrl <- bartisan_control(num_trees = 50L, num_burn = NUM_BURN,
                           num_draws = NUM_DRAWS, chains = CHAINS)
  t <- system.time(fit <- bartisan(y ~ ., data = d, family = fam,
                                   control = ctrl))
  c(elapsed = t[["elapsed"]], worst(fit))
}

quiet <- function(expr) {
  invisible(utils::capture.output(out <- suppressMessages(expr)))
  out
}

grid <- expand.grid(model = c("gaussian_ls", "zi_poisson"),
                    n = c(400L, 2000L), stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid) * 2L * REPS, title = "Shared forests: ESS per second",
                unit = "fit", kind = "benchmark")
on.exit(prog_end(pr, "failed"), add = TRUE)
fit_no <- 0L

res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]

  do.call(rbind, lapply(c(FALSE, TRUE), function(shared) {
    got <- do.call(rbind, lapply(seq_len(REPS), function(r) {
      out <- quiet(one(g$model, g$n, 20L, shared, seed = 6000L + 31L * i + r))
      fit_no <<- fit_no + 1L
      prog_tick(pr, i = fit_no, secs = out[["elapsed"]],
                label = sprintf("%s n=%d %s rep %d", g$model, g$n,
                                if (shared) "shared" else "separate", r))
      out
    }))

    out <- data.frame(model = g$model, n = g$n, shared = shared,
                      elapsed = mean(got[, "elapsed"]),
                      ess_min = mean(got[, "ess_min"]),
                      ess_median = mean(got[, "ess_median"]),
                      rhat_max = max(got[, "rhat_max"]),
                      stringsAsFactors = FALSE)
    out$ess_per_sec <- out$ess_min / out$elapsed

    cat(sprintf("%-12s n=%4d shared=%-5s  %6.1fs  ess_min %6.0f  ess/s %5.1f  rhat_max %.3f\n",
                out$model, out$n, out$shared, out$elapsed, out$ess_min,
                out$ess_per_sec, out$rhat_max))
    out
  }))
}))

saveRDS(list(res = res, chains = CHAINS, reps = REPS, num_burn = NUM_BURN,
             num_draws = NUM_DRAWS), OUT)

on.exit()
prog_end(pr, "done")

cat("\nwrote", OUT, "\n")
print(res)
