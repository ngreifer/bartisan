# Do R-hat and ESS behave the way the chains-versus-draws argument says they do,
# on a process where the answer is known?
#
# An AR(1) has an integrated autocorrelation time of (1 + rho) / (1 - rho), so
# the effective sample size of m chains of length N is mN divided by that, and it
# depends only on the product. Any dependence on how mN is split between m and N
# is an artifact of the estimator, not information. R-hat has no such invariance
# and is not supposed to: it is a statement about agreement, and short chains
# agree less.
#
# Two starts. From stationarity there is nothing for R-hat to find, so whatever
# it reports is its null behavior. From a dispersed start there is, and how fast
# it goes away is the thing the user is choosing between.
#
# Usage: Rscript _dev/ess-rhat-calibration.R [reps]
ROOT <- "/Users/NoahGreifer/Dropbox/Research/R/bartisan"
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

suppressMessages(library(bartisan))
args <- commandArgs(TRUE)
REPS <- if (length(args) > 0) as.integer(args[[1]]) else 200L
RHO  <- 0.99

# `filter()` rather than a loop, because the ladder below runs to 80000 draws a
# chain and an R-level recursion at that length would cost more than everything
# else here put together. The recursion is the same one: x[i] = rho x[i-1] + e,
# started at `start`, so a start drawn from N(0, 1) is a start from stationarity.
ar1 <- function(n, rho, start) {
  as.numeric(stats::filter(stats::rnorm(n, sd = sqrt(1 - rho^2)), rho,
                           method = "recursive", init = start))
}

tau <- (1 + RHO) / (1 - RHO)
cells <- data.frame(chains = c(2L, 4L, 8L, 16L), draws = c(8000L, 4000L, 2000L, 1000L))

OUT <- file.path(ROOT, "_dev/ess-rhat-calibration.rds")
LADDER <- c(25, 50, 100, 200, 400, 800, 1600)

# Two dispersions over the four splits, then the seven rungs of the ladder.
n_total <- (2L * nrow(cells) + length(LADDER)) * REPS

pr <- prog_init(total = n_total, unit = "replicate", kind = "simulation",
                title = "R-hat and ESS on a known process")
on.exit(prog_end(pr, "failed", "aborted before the last cell"), add = TRUE)

out <- list()
ladder <- list()
done <- 0L

# Written after every cell rather than at the end, so a run that is killed
# leaves its finished cells readable. `complete` is what tells a reader which of
# the two it is looking at.
checkpoint <- function(complete) {
  split <- do.call(rbind, out)
  if (!is.null(split)) {
    split[["ess_over_truth"]] <- split[["ess"]] / split[["truth"]]
  }
  saveRDS(list(split = split,
               ladder = if (length(ladder)) do.call(rbind, ladder),
               complete = complete, done = done, total = n_total,
               reps = REPS, rho = RHO), OUT)
}

set.seed(20260905)

for (disp in c(FALSE, TRUE)) {
  for (i in seq_len(nrow(cells))) {
    m <- cells[["chains"]][i]
    n <- cells[["draws"]][i]
    r <- replicate(REPS, {
      starts <- if (disp) stats::rnorm(m, sd = 3) else stats::rnorm(m)
      x <- vapply(seq_len(m), function(j) ar1(n, RHO, starts[j]), numeric(n))
      got <- bartisan:::diagnosis_stats(x)[c("rhat", "ess_bulk")]
      done <<- done + 1L
      prog_tick(pr, i = done,
                label = sprintf("%d chains x %d draws%s", m, n,
                                if (disp) ", dispersed" else ""))
      got
    })
    out[[length(out) + 1L]] <- data.frame(
      dispersed = disp, chains = m, draws = n, total = m * n,
      truth = m * n / tau,
      rhat = mean(r["rhat", ]), rhat_p90 = stats::quantile(r["rhat", ], 0.9, names = FALSE),
      over = mean(r["rhat", ] > 1.01), ess = mean(r["ess_bulk", ]))

    checkpoint(FALSE)
  }
}

res <- do.call(rbind, out)
res[["ess_over_truth"]] <- res[["ess"]] / res[["truth"]]

# How much effective sample size R-hat needs before 1.01 means anything. Four
# chains throughout, from stationarity, with the length chosen to put the true
# effective sample size on a ladder. If the two thresholds in `diagnose()` are
# calibrated to each other, the null R-hat should cross 1.01 at about 400.
for (target in LADDER) {
  n <- ceiling(target * tau / 4)
  r <- replicate(REPS, {
    x <- vapply(seq_len(4L), function(j) ar1(n, RHO, stats::rnorm(1L)), numeric(n))
    got <- bartisan:::diagnosis_stats(x)[c("rhat", "ess_bulk")]
    done <<- done + 1L
    prog_tick(pr, i = done, label = sprintf("ladder, true ESS %d", target))
    got
  })
  ladder[[length(ladder) + 1L]] <- data.frame(
    truth = target, draws = n, rhat = mean(r["rhat", ]),
    rhat_p90 = stats::quantile(r["rhat", ], 0.9, names = FALSE),
    over = mean(r["rhat", ] > 1.01), ess = mean(r["ess_bulk", ]))

  checkpoint(FALSE)
}

checkpoint(TRUE)
ladder <- do.call(rbind, ladder)

on.exit()
prog_end(pr, "done", sprintf("%d replicates", done))

cat(sprintf("\nAR(1) with rho = %.2f, so tau = %.0f and the true ESS of %d draws is %.0f.\n",
            RHO, tau, 16000L, 16000 / tau))
cat(sprintf("%d replicates per cell.\n\n", REPS))
for (disp in c(FALSE, TRUE)) {
  cat(sprintf("--- chains started %s\n", if (disp) "far apart (sd 3)" else "from stationarity"))
  cat(sprintf("%7s %7s | %8s %8s %8s | %8s %8s %8s\n",
              "chains", "draws", "mean", "90th pct", "P(>1.01)", "true", "est.", "est/true"))
  cat(sprintf("%7s %7s | %8s %8s %8s | %8s %8s %8s\n",
              "", "", "R-hat", "R-hat", "", "ESS", "ESS", ""))
  d <- res[res[["dispersed"]] == disp, ]
  for (i in seq_len(nrow(d)))
    cat(sprintf("%7d %7d | %8.3f %8.3f %7.0f%% | %8.0f %8.0f %8.2f\n",
                d[["chains"]][i], d[["draws"]][i], d[["rhat"]][i], d[["rhat_p90"]][i],
                100 * d[["over"]][i], d[["truth"]][i], d[["ess"]][i],
                d[["ess_over_truth"]][i]))
  cat("\n")
}
cat("--- four chains from stationarity, lengthened: what R-hat's null looks like\n")
cat(sprintf("%9s %8s | %8s %8s %8s %9s\n", "true ESS", "draws", "est. ESS", "mean", "90th pct", "P(>1.01)"))
cat(sprintf("%9s %8s | %8s %8s %8s %9s\n", "", "per chain", "", "R-hat", "R-hat", ""))
for (i in seq_len(nrow(ladder)))
  cat(sprintf("%9.0f %8d | %8.0f %8.3f %8.3f %8.0f%%\n", ladder[["truth"]][i],
              ladder[["draws"]][i], ladder[["ess"]][i], ladder[["rhat"]][i],
              ladder[["rhat_p90"]][i], 100 * ladder[["over"]][i]))

cat("\nDONE\n")
