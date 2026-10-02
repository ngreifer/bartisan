# Does the bandwidth move leave behind the predictor its tree encodes?
#
# What is tested. The soft-rule Poisson arm of `_dev/sbc.R` gives posteriors
# slightly too narrow whenever the bandwidth is drawn, and is uniform when it is
# fixed (`_dev/TASKS.md`, 2026-10-01 and 2026-10-02). The move's acceptance
# ratio reads correctly and its likelihood difference matches two full
# evaluations. What was never checked is its bookkeeping: `update_bandwidth()`
# in `src/mcmc.cpp` subtracts the tree from the predictor, rebuilds the gate
# supports at the proposed bandwidth, adds the tree back, and on a rejection
# restores the supports from a snapshot. The claim is that after every move the
# predictor the sampler carries equals a fresh evaluation of the stored trees at
# their stored bandwidths. It is false if the snapshot misses state, leaving a
# rejected tree's weights out of step with its bandwidth, or if the rebuilt
# predictor of an accepted move differs from a fresh evaluation. The engine
# only ever updates the predictor incrementally, so either error persists into
# every recorded draw after it.
#
# What is measured. The SBC design -- the same fixed design of n = 400 on two
# predictors, 20 trees, `sigma_mu` fixed at 0.35, bandwidth prior mean 0.1,
# `x_transform = "range"`, smoothstep gates, the bandwidth drawn, 400 warmup
# sweeps and 1000 kept, one chain -- with a Poisson response at two levels,
# mean counts near 1 and near 55 (the SBC's offsets 0 and 4), four fits each,
# every sweep recorded. Per fit: the largest absolute gap between the recorded
# predictor and `predict(type = "link", draws = TRUE)`, which replays the stored
# trees through the prediction path rather than the sampler's; and the number
# of accepted bandwidth moves among the kept draws, counted as changes in a
# tree's bandwidth between consecutive draws, to show the move was exercised.
#
# What each outcome would mean. Gaps at rounding (1e-10 or below) across
# thousands of accepted moves clear the move's bookkeeping, and the deviation
# lies in how the sampler treats a skewed Poisson target once the bandwidth
# moves rather than in the implementation. Gaps above that, growing with the
# number of accepted moves, put the bug in the move's bookkeeping.
#
# Run with: Rscript _dev/bandwidth-predictor.R
# Writes:   _dev/bandwidth-predictor.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

N <- 400L
set.seed(20260905)
u <- matrix(stats::runif(N * 2L), N, 2L)
colnames(u) <- c("x1", "x2")
d0 <- as.data.frame(u)

# A smooth truth with a range of about one on the log scale, as the SBC's
# forests have; its exact form does not matter to a bookkeeping check.
truth <- 0.6 * sin(2 * pi * d0$x1) + 0.4 * d0$x2

control <- bartisan_control(num_trees = 20L, num_burn = 400L, num_draws = 1000L,
                            chains = 1L, gate = "smoothstep", sigma_mu = 0.35,
                            update_sigma_mu = FALSE, sparsity = FALSE,
                            x_transform = "range", bandwidth = 0.1,
                            update_bandwidth = TRUE, augment = FALSE,
                            verbose = FALSE)

grid <- expand.grid(level = c(0, 4), rep = 1:4)

pr <- prog_init(total = nrow(grid), title = "Bandwidth move bookkeeping",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  set.seed(1000L * g$level + g$rep)
  d <- d0
  d$y <- stats::rpois(N, exp(g$level + truth))

  fit <- bartisan(y ~ ., data = d, family = stats::poisson(), control = control)

  recorded <- fit[["eta"]][[1L]]
  replayed <- stats::predict(fit, newdata = d0, type = "link", draws = TRUE)
  gap <- abs(as.vector(replayed) - as.vector(recorded))

  bw <- fit[["bandwidth"]]
  accepted <- sum(diff(bw) != 0)

  rows[[k]] <- data.frame(level = g$level, rep = g$rep,
                          mean_count = mean(d$y), max_gap = max(gap),
                          gap_last_draw = max(abs(replayed[nrow(replayed), ] -
                                                    recorded[nrow(recorded), ])),
                          accepted = accepted,
                          accept_rate = accepted / length(diff(bw)))
  prog_tick(pr, label = sprintf("level %g, fit %d: largest gap %.1e", g$level,
                                g$rep, max(gap)))
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/bandwidth-predictor.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
print(res, row.names = FALSE, digits = 3)
cat(sprintf("\nlargest gap over all fits %.2e, after %d accepted moves in all\n",
            max(res$max_gap), sum(res$accepted)))
