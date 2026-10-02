# Does the bandwidth move's step-size floor bind when the likelihood is sharp?
#
# What is tested. The soft-rule Poisson SBC arm deviates from uniform ranks,
# and `_dev/sbc.R` with a known offset added to both the generator and the fit
# shows the deviation growing with the level of the response rather than
# shrinking: the dispersion contrast is +2.8 at a mean count near one, +2.9 at
# seven and +4.3 at fifty-five, where coverage falls to 0.888 and the
# posterior for the calibrated contrast narrows from 1.46 to 0.36. That rules
# out the Laplace approximation's accuracy, which improves as a Poisson
# approaches a Gaussian. What tracks the damage is how narrow the posterior is.
#
# `update_bandwidth()` in `src/mcmc.cpp` proposes a multiplicative random walk,
# `b_new = b_old * exp(log_step * U(-1, 1))`, adapts `log_step` toward an
# acceptance of 0.44 by Robbins-Monro during warmup only, and clamps it below
# at `log(1.02)`. That floor is a two percent minimum step. The claim: where
# the likelihood is sharp enough that a tree's bandwidth posterior is narrower
# than a couple of percent, the floor binds, the move cannot propose inside the
# posterior, acceptance collapses and the chain sticks, which is why neither
# four times the draws nor ten times the warmup moved the deviation.
#
# What is measured. The SBC design of `_dev/sbc.R` -- 20 trees, n = 400, two
# predictors, `sigma_mu` fixed at 0.35, bandwidth prior mean 0.1,
# `x_transform = "range"`, 400 warmup and 2000 kept draws, one chain -- with a
# Poisson response at offsets of 0, 2 and 4, five replicates each. Per fit,
# over the stored per-tree bandwidths:
#
#   moved        the share of consecutive retained draws in which a tree's
#                bandwidth changed, median over trees. This is the acceptance
#                rate of the move as the stored draws see it, against the 0.44
#                the adaptation targets.
#   ess          the bandwidth's bulk effective sample size, median over trees
#   rel_sd       the bandwidth's posterior standard deviation over its mean,
#                median over trees. The floor is a step of 2%, so a posterior
#                whose relative spread is near or below that is one the move
#                cannot resolve.
#   stuck        the share of trees whose bandwidth never moved at all
#
# What each outcome would mean. Acceptance falling well below 0.44 as the
# offset rises, with `rel_sd` approaching or crossing two percent, identifies
# the clamp: a real defect with a contained fix, either lowering the floor or
# letting the adaptation run against a floor set from the data. Acceptance
# holding near its target with a comfortable `rel_sd` leaves the floor
# innocent and the move's internals as the thing to instrument.
#
# Run with: Rscript _dev/bandwidth-step-floor.R
# Writes:   _dev/bandwidth-step-floor.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

P <- 2L
N <- 400L
TREES <- 20L
SIGMA_MU <- 0.35
GAMMA <- 0.95
BETA <- 2
BANDWIDTH <- 0.1
HALF_WIDTH <- 4.055935661788187
REPS <- 5L

draw_tree <- function(depth = 0L, limits = NULL) {
  if (is.null(limits)) limits <- lapply(seq_len(P), function(i) c(0, 1))
  if (runif(1L) >= GAMMA * (1 + depth)^(-BETA)) {
    return(list(leaf = TRUE, mu = rnorm(1L, 0, SIGMA_MU)))
  }
  v <- sample.int(P, 1L)
  lim <- limits[[v]]
  val <- lim[1L] + (lim[2L] - lim[1L]) * runif(1L)
  left <- limits
  left[[v]] <- c(lim[1L], val)
  right <- limits
  right[[v]] <- c(val, lim[2L])
  list(leaf = FALSE, var = v, val = val,
       left = draw_tree(depth + 1L, left), right = draw_tree(depth + 1L, right))
}

left_prob <- function(x, val, bandwidth) {
  t <- 0.5 + 0.5 * (val - x) / (bandwidth * HALF_WIDTH)
  t <- pmin(pmax(t, 0), 1)
  t * t * (3 - 2 * t)
}

eval_tree <- function(node, u, bandwidth, w = rep.int(1, nrow(u))) {
  if (isTRUE(node$leaf)) return(node$mu * w)
  p <- left_prob(u[, node$var], node$val, bandwidth)
  eval_tree(node$left, u, bandwidth, w * p) +
    eval_tree(node$right, u, bandwidth, w * (1 - p))
}

draw_forest <- function(u) {
  Reduce(`+`, lapply(seq_len(TREES), function(i) {
    eval_tree(draw_tree(), u, rexp(1L, rate = 1 / BANDWIDTH))
  }), numeric(nrow(u)))
}

set.seed(20260905)
u <- matrix(runif(N * P), N, P)
colnames(u) <- paste0("x", seq_len(P))
d0 <- as.data.frame(u)

grid <- expand.grid(offset = c(0, 2, 4), rep = seq_len(REPS))

pr <- prog_init(total = nrow(grid), title = "Bandwidth move: the step floor",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

ess1 <- function(x) posterior::ess_bulk(matrix(x, ncol = 1L))

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]

  set.seed(5000L + g$rep)
  eta <- draw_forest(u)
  d <- d0
  d$y <- rpois(N, exp(eta + g$offset))

  fit <- bartisan(y ~ ., data = d, family = poisson(),
                  offset = if (g$offset != 0) rep(g$offset, N),
                  control = bartisan_control(num_trees = TREES, num_burn = 400L,
                                             num_draws = 2000L, chains = 1L,
                                             gate = "smoothstep",
                                             sigma_mu = SIGMA_MU,
                                             update_sigma_mu = FALSE,
                                             sparsity = FALSE,
                                             x_transform = "range",
                                             bandwidth = BANDWIDTH,
                                             augment = FALSE, verbose = FALSE))

  bw <- fit[["bandwidth"]]
  moved <- apply(bw, 2L, function(x) mean(diff(x) != 0))
  rel_sd <- apply(bw, 2L, function(x) sd(x) / mean(x))

  rows[[k]] <- data.frame(
    offset = g$offset, rep = g$rep, mean_count = mean(d$y),
    moved = median(moved), ess = median(apply(bw, 2L, ess1)),
    rel_sd = median(rel_sd), stuck = mean(moved == 0),
    bw_mean = median(bw))

  prog_tick(pr, label = sprintf("offset %g / rep %d", g$offset, g$rep))
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/bandwidth-step-floor.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- do.call(rbind, lapply(split(res, res$offset), function(x) {
  data.frame(offset = x$offset[1L], mean_count = round(mean(x$mean_count), 1),
             moved = round(median(x$moved), 3), ess = round(median(x$ess)),
             rel_sd = round(median(x$rel_sd), 4),
             stuck = round(mean(x$stuck), 3))
}))
cat("\nthe move's step floor is 2%; `rel_sd` is the posterior spread it must resolve\n")
print(out, row.names = FALSE)
