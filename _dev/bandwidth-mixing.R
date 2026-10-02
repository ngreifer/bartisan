# Does the per-tree bandwidth mix, under a Poisson likelihood and under a logit?
#
# What is tested. Pooled over 1200 independent replicates in two disjoint seed
# blocks, the soft-rule Poisson arm of `_dev/sbc.R` deviates from uniform
# ranks: chi-square 29.3 on 9 df (p = 0.0006), the top rank bin 4.3 standard
# errors high, the mean rank 2.5 standard errors above the middle, and coverage
# 0.942 against a nominal 0.95. Holding every tree's bandwidth fixed removes it
# (mean-rank z -0.20, top-bin z -0.27), the same arm under hard rules shows
# nothing, and the logit under soft rules with the bandwidth drawn shows
# nothing, which puts it in the per-tree bandwidth move under a Poisson
# likelihood. Reading that move leaves its acceptance ratio looking right: a
# multiplicative random walk whose Jacobian is `log(new / old)`, an exponential
# prior ratio, a likelihood difference, and a rollback from a snapshot. The
# claim to test is the benign alternative: that the move is correct and slow,
# so that the bandwidth is an under-explored nuisance dimension at this chain
# length rather than a wrongly targeted one.
#
# What is measured. The design of `_dev/sbc.R` -- 20 trees, n = 400, two
# predictors, `sigma_mu` fixed at 0.35, bandwidth prior mean 0.1, 400 warmup
# and 1000 kept draws -- for soft-rule Poisson and soft-rule logit, 25
# replicates each. `fit$bandwidth` holds one column per tree, so per fit:
#
#   ess_median    the median over trees of the bulk effective sample size of a
#                 tree's bandwidth over the 1000 retained draws
#   ess_min       the smallest of them
#   moved         the share of consecutive draws in which a tree's bandwidth
#                 changed, median over trees, which is the acceptance rate of
#                 the move as the stored draws see it
#   bw_median     the median drawn bandwidth, against a prior mean of 0.1
#
# What each outcome would mean. A Poisson bandwidth effective sample size far
# below the logit's says the move is correct and slow: the deviation is then an
# under-explored nuisance dimension, a longer chain should remove it, and the
# finding is a limitation to record rather than a defect to fix. Comparable
# effective sample sizes in the two arms would mean the move explores the
# bandwidth perfectly well and still leaves the posterior wrong under a Poisson
# likelihood, which makes the acceptance ratio or the likelihood difference it
# uses the thing to re-derive.
#
# Run with: Rscript _dev/bandwidth-mixing.R [replicates]
# Writes:   _dev/bandwidth-mixing.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) > 0L) as.integer(args[1L]) else 25L

P <- 2L
N <- 400L
TREES <- 20L
SIGMA_MU <- 0.35
GAMMA <- 0.95
BETA <- 2
BANDWIDTH <- 0.1
HALF_WIDTH <- 4.055935661788187

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

arms <- list(
  list(name = "poisson", family = poisson(),
       draw = function(eta) rpois(N, exp(eta))),
  list(name = "logit", family = binomial(),
       draw = function(eta) rbinom(N, 1L, plogis(eta)))
)

grid <- expand.grid(arm = seq_along(arms), rep = seq_len(REPS))

pr <- prog_init(total = nrow(grid), title = "Bandwidth mixing by family",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

ess1 <- function(x) posterior::ess_bulk(matrix(x, ncol = 1L))

rows <- list()
for (k in seq_len(nrow(grid))) {
  arm <- arms[[grid$arm[k]]]
  r <- grid$rep[k]

  set.seed(5000L + r)
  eta <- draw_forest(u)
  d <- d0
  d$y <- arm$draw(eta)

  if (length(unique(d$y)) < 2L) {
    prog_tick(pr, ok = FALSE, label = sprintf("%s rep %d", arm$name, r),
              msg = "response had no variation")
    next
  }

  fit <- bartisan(y ~ ., data = d, family = arm$family,
                  control = bartisan_control(num_trees = TREES, num_burn = 400L,
                                             num_draws = 1000L, chains = 1L,
                                             gate = "smoothstep",
                                             sigma_mu = SIGMA_MU,
                                             update_sigma_mu = FALSE,
                                             sparsity = FALSE,
                                             x_transform = "range",
                                             bandwidth = BANDWIDTH,
                                             augment = FALSE, verbose = FALSE))

  bw <- fit[["bandwidth"]]
  ess <- apply(bw, 2L, ess1)
  moved <- apply(bw, 2L, function(x) mean(diff(x) != 0))

  rows[[length(rows) + 1L]] <- data.frame(
    arm = arm$name, rep = r, ess_median = median(ess), ess_min = min(ess),
    moved = median(moved), bw_median = median(bw))
  prog_tick(pr, label = sprintf("%s rep %d", arm$name, r))
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/bandwidth-mixing.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- do.call(rbind, lapply(split(res, res$arm), function(x) {
  data.frame(arm = x$arm[1L], fits = nrow(x),
             ess_median = round(median(x$ess_median)),
             ess_min = round(median(x$ess_min)),
             moved = round(median(x$moved), 3),
             bw_median = round(median(x$bw_median), 3))
}))
print(out, row.names = FALSE)
