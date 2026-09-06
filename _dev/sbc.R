# Simulation-based calibration: does the sampler draw from the posterior it
# defines?
#
# Run with: Rscript _dev/sbc.R [replicates]
# Writes:   _dev/sbc.rds
#
# ---- what this can and cannot test -------------------------------------------
#
# SBC draws a parameter from the prior, simulates data from it, fits, and records
# where the truth falls among the posterior draws. Over many replicates those
# ranks are uniform if and only if the sampler targets the right posterior. It is
# the only check here that tests the *shape* of a posterior rather than where its
# center lands, which is what an interval is.
#
# It needs the generating prior to be the model's prior, and this model's prior
# is **empirical**, as every BART implementation's is:
#
#   * the leaf scale comes from the response's spread, though `sigma_mu` in
#     `bartisan_control()` fixes it, which is what this uses;
#   * the residual scale for a Gaussian is `residual_scale(y, x)`, which nothing
#     fixes, so this uses a `binomial()` instead, which has no such parameter;
#   * the predictor is centered at an empirical value, `qlogis(mean(y))` for a
#     binomial, which no prior of the model's own produced.
#
# The last one is why the target here is a **contrast**, one observation's
# predictor minus another's. The centering is common to both and cancels, so the
# quantity being calibrated is one the empirical part of the prior does not
# touch. Calibrating the level instead would be measuring the centering, which
# absorbs it entirely and would tell us nothing.
#
# The tree prior itself is data-free and is replicated exactly below: branch with
# probability `gamma * (1 + depth)^-beta`, split on a variable drawn uniformly,
# and take a cutpoint uniform on what the nearest ancestor splitting the same
# variable left available, starting from the whole of [0, 1]. Predictors are
# mapped to [0, 1] before fitting, so the simulation works on that scale
# directly and the mapping introduces nothing.

library(bartisan)

# `Rscript _dev/sbc.R <reps> <n> <gate>`. Naming an n and a gate runs that cell
# and writes its own file, so the grid can go to the queue at once.
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) > 0L) as.integer(args[1L]) else 200L
N <- if (length(args) > 1L) as.integer(args[2L]) else 400L
GATE <- if (length(args) > 2L) args[3L] else "hard"

P <- 2L
TREES <- 20L
SIGMA_MU <- 0.35
GAMMA <- 0.95
BETA <- 2
BANDWIDTH <- 0.1   # the mean of the exponential prior on a tree's bandwidth
L <- 100L      # thinned draws per fit, so a rank is one of 0..L

SOFT <- !identical(GATE, "hard")
HALF_WIDTH <- 4.055935661788187   # smoothstep, from `node.h`

# ---- the prior ---------------------------------------------------------------

# One tree, as a nested list. `limits` carries what the ancestors on each
# variable have left, which is exactly what `Node::get_limits()` reconstructs by
# walking up to the nearest ancestor that split the same variable.
draw_tree <- function(depth = 0L, limits = NULL) {
  if (is.null(limits)) {
    limits <- lapply(seq_len(P), function(i) c(0, 1))
  }

  if (stats::runif(1L) >= GAMMA * (1 + depth)^(-BETA)) {
    return(list(leaf = TRUE, mu = stats::rnorm(1L, 0, SIGMA_MU)))
  }

  v <- sample.int(P, 1L)
  lim <- limits[[v]]
  val <- lim[1L] + (lim[2L] - lim[1L]) * stats::runif(1L)

  left <- limits
  left[[v]] <- c(lim[1L], val)
  right <- limits
  right[[v]] <- c(val, lim[2L])

  list(leaf = FALSE, var = v, val = val,
       left = draw_tree(depth + 1L, left),
       right = draw_tree(depth + 1L, right))
}

# `left_prob()` from `node.h`. A hard rule sends `u <= val` left; a smoothstep
# gate sends a fraction of it, and every observation reaches every leaf with a
# weight, which is why the soft branch below carries weights rather than a
# partition.
left_prob <- function(x, val, bandwidth) {
  if (!SOFT) {
    return(as.numeric(x <= val))
  }

  t <- 0.5 + 0.5 * (val - x) / (bandwidth * HALF_WIDTH)
  t <- pmin(pmax(t, 0), 1)
  t * t * (3 - 2 * t)
}

# The tree's contribution, as a weighted sum over its leaves. `w` is how much of
# each observation has reached this node.
eval_tree <- function(node, u, bandwidth, w = rep.int(1, nrow(u))) {
  if (isTRUE(node$leaf)) {
    return(node$mu * w)
  }

  p <- left_prob(u[, node$var], node$val, bandwidth)

  eval_tree(node$left, u, bandwidth, w * p) +
    eval_tree(node$right, u, bandwidth, w * (1 - p))
}

draw_forest <- function(u) {
  # One bandwidth per tree, from the exponential prior the bandwidth move uses.
  Reduce(`+`, lapply(seq_len(TREES), function(i) {
    b <- if (SOFT) stats::rexp(1L, rate = 1 / BANDWIDTH) else 0
    eval_tree(draw_tree(), u, b)
  }), numeric(nrow(u)))
}

# ---- the run -----------------------------------------------------------------

# The predictors are fixed across replicates: SBC calibrates the parameter given
# the design, and re-drawing the design each time only adds noise.
set.seed(20260905)
u <- matrix(stats::runif(N * P), N, P)
colnames(u) <- paste0("x", seq_len(P))
d0 <- as.data.frame(u)

# Two observations far apart in the first predictor, so their contrast is a
# quantity the forest actually has to represent.
A <- which.min(u[, 1L])
B <- which.max(u[, 1L])

control <- bartisan_control(num_trees = TREES, num_burn = 400L,
                            num_draws = 1000L, chains = 1L, gate = GATE,
                            sigma_mu = SIGMA_MU, update_sigma_mu = FALSE,
                            sparsity = FALSE, x_transform = "range",
                            bandwidth = BANDWIDTH, augment = FALSE)

ranks <- integer(0)
truths <- widths <- covered <- numeric(0)

for (r in seq_len(reps)) {
  set.seed(5000L + r)

  eta <- draw_forest(u)
  truth <- eta[A] - eta[B]

  d <- d0
  d$y <- stats::rbinom(N, 1L, stats::plogis(eta))

  # A response with no variation carries no likelihood and the fit refuses it.
  if (length(unique(d$y)) < 2L) next

  fit <- bartisan(y ~ ., data = d, family = stats::binomial(), control = control)

  e <- fit[["eta"]][[1L]]
  contrast <- e[, A] - e[, B]
  thin <- contrast[seq(1L, length(contrast), length.out = L)]

  ranks <- c(ranks, sum(thin < truth))
  truths <- c(truths, truth)
  ci <- stats::quantile(contrast, c(0.025, 0.975), names = FALSE)
  widths <- c(widths, ci[2L] - ci[1L])
  covered <- c(covered, as.numeric(ci[1L] <= truth && truth <= ci[2L]))

  if (r %% 25L == 0L) {
    cat(sprintf("  %d of %d\n", r, reps))
    utils::flush.console()
  }
}

out <- data.frame(rank = ranks, truth = truths, width = widths,
                  covered = covered, n = N, gate = GATE)
saveRDS(out, sprintf("_dev/sbc-%s-%d.rds", GATE, N))

# ---- the report --------------------------------------------------------------

bins <- 10L
counts <- table(cut(out$rank, breaks = seq(0, L + 1, length.out = bins + 1L),
                    include.lowest = TRUE))
expected <- nrow(out) / bins
chisq <- sum((as.numeric(counts) - expected)^2) / expected

cat(sprintf("\nn = %d, %s rules: %d replicates, %d thinned draws each\n\n",
            N, GATE, nrow(out), L))
cat("rank histogram, 10 bins (uniform is what a correct sampler gives):\n")
cat(sprintf("  %s\n", paste(sprintf("%4d", as.numeric(counts)), collapse = "")))
cat(sprintf("  expected %.1f per bin\n", expected))
cat(sprintf("\nchi-square %.1f on %d df, p = %.3f\n", chisq, bins - 1L,
            stats::pchisq(chisq, bins - 1L, lower.tail = FALSE)))

# A U shape is a posterior too narrow, a hump in the middle one too wide, and a
# slope is a bias. Reporting the three separately says which, where a single
# p-value does not.
mid <- mean(out$rank > L * 0.25 & out$rank < L * 0.75)
cat(sprintf("\nshape: %.0f%% of ranks in the middle half (50%% is uniform)\n",
            100 * mid))
cat(sprintf("       mean rank %.1f of %d (%.1f is uniform)\n",
            mean(out$rank), L, L / 2))
cat(sprintf("\n95%% interval coverage of the contrast: %.3f (SE %.3f)\n",
            mean(out$covered),
            sqrt(mean(out$covered) * (1 - mean(out$covered)) / nrow(out))))
