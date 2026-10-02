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

source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

library(bartisan)

# `Rscript _dev/sbc.R <reps> <n> <gate> <family>`. Naming an n and a gate runs
# that cell and writes its own file, so the grid can go to the queue at once.
#
# The family (added 2026-10-01) is one of "logit", the original and the
# default; "probit", which fits binomial("probit") with its augmentation on,
# so that the latent-variable sampler rather than the Laplace one is what is
# calibrated; and "poisson". Each has a likelihood with no parameter beyond the
# predictor, so the prior is the tree prior alone and SBC is well posed; the
# Gaussian's residual scale and the negative binomial's dispersion are
# calibrated from the data and are not. The target stays the contrast between
# the two extreme observations, for the reason above.
#
# What the two new arms test: that the probit augmentation and the Poisson
# target draw from the posterior they define, which the identity and recovery
# matrices cannot say, since an interval that is too narrow or too wide passes
# both. Uniform ranks and 95% coverage near 0.95 say they do; a U-shaped
# histogram says the posterior is too narrow, a hump too wide, a slope a bias.
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) > 0L) as.integer(args[1L]) else 200L
N <- if (length(args) > 1L) as.integer(args[2L]) else 400L
GATE <- if (length(args) > 2L) args[3L] else "hard"
FAMILY <- if (length(args) > 3L) args[4L] else "logit"

# `SBC_SEED_BASE` shifts the replicate seeds, for an independent replication of
# a run that has already been read; the default block is 5000. `SBC_FIX_BANDWIDTH`
# holds every tree's bandwidth at `BANDWIDTH` in the generator and switches
# `update_bandwidth` off in the fit, which removes the one piece of machinery
# that exists only under soft rules. `SBC_TAG` names the output file.
SEED_BASE <- as.integer(Sys.getenv("SBC_SEED_BASE", "5000"))
# `SBC_DRAWS` lengthens the chain without changing anything else, which tells a
# nuisance dimension that is merely under-explored from a target that is wrong.
DRAWS <- as.integer(Sys.getenv("SBC_DRAWS", "1000"))
# `SBC_WARMUP` lengthens warmup without changing the retained draws, which
# separates a chain that starts the retained draws in the wrong place from one
# whose draws are merely autocorrelated.
WARMUP <- as.integer(Sys.getenv("SBC_WARMUP", "400"))
# `SBC_OFFSET` adds a known constant to the generating predictor and supplies
# the same constant to the fit, which leaves the model unchanged and raises the
# level of the response. For a Poisson it is how far the log density is pushed
# toward a quadratic: at a mean near one it is sharply skewed and at a mean of
# fifty it is nearly normal, so the Laplace approximation the tree moves accept
# on gets better as this grows.
OFFSET <- as.numeric(Sys.getenv("SBC_OFFSET", "0"))
FIX_BANDWIDTH <- nzchar(Sys.getenv("SBC_FIX_BANDWIDTH"))
TAG <- Sys.getenv("SBC_TAG", "")

family <- switch(FAMILY,
                 logit = stats::binomial(),
                 probit = stats::binomial("probit"),
                 poisson = stats::poisson(),
                 stop("family must be logit, probit or poisson"))

draw_response <- switch(FAMILY,
                        logit = function(eta) stats::rbinom(N, 1L, stats::plogis(eta)),
                        probit = function(eta) stats::rbinom(N, 1L, stats::pnorm(eta)),
                        poisson = function(eta) stats::rpois(N, exp(eta)))

P <- 2L
TREES <- 20L
SIGMA_MU <- 0.35
GAMMA <- 0.95
BETA <- 2
# The mean of the exponential prior on a tree's bandwidth. `SBC_BANDWIDTH`
# varies it, which is the valid way to ask where in the prior a deviation
# lives: each value is a model of its own and so a marginal test of its own,
# where conditioning one run's ranks on the bandwidth it drew is not, since
# SBC's uniformity is marginal over the prior and conditioning on a component
# of the generating parameter breaks it for a correct sampler too.
BANDWIDTH <- as.numeric(Sys.getenv("SBC_BANDWIDTH", "0.1"))
# Thinned draws per fit, so a rank is one of 0..L. SBC wants the posterior
# draws independent of one another: a chain thinned by less than its own
# autocorrelation time gives a U-shaped histogram whatever the sampler does,
# so `SBC_L` exists to thin harder than the default 10-fold.
L <- as.integer(Sys.getenv("SBC_L", "100"))

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

# The bandwidths the last forest drew, kept so that a replicate's rank can be
# read against how smooth its truth was. Stashing them consumes no draw.
LAST_BW <- numeric(0)

draw_forest <- function(u) {
  # One bandwidth per tree, from the exponential prior the bandwidth move uses.
  bw <- numeric(TREES)

  out <- Reduce(`+`, lapply(seq_len(TREES), function(i) {
    b <- {
      if (!SOFT) 0
      else if (FIX_BANDWIDTH) BANDWIDTH
      else stats::rexp(1L, rate = 1 / BANDWIDTH)
    }
    bw[i] <<- b
    eval_tree(draw_tree(), u, b)
  }), numeric(nrow(u)))

  LAST_BW <<- bw
  out
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

control <- bartisan_control(num_trees = TREES, num_burn = WARMUP,
                            num_draws = DRAWS, chains = 1L, gate = GATE,
                            sigma_mu = SIGMA_MU, update_sigma_mu = FALSE,
                            sparsity = FALSE, x_transform = "range",
                            bandwidth = BANDWIDTH,
                            update_bandwidth = !FIX_BANDWIDTH,
                            augment = identical(FAMILY, "probit"))

ranks <- integer(0)
truths <- widths <- covered <- numeric(0)

# The seed block is part of the path whenever it is not the default, because a
# run under a different block is a different run: on 2026-10-02 one was
# launched without a tag and overwrote a finished block's results, which were
# recoverable only because every replicate seeds itself.
OUT <- {
  base <- if (identical(FAMILY, "logit")) sprintf("_dev/sbc-%s-%d", GATE, N)
          else sprintf("_dev/sbc-%s-%s-%d", FAMILY, GATE, N)
  parts <- c(if (nzchar(TAG)) TAG,
             if (SEED_BASE != 5000L) sprintf("seed%d", SEED_BASE),
             if (BANDWIDTH != 0.1) sprintf("bw%s", format(BANDWIDTH)),
             if (WARMUP != 400L) sprintf("warm%d", WARMUP),
             if (OFFSET != 0) sprintf("off%s", format(OFFSET)),
             if (L != 100L) sprintf("L%d", L))
  sprintf("%s%s.rds", base,
          if (length(parts)) paste0("-", paste(parts, collapse = "-")) else "")
}

# Resume. Every replicate seeds itself from `SEED_BASE + r` before it draws
# anything, so picking up at the next index reproduces exactly what an
# uninterrupted run would have produced; a killed run loses at most the
# replicate it was in. Delete the output file to start over.
START_AT <- 1L

if (file.exists(OUT)) {
  prev <- readRDS(OUT)

  if (!isTRUE(prev$complete) && isTRUE(prev$done >= 1L) &&
      isTRUE(prev$total == reps)) {
    ranks <- prev$res$rank
    truths <- prev$res$truth
    widths <- prev$res$width
    covered <- prev$res$covered
    START_AT <- prev$done + 1L
    cat(sprintf("resuming at replicate %d of %d\n", START_AT, reps))
  }
}

todo <- if (START_AT > reps) integer(0) else seq(START_AT, reps)

# `SBC_REPLAY_ONLY` redraws each replicate's forest and records the bandwidths
# it drew, with no fitting at all. The generating forest is a deterministic
# function of `SEED_BASE + r`, so this recovers for a finished run what it did
# not store, and the truth it also records is the check that the replay is
# faithful: it has to equal the truth that run saved.
if (nzchar(Sys.getenv("SBC_REPLAY_ONLY"))) {
  out <- data.frame(rep = seq_len(reps), truth = NA_real_,
                    bw_mean = NA_real_, bw_max = NA_real_)

  for (r in seq_len(reps)) {
    set.seed(SEED_BASE + r)
    eta <- draw_forest(u)
    out$truth[r] <- eta[A] - eta[B]
    out$bw_mean[r] <- mean(LAST_BW)
    out$bw_max[r] <- max(LAST_BW)
  }

  file <- sprintf("_dev/sbc-replay-%s-%s-%d-%d.rds", FAMILY, GATE, N, SEED_BASE)
  saveRDS(list(res = out, complete = TRUE, seed_base = SEED_BASE), file)
  cat(sprintf("replayed %d forests to %s\n", reps, file))
  quit(save = "no")
}

pr <- prog_init(total = max(length(todo), 1L),
                title = sprintf("SBC: %s, %s rules, n = %d%s", FAMILY, GATE, N,
                                if (START_AT > 1L)
                                  sprintf(" (resumed at %d)", START_AT) else ""),
                unit = "replicate", kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last replicate"),
        add = TRUE)

# Written after every replicate rather than after the last one, so a run that
# is killed leaves its finished replicates readable. `complete` is what tells a
# reader which of the two it is looking at.
checkpoint <- function(complete, done) {
  saveRDS(list(res = data.frame(rank = ranks, truth = truths, width = widths,
                                covered = covered, n = N, gate = GATE),
               complete = complete, done = done, total = reps, draws = L),
          OUT)
}

for (r in todo) {
  set.seed(SEED_BASE + r)
  t0 <- Sys.time()

  eta <- draw_forest(u)
  truth <- eta[A] - eta[B]

  d <- d0
  d$y <- draw_response(eta + OFFSET)

  # A response with no variation carries no likelihood and the fit refuses it.
  if (length(unique(d$y)) < 2L) {
    prog_tick(pr, ok = FALSE, label = sprintf("replicate %d", r),
              msg = "response had no variation")
    next
  }

  fit <- bartisan(y ~ ., data = d, family = family, control = control,
                  offset = if (OFFSET != 0) rep(OFFSET, N))

  e <- fit[["eta"]][[1L]]
  contrast <- e[, A] - e[, B]
  thin <- contrast[seq(1L, length(contrast), length.out = L)]

  ranks <- c(ranks, sum(thin < truth))
  truths <- c(truths, truth)
  ci <- stats::quantile(contrast, c(0.025, 0.975), names = FALSE)
  widths <- c(widths, ci[2L] - ci[1L])
  covered <- c(covered, as.numeric(ci[1L] <= truth && truth <= ci[2L]))

  prog_tick(pr, label = sprintf("replicate %d", r),
            secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  checkpoint(FALSE, r)
}

checkpoint(TRUE, reps)
out <- data.frame(rank = ranks, truth = truths, width = widths,
                  covered = covered, n = N, gate = GATE)

on.exit()
prog_end(pr, "done", sprintf("%d replicates", nrow(out)))

# ---- the report --------------------------------------------------------------

bins <- 10L
counts <- table(cut(out$rank, breaks = seq(0, L + 1, length.out = bins + 1L),
                    include.lowest = TRUE))
expected <- nrow(out) / bins
chisq <- sum((as.numeric(counts) - expected)^2) / expected

cat(sprintf("\n%s, n = %d, %s rules: %d replicates, %d thinned draws each\n\n",
            FAMILY, N, GATE, nrow(out), L))
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
