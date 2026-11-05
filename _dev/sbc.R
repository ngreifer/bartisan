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
# predictor minus another's. The centering is common to both and cancels from
# the estimand. Calibrating the level instead would be measuring the centering,
# which absorbs it entirely and would tell us nothing.
#
# It does not cancel from the posterior, which an earlier version of this
# comment claimed (corrected 2026-11-04). The leaf prior does not separate the
# level from the contrast: a leaf holding A and most of the sample moves both,
# so the two covary a priori. The fit's prior is the tree prior on the forest
# around the empirical anchor, the generator's is the tree prior around zero,
# and the data pin the level, so the two give the contrast different
# conditional priors. `SBC_NO_CENTER` below removes the difference exactly.
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
#
# "binomial" (added 2026-11-04) is a binomial logit with `SBC_TRIALS` trials
# per observation (default 20), fit as `cbind(successes, failures)`.
# `SBC_AUGMENT` sets `augment` for any family ("1" on, "0" off); unset, only
# probit is augmented. Together they split the soft-gate Poisson narrowing:
#
# What is tested. Whether the narrowing comes from the general leaf path (a
# Laplace proposal with a Metropolis correction) meeting a skewed likelihood
# while the bandwidth moves, or from something specific to `poisson()`. A
# binomial with many trials and a small probability is shaped like a Poisson,
# and it can be fit directly on the general path or through the Polya-Gamma
# augmentation, exact at an integer number of trials, which makes each leaf
# target quadratic as *SoftBart*'s Gaussian targets are.
#
# What is measured. 20 trials, `SBC_OFFSET=-2.7` for counts at about the
# Poisson arm's level (prior median 1.3 and mean 2.5, against 1 and 3.4), `SBC_NO_CENTER` so the fit's prior is the generator's, soft
# rules with the bandwidth drawn, 600 replicates, augmentation on and off on
# the same data sets. The dispersion z of the contrast's ranks in each arm and
# their paired difference, against soft Poisson's +2.61.
#
# What each outcome would mean. Only the direct arm narrows: the general path
# under a moving bandwidth is the cause. Both narrow: the leaf path is cleared
# and the bandwidth move under a skewed likelihood is the suspect. Neither: a
# Poisson-shaped likelihood is not enough, and the cause is in `poisson()`
# itself. `_dev/TASKS.md`, 2026-11-04, has the full statement.
#
# "negbin" (added 2026-11-04) is `negbin(theta = SBC_THETA)`, default 5, with
# the same theta in the generator, so the predictor is the only parameter and
# SBC is valid. It asks whether the soft-gate narrowing reaches `negbin()`,
# augmented by default (a Poisson-shaped leaf target given the gamma rates)
# and direct (`SBC_AUGMENT=0`); the statement is in `_dev/TASKS.md`.
#
# `SBC_COLLAPSED` (added 2026-11-04) sets the internal flag
# `the$collapsed_bandwidth`, which proposes each tree's bandwidth with its
# leaves integrated out where the target is quadratic and then draws the
# leaves jointly. It tests whether holding the leaves fixed in the bandwidth
# move is what narrows the soft-rule posteriors; run it on the augmented
# binomial arm, whose target is quadratic. The statement is in
# `_dev/TASKS.md`. It needs a build that has the flag, and stops otherwise.
#
# A single tree and several chains (added 2026-11-04), from the BART mixing
# literature (Ronen et al. 2022, Tan et al. 2026, Kim and Rockova 2025): a
# chain held in one mode of the tree posterior reports only that mode's spread.
#
# What is tested. Whether the narrowing needs an ensemble (`SBC_TREES=1`), and
# whether pooling separately started chains removes it (`SBC_CHAINS=4`).
#
# What is measured. One soft tree with the bandwidth drawn and one hard tree,
# `SBC_SIGMA_MU=1.565` so the function's prior spread matches 20 trees at .35,
# soft Poisson with `SBC_NO_CENTER`, 1200 replicates each. Four chains of 250
# draws against the single-chain augmented binomial and soft Poisson runs on
# the same seeds, 200 replicates each, read paired.
#
# What each outcome would mean. Only the soft single tree narrows: it lives in
# one soft tree. Both single trees narrow: tree-structure mixing, which more
# trees relieve. Neither: it needs several trees. Four chains remove it: mode
# confinement is the cause and more chains or better structure moves the fix;
# they do not: what remains is the Stan reference. `_dev/TASKS.md` has the
# full statement.
#
# A single tree with four chains (added 2026-11-04, after both single-tree arms
# narrowed). Tested: whether that narrowing is a chain confined to one mode of
# the tree posterior. Measured: `SBC_TREES=1` with `SBC_CHAINS=4`, both gates,
# 1200 replicates, paired against the single-chain arms. Meaning: if it goes,
# the single tree is a mixing problem separate from the 20-tree soft one; if it
# stays, one tree is the smallest case of a too-narrow posterior and the place
# to build the Stan reference.
#
# `SBC_INIT_PRIOR` (added 2026-11-04) sets the internal flag
# `the$init_from_prior`, which starts each chain's trees from a draw of the
# tree prior rather than from stumps. Tested: whether the 20-tree soft
# narrowing is a region all chains reach from the common stump start.
# Measured: four chains from prior draws against four from stumps, same 200
# seeds, soft Poisson and the augmented binomial. Meaning: if it shrinks, the
# start is the cause; if not, the Stan reference is what is left.
#
# `SBC_NO_CHANGE` (added 2026-11-04) sets the internal flag
# `the$no_change_move`, which drops the change move and leaves birth and death,
# a valid sampler on their own. Tested: whether the change move is what leaves
# a too-narrow distribution invariant. Measured: one chain, soft Poisson and the
# augmented binomial, 200 replicates paired against the single-chain runs.
# Meaning: if the narrowing goes, the change move is the defect; if not, birth
# and death or the problem remain.
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
# the same constant to the fit, which raises the level of the response. For a
# Poisson it is how far the log density is pushed toward a quadratic: at a mean
# near one it is sharply skewed and at a mean of fifty it is nearly normal, so
# the Laplace approximation the tree moves accept on gets better as this grows.
#
# It does not leave the model unchanged, which this comment said until
# 2026-11-04. The fit's anchor is `log(mean(y))`, which ignores the offset, and
# the offset is added on top of it, so the anchor sits near twice the offset and
# the forest has to carry about minus the offset by itself, against a prior on
# its level with standard deviation sqrt(20) * 0.35, about 1.6. The offset runs
# of 2026-10-02 therefore raised the counts and the prior mismatch together.
OFFSET <- as.numeric(Sys.getenv("SBC_OFFSET", "0"))
# `SBC_NO_CENTER` makes the fit's predictor `OFFSET + forest` exactly, which is
# the generator's. The fit's prior is then the generating prior with no
# empirical part, and the level becomes a valid target too; its rank is recorded
# as `rank_level` (the predictor at A).
#
# How depends on the installed build. Since 2026-11-04 the package computes the
# intercept-only fit with the offset included, so a constant offset is absorbed
# and can no longer cancel it; there, the fit is given `anchor = 0` in
# `bartisan_control()` and the offset `OFFSET`. The CRAN 0.1.0 build has no
# `anchor` and ignores the offset when it computes the start, so there the fit
# is given `OFFSET - log(mean(y))` as its offset, which the start adds back to.
# The four runs of 2026-11-04 used the second route.
#
# What is tested (2026-11-04). The claim, raised in an outside comment, that
# the empirical anchor is what narrows the soft-rule Poisson posterior: that
# the dispersion contrast near +2 with the bandwidth drawn comes from the
# mismatch between the fit's prior and the generator's rather than from the
# sampler. Against it already: the anchor is identical in the arms that pass
# (bandwidth fixed, hard rules, soft logit). Separately, the claim that the
# excess at offset 4 and 1000 draws (+4.3, `_dev/TASKS.md` 2026-10-02) came
# from the anchor sitting near twice the offset rather than from the counts.
#
# What is measured. Four arms of 600 replicates on the default seed block,
# soft rules, bandwidth drawn, 1000 draws, so replicate r is the same data set
# and the same truth in every arm: offsets 0 and 4, each with the anchor as
# the package sets it and with `SBC_NO_CENTER`. The pairs differ only in the
# anchor. Read the dispersion z of the contrast's ranks (positive is too
# narrow; the default arm read +2.78), with the mean-rank z beside it and
# coverage as the noisier proxy. In the cancelled arms, read `rank_level` too.
#
# What each outcome would mean. At offset 0, a dispersion z that stays near +2
# with the anchor cancelled clears the anchor, and the bandwidth-drawn reading
# of the record stands; one that falls to about zero, the fixed-bandwidth
# arm's -0.24, means the anchor was the cause and every soft-Poisson
# conclusion since 2026-10-01 has to be re-read. At offset 4, the uncancelled
# arm should reproduce about +4.3. If the cancelled arm comes back to the
# offset-0 level, the growth with the counts was the anchor and the Laplace
# reading of it is withdrawn for good; if it stays near +4.3, the growth is
# real and the counts are what track it.
NO_CENTER <- nzchar(Sys.getenv("SBC_NO_CENTER"))
FIX_BANDWIDTH <- nzchar(Sys.getenv("SBC_FIX_BANDWIDTH"))
TAG <- Sys.getenv("SBC_TAG", "")
TRIALS <- as.integer(Sys.getenv("SBC_TRIALS", "20"))
THETA <- as.numeric(Sys.getenv("SBC_THETA", "5"))
COLLAPSED <- nzchar(Sys.getenv("SBC_COLLAPSED"))
INIT_PRIOR <- nzchar(Sys.getenv("SBC_INIT_PRIOR"))
NO_CHANGE <- nzchar(Sys.getenv("SBC_NO_CHANGE"))
AUGMENT <- {
  if (nzchar(Sys.getenv("SBC_AUGMENT"))) Sys.getenv("SBC_AUGMENT") == "1"
  else identical(FAMILY, "probit")
}

family <- switch(FAMILY,
                 logit = stats::binomial(),
                 probit = stats::binomial("probit"),
                 poisson = stats::poisson(),
                 binomial = stats::binomial(),
                 negbin = negbin(theta = THETA),
                 stop("family must be logit, probit, poisson, binomial or negbin"))

draw_response <- switch(FAMILY,
                        logit = function(eta) stats::rbinom(N, 1L, stats::plogis(eta)),
                        probit = function(eta) stats::rbinom(N, 1L, stats::pnorm(eta)),
                        poisson = function(eta) stats::rpois(N, exp(eta)),
                        binomial = function(eta) {
                          stats::rbinom(N, TRIALS, stats::plogis(eta))
                        },
                        negbin = function(eta) {
                          stats::rnbinom(N, size = THETA, mu = exp(eta))
                        })

FORMULA <- {
  if (identical(FAMILY, "binomial")) cbind(y, fail) ~ x1 + x2
  else y ~ x1 + x2
}

P <- 2L
# `SBC_TREES` and `SBC_SIGMA_MU` (added 2026-11-04) set the number of trees and
# the leaf scale, in the generator and the fit alike, so a single tree (Bayesian
# CART) can be calibrated with the same harness. `SBC_CHAINS` fits each
# replicate with that many chains and ranks the truth among their pooled draws;
# `SBC_DRAWS` is then the draws per chain.
TREES <- as.integer(Sys.getenv("SBC_TREES", "20"))
SIGMA_MU <- as.numeric(Sys.getenv("SBC_SIGMA_MU", "0.35"))
CHAINS <- as.integer(Sys.getenv("SBC_CHAINS", "1"))

# `SBC_WORKERS` (added 2026-11-04) sets how many processes fit replicates at
# once. Unset, `prog_parallel()` asks `prog_cores()`, which reads the cores and
# the load. One run either way: one file, one widget run, one report.
WORKERS <- {
  w <- Sys.getenv("SBC_WORKERS")
  if (nzchar(w)) as.integer(w) else NULL
}
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

HAS_ANCHOR <- "anchor" %in% names(formals(bartisan_control))

if (COLLAPSED) {
  if (!"collapsed_bandwidth" %in% names(bartisan_control())) {
    stop("`SBC_COLLAPSED` needs a build with `the$collapsed_bandwidth`")
  }

  flags <- bartisan:::the
  flags$collapsed_bandwidth <- TRUE
}

if (INIT_PRIOR) {
  if (!"init_from_prior" %in% names(bartisan_control())) {
    stop("`SBC_INIT_PRIOR` needs a build with `the$init_from_prior`")
  }

  flags <- bartisan:::the
  flags$init_from_prior <- TRUE
}

if (NO_CHANGE) {
  if (!"no_change_move" %in% names(bartisan_control())) {
    stop("`SBC_NO_CHANGE` needs a build with `the$no_change_move`")
  }

  flags <- bartisan:::the
  flags$no_change_move <- TRUE
}

control <- do.call(bartisan_control, c(
  list(num_trees = TREES, num_burn = WARMUP, num_draws = DRAWS, chains = CHAINS,
       gate = GATE, sigma_mu = SIGMA_MU, update_sigma_mu = FALSE,
       sparsity = FALSE, x_transform = "range", bandwidth = BANDWIDTH,
       update_bandwidth = !FIX_BANDWIDTH,
       augment = AUGMENT),
  if (NO_CENTER && HAS_ANCHOR) list(anchor = 0)))

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
             if (DRAWS != 1000L) sprintf("draws%d", DRAWS),
             if (OFFSET != 0) sprintf("off%s", format(OFFSET)),
             if (NO_CENTER) "nocenter",
             if (identical(FAMILY, "binomial")) sprintf("trials%d", TRIALS),
             if (identical(FAMILY, "negbin")) sprintf("theta%s", format(THETA)),
             if (COLLAPSED) "collapsed",
             if (nzchar(Sys.getenv("SBC_AUGMENT")))
               if (AUGMENT) "aug" else "direct",
             if (TREES != 20L) sprintf("trees%d", TREES),
             if (SIGMA_MU != 0.35) sprintf("smu%s", format(SIGMA_MU)),
             if (CHAINS != 1L) sprintf("chains%d", CHAINS),
             if (INIT_PRIOR) "initprior",
             if (NO_CHANGE) "nochange",
             if (L != 100L) sprintf("L%d", L))
  sprintf("%s%s.rds", base,
          if (length(parts)) paste0("-", paste(parts, collapse = "-")) else "")
}

LABEL <- sprintf("%s, n = %d, %s rules%s%s", FAMILY, N, GATE,
                 if (TREES != 20L) sprintf(", %d tree%s", TREES,
                                           if (TREES == 1L) "" else "s") else "",
                 if (CHAINS != 1L) sprintf(", %d chains", CHAINS) else "")

todo <- seq_len(reps)

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

# The rank of `value` among `draws`, with ties broken uniformly at random. A
# quantity with an atom needs it: with a single tree the contrast is exactly
# zero whenever no split separates A from B, and `sum(draws < value)` then puts
# every such replicate at the bottom, which read as a dispersion z of +21 for
# one hard tree (2026-11-04). With no ties it is the old rank exactly.
tied_rank <- function(draws, value) {
  sum(draws < value) + sample.int(sum(draws == value) + 1L, 1L) - 1L
}

# One replicate, start to finish, so that the replicates can be handed to
# worker processes. Seeded here, from `SEED_BASE + r` with R's default
# generators named explicitly, rather than by `prog_parallel()`: replicate r
# then draws the same truth and the same data in every run this script has
# made, sequential or parallel, which is what pairs one run against another.
# A worker starts under `future`'s L'Ecuyer generator, so the kinds have to be
# named for that to hold.
one_replicate <- function(r) {
  set.seed(SEED_BASE + r, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")

  eta <- draw_forest(u)
  truth <- eta[A] - eta[B]

  d <- d0
  d$y <- draw_response(eta + OFFSET)
  d$fail <- TRIALS - d$y

  # A response with no variation carries no likelihood and the fit refuses it.
  # Recorded as a failed unit, as the sequential loop recorded it.
  if (length(unique(d$y)) < 2L) {
    stop("response had no variation")
  }

  # The anchor the package will add, computed as `prepare_response()` does, so
  # that subtracting it here leaves the fit's predictor at `OFFSET + forest`.
  anchor <- log(max(mean(d$y), 1 / (2 * N)))

  if (NO_CENTER && !HAS_ANCHOR && !identical(FAMILY, "poisson")) {
    stop("the offset route to `SBC_NO_CENTER` is written for poisson only; ",
         "install a build with `anchor`")
  }

  fit_offset <- {
    if (NO_CENTER && !HAS_ANCHOR) rep(OFFSET - anchor, N)
    else if (OFFSET != 0) rep(OFFSET, N)
  }

  # Through `do.call()` so the offset reaches the model frame as a value: the
  # model frame looks for an `offset` expression in the formula's environment,
  # which is the global one, and `fit_offset` lives here.
  fit <- do.call(bartisan, list(FORMULA, data = d, family = family,
                                control = control, offset = fit_offset))

  # Either route has to leave the start at the generator's: zero before the
  # offset under `anchor`, the empirical anchor under the offset route.
  if (NO_CENTER) {
    stopifnot(abs(fit[["intercept"]] - if (HAS_ANCHOR) 0 else anchor) < 1e-12)
  }

  e <- fit[["eta"]][[1L]]
  contrast <- e[, A] - e[, B]
  keep <- seq(1L, length(contrast), length.out = L)
  ci <- stats::quantile(contrast, c(0.025, 0.975), names = FALSE)

  data.frame(rep = r,
             rank = tied_rank(contrast[keep], truth),
             rank_level = tied_rank(e[keep, A], eta[A] + OFFSET),
             truth = truth,
             width = ci[2L] - ci[1L],
             covered = as.numeric(ci[1L] <= truth && truth <= ci[2L]))
}

pr <- prog_init(total = length(todo), title = sprintf("SBC: %s", LABEL),
                unit = "replicate", kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last replicate"),
        add = TRUE)

# Each worker checkpoints the replicates it finishes beside the output file,
# and rerunning the same command resumes from them, so a killed run loses at
# most the replicates in flight. `prog_collect()` on the stem reads a partial
# run back.
rows <- prog_parallel(pr, todo, one_replicate, workers = WORKERS,
                      label = function(r) sprintf("replicate %d", r),
                      seed = FALSE, checkpoint = sub("\\.rds$", "", OUT),
                      packages = "bartisan")

out <- do.call(rbind, rows)
out$n <- N
out$gate <- GATE

saveRDS(list(res = out, complete = TRUE, done = reps, total = reps, draws = L,
             label = LABEL, level_valid = NO_CENTER),
        OUT)

on.exit()
prog_end(pr, "done", sprintf("%d replicates", nrow(out)))

# ---- the report --------------------------------------------------------------

# In `_dev/sbc-report.R`, shared with anything that rereads a run's file.
source("_dev/sbc-report.R")
sbc_report(out, L, LABEL, level = NO_CENTER)
