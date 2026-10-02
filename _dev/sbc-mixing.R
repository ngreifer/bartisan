# Is the soft-rule Poisson SBC failure mixing, or a wrong target?
#
# What is tested. `_dev/sbc.R` at 600 replicates is uniform for hard-rule
# Poisson and for soft-rule logit, and not for soft-rule Poisson: chi-square
# 28.4 on 9 df (p = 0.001), with a quadratic contrast of +2.8 standard errors
# and coverage 0.930 against a nominal 0.95. A U-shaped rank histogram with
# below-nominal coverage is what a chain that has not explored the posterior
# produces, since draws thinned from an autocorrelated chain underrepresent the
# tails; it is also what a posterior that is genuinely too narrow produces. The
# claim here is the first: that the soft-rule Poisson chain mixes far worse on
# this design than the arms that pass, so the failure is draws rather than a
# wrong target. Both soft arms share the general Laplace leaf path, since
# `exponential_usable()` declines under soft rules and the binomial logit with
# its augmentation off has no shape shortcut either, so the target code is
# common to the arm that fails and an arm that passes.
#
# What is measured. The SBC design of `_dev/sbc.R` exactly -- 20 trees, n = 400,
# two predictors, `sigma_mu` fixed at 0.35, 400 warmup and 1000 kept draws --
# for four arms, 25 replicates each. Per fit, over the 1000 retained draws of
# the contrast eta[A] - eta[B] between the two observations extreme in the first
# predictor, which is the quantity SBC calibrates:
#
#   ess        bulk effective sample size (posterior::ess_bulk)
#   acf1       lag-1 autocorrelation
#   sd         posterior standard deviation of the contrast
#
# Read the median ESS per arm.
#
# What each outcome would mean. If soft-rule Poisson's ESS is a small fraction
# of the other arms', the failure is mixing: the remedy is more draws, the
# target is not implicated, and a longer-chain SBC run should restore
# uniformity. If all four arms have comparable ESS, mixing is excluded and the
# soft-rule Poisson target is wrong, which is a defect in the Laplace
# proposal's soft-rule information term (`sum w_i^2 info_i`) and is where to
# read next.
#
# Run with: Rscript _dev/sbc-mixing.R [replicates]
# Writes:   _dev/sbc-mixing.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) > 0L) as.integer(args[1L]) else 25L

# The generator of `_dev/sbc.R`, reproduced so that the fits see the same kind
# of data; only the forest draw matters here, not its exactness as a prior.
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

left_prob <- function(x, val, bandwidth, soft) {
  if (!soft) return(as.numeric(x <= val))
  t <- 0.5 + 0.5 * (val - x) / (bandwidth * HALF_WIDTH)
  t <- pmin(pmax(t, 0), 1)
  t * t * (3 - 2 * t)
}

eval_tree <- function(node, u, bandwidth, soft, w = rep.int(1, nrow(u))) {
  if (isTRUE(node$leaf)) return(node$mu * w)
  p <- left_prob(u[, node$var], node$val, bandwidth, soft)
  eval_tree(node$left, u, bandwidth, soft, w * p) +
    eval_tree(node$right, u, bandwidth, soft, w * (1 - p))
}

draw_forest <- function(u, soft) {
  Reduce(`+`, lapply(seq_len(TREES), function(i) {
    b <- if (soft) rexp(1L, rate = 1 / BANDWIDTH) else 0
    eval_tree(draw_tree(), u, b, soft)
  }), numeric(nrow(u)))
}

set.seed(20260905)
u <- matrix(runif(N * P), N, P)
colnames(u) <- paste0("x", seq_len(P))
d0 <- as.data.frame(u)
A_i <- which.min(u[, 1L])
B_i <- which.max(u[, 1L])

arms <- list(
  list(name = "poisson hard", family = poisson(), gate = "hard",
       draw = function(eta) rpois(N, exp(eta))),
  list(name = "poisson soft", family = poisson(), gate = "smoothstep",
       draw = function(eta) rpois(N, exp(eta))),
  list(name = "logit hard", family = binomial(), gate = "hard",
       draw = function(eta) rbinom(N, 1L, plogis(eta))),
  list(name = "logit soft", family = binomial(), gate = "smoothstep",
       draw = function(eta) rbinom(N, 1L, plogis(eta)))
)

grid <- expand.grid(arm = seq_along(arms), rep = seq_len(REPS))

pr <- prog_init(total = nrow(grid), title = "SBC arms: mixing of the contrast",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()

for (k in seq_len(nrow(grid))) {
  arm <- arms[[grid$arm[k]]]
  r <- grid$rep[k]
  soft <- !identical(arm$gate, "hard")

  set.seed(5000L + r)
  eta <- draw_forest(u, soft)
  d <- d0
  d$y <- arm$draw(eta)

  if (length(unique(d$y)) < 2L) {
    prog_tick(pr, ok = FALSE, label = sprintf("%s rep %d", arm$name, r),
              msg = "response had no variation")
    next
  }

  control <- bartisan_control(num_trees = TREES, num_burn = 400L,
                              num_draws = 1000L, chains = 1L, gate = arm$gate,
                              sigma_mu = SIGMA_MU, update_sigma_mu = FALSE,
                              sparsity = FALSE, x_transform = "range",
                              bandwidth = BANDWIDTH, augment = FALSE,
                              verbose = FALSE)
  t0 <- Sys.time()
  fit <- bartisan(y ~ ., data = d, family = arm$family, control = control)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  e <- fit[["eta"]][[1L]]
  contrast <- e[, A_i] - e[, B_i]

  rows[[length(rows) + 1L]] <- data.frame(
    arm = arm$name, rep = r,
    ess = posterior::ess_bulk(matrix(contrast, ncol = 1L)),
    acf1 = stats::acf(contrast, lag.max = 1L, plot = FALSE)$acf[2L],
    sd = sd(contrast), truth = eta[A_i] - eta[B_i], secs = secs)

  prog_tick(pr, label = sprintf("%s rep %d", arm$name, r), secs = secs)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/sbc-mixing.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- do.call(rbind, lapply(split(res, res$arm), function(x) {
  data.frame(arm = x$arm[1L], fits = nrow(x),
             ess_median = median(x$ess), ess_min = min(x$ess),
             acf1_median = median(x$acf1), sd_median = median(x$sd),
             secs = mean(x$secs))
}))
print(out, digits = 3, row.names = FALSE)
