# Do the compiled Poisson family and a hand-written one realize one posterior?
#
# What is tested. The soft-rule Poisson arm of `_dev/sbc.R` deviates from
# uniform ranks with a dispersion effect size of about 0.09 and coverage 0.937
# to 0.942 against a nominal 0.95, and quadrupling the chain does not shrink
# it, so the sampler converges to a posterior that is slightly too narrow for
# the calibrated contrast rather than merely failing to explore one. Holding
# every tree's bandwidth fixed removes the deviation. The question this answers
# is which code the fault is in. `poisson()` and
# `custom_family(function(y, eta) dpois(y, exp(eta[, 1]), log = TRUE))`
# specify the same density, so the posterior they define is the same one; the
# paths differ in the blocked density, in the derivatives (analytic against the
# base class's central differences) and in the declared target form
# (`TARGET_EXP_UP`, which soft rules decline, against `TARGET_GENERAL`), and
# they share the bandwidth move, the Laplace tree machinery and the chunked
# likelihood difference the bandwidth move uses.
#
# What is measured. Four data sets from the `_dev/sbc.R` generator at n = 400
# with two predictors, each fit both ways under soft rules at that design -- 20
# trees, `sigma_mu` fixed at 0.35, bandwidth prior mean 0.1,
# `x_transform = "range"` -- with 400 warmup and 4000 kept draws, one chain.
# Per fit:
#
#   contrast_mean, contrast_sd   the posterior of eta[A] - eta[B] between the
#                                two observations extreme in the first
#                                predictor, which is what SBC calibrates
#   bw_mean, bw_sd               the per-tree bandwidth pooled over trees
#   loglik_mean                  the reported log likelihood
#   ess_contrast, ess_bw         what each carries, so that a difference can be
#                                read against its own Monte Carlo error
#
# What each outcome would mean. Agreement within Monte Carlo error says both
# paths realize the same posterior, which exonerates the compiled family's own
# code and leaves the machinery the two share, the bandwidth move foremost. A
# contrast or a bandwidth posterior that differs by more than its Monte Carlo
# error localizes the fault to whichever path differs, and the derivatives are
# the first place to read, since the hand-written family differences them
# numerically where the compiled one supplies them in closed form.
#
# Run with: Rscript _dev/poisson-paths.R [datasets]
# Writes:   _dev/poisson-paths.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

args <- commandArgs(trailingOnly = TRUE)
SETS <- if (length(args) > 0L) as.integer(args[1L]) else 4L
OUT <- "_dev/poisson-paths.rds"

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
A_i <- which.min(u[, 1L])
B_i <- which.max(u[, 1L])

control <- bartisan_control(num_trees = TREES, num_burn = 400L,
                            num_draws = 4000L, chains = 1L, gate = "smoothstep",
                            sigma_mu = SIGMA_MU, update_sigma_mu = FALSE,
                            sparsity = FALSE, x_transform = "range",
                            bandwidth = BANDWIDTH, augment = FALSE,
                            verbose = FALSE)

paths <- list(compiled = function() poisson(),
              custom = function() {
                custom_family(function(y, eta) {
                  stats::dpois(y, exp(eta[, 1]), log = TRUE)
                })
              })

grid <- expand.grid(path = names(paths), set = seq_len(SETS),
                    stringsAsFactors = FALSE)

rows <- list()
START_AT <- 1L

if (file.exists(OUT)) {
  prev <- readRDS(OUT)
  if (!isTRUE(prev$complete) && !is.null(prev$rows) && nrow(prev$rows) > 0L) {
    rows <- lapply(seq_len(nrow(prev$rows)), function(i) prev$rows[i, , drop = FALSE])
    START_AT <- nrow(prev$rows) + 1L
    cat(sprintf("resuming at fit %d of %d\n", START_AT, nrow(grid)))
  }
}

todo <- if (START_AT > nrow(grid)) integer(0) else seq(START_AT, nrow(grid))

pr <- prog_init(total = max(length(todo), 1L),
                title = "Poisson: compiled against hand-written",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

ess1 <- function(x) posterior::ess_bulk(matrix(x, ncol = 1L))

for (k in todo) {
  g <- grid[k, ]

  set.seed(5000L + g$set)
  eta <- draw_forest(u)
  d <- d0
  d$y <- rpois(N, exp(eta))

  t0 <- Sys.time()
  set.seed(77L + g$set)
  fit <- bartisan(y ~ ., data = d, family = paths[[g$path]](),
                  control = control)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  e <- fit[["eta"]][[1L]]
  contrast <- e[, A_i] - e[, B_i]
  bw <- as.vector(fit[["bandwidth"]])

  rows[[k]] <- data.frame(
    path = g$path, set = g$set, truth = eta[A_i] - eta[B_i],
    contrast_mean = mean(contrast), contrast_sd = sd(contrast),
    ess_contrast = ess1(contrast),
    bw_mean = mean(bw), bw_sd = sd(bw),
    ess_bw = median(apply(fit[["bandwidth"]], 2L, ess1)),
    loglik_mean = mean(fit[["loglik"]]), secs = secs)

  prog_tick(pr, label = sprintf("%s / set %d", g$path, g$set), secs = secs)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)), OUT)
}

if (!length(todo)) {
  saveRDS(list(rows = do.call(rbind, rows), complete = TRUE), OUT)
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
w <- reshape(res[, c("set", "path", "contrast_mean", "contrast_sd", "bw_mean",
                     "bw_sd", "loglik_mean")],
             idvar = "set", timevar = "path", direction = "wide")
print(w, digits = 4, row.names = FALSE)
cat("\nmean difference, custom minus compiled, with the Monte Carlo error of each fit's mean:\n")
for (col in c("contrast_mean", "contrast_sd", "bw_mean", "bw_sd", "loglik_mean")) {
  a <- res[res$path == "custom", col]
  b <- res[res$path == "compiled", col]
  cat(sprintf("  %-14s %+.4f  (sd of the four differences %.4f)\n", col,
              mean(a - b), sd(a - b)))
}
mc <- res$contrast_sd / sqrt(res$ess_contrast)
cat(sprintf("\nMonte Carlo error of one fit's contrast mean: %.4f to %.4f\n",
            min(mc), max(mc)))
