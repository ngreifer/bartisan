# Does `_dev/sbc.R`'s R replication of the forest prior match the engine's?
#
# What is tested. The soft-rule Poisson arm of `_dev/sbc.R` fails SBC
# (chi-square 28.4 on 9 df, p = 0.001, quadratic contrast +2.8 standard errors,
# coverage 0.930) where hard-rule Poisson and soft-rule logit pass, and
# `_dev/sbc-mixing.R` excluded mixing: the failing arm has a median effective
# sample size of 249 against 101 for the hard-rule Poisson arm that passes.
# SBC fails whenever the generating prior differs from the model's, with no
# defect in the sampler at all, and the script's soft-rule prior is the part
# with real content in it: an exponential bandwidth per tree and the smoothstep
# gate, both reproduced in R. The claim: the replication is exact under hard
# rules and may not be under soft rules.
#
# What is measured. The prior distribution of the calibrated contrast,
# eta[A] - eta[B] between the two observations extreme in the first predictor,
# by two routes at identical settings (20 trees, `sigma_mu` fixed at 0.35,
# gamma 0.95, beta 2, bandwidth 0.1, `x_transform = "range"`, n = 400, two
# predictors): the script's own generator, and the engine itself through
# `prior_only = TRUE`, which ignores the response and draws from the model's
# prior. 2000 generator draws against 2000 engine draws per gate. Reported per
# gate: the standard deviation of the contrast by each route, their ratio, the
# 5th, 50th and 95th percentiles, and a two-sample Kolmogorov-Smirnov test. The
# same for the marginal predictor at a single observation, which separates a
# disagreement about the forest's overall scale from one about its smoothness.
#
# What each outcome would mean. Agreement under hard rules and disagreement
# under soft rules says the failure is in this test rather than in the sampler:
# the soft arms of `_dev/sbc.R` would then say nothing until the generator is
# fixed, and the package's soft-rule prior would need reading against the
# generator line by line. Agreement under both gates leaves the soft-rule
# Poisson target as the explanation, to be read in the Laplace proposal's
# soft-rule information term.
#
# Run with: Rscript _dev/sbc-prior-check.R
# Writes:   _dev/sbc-prior-check.rds

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
DRAWS <- 2000L

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
# A third point, near the middle, for the marginal.
M_i <- which.min(abs(u[, 1L] - 0.5))

pr <- prog_init(total = 4L, title = "SBC prior replication", unit = "block",
                kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
draws_kept <- list()

for (gate in c("hard", "smoothstep")) {
  soft <- !identical(gate, "hard")

  set.seed(99L)
  gen <- t(vapply(seq_len(DRAWS), function(i) {
    eta <- draw_forest(u, soft)
    c(eta[A_i] - eta[B_i], eta[M_i])
  }, numeric(2L)))
  prog_tick(pr, label = paste("generator", gate))

  # The engine's own prior: no response information at all, so the draws are
  # from the model's prior over forests.
  d <- d0
  d$y <- rpois(N, 1)
  set.seed(99L)
  fit <- bartisan(y ~ ., data = d, family = poisson(), prior_only = TRUE,
                  control = bartisan_control(num_trees = TREES, num_burn = 1000L,
                                             num_draws = DRAWS, chains = 1L,
                                             gate = gate, sigma_mu = SIGMA_MU,
                                             update_sigma_mu = FALSE,
                                             sparsity = FALSE,
                                             x_transform = "range",
                                             bandwidth = BANDWIDTH,
                                             augment = FALSE, verbose = FALSE))
  e <- fit[["eta"]][[1L]]
  eng <- cbind(e[, A_i] - e[, B_i], e[, M_i] - fit[["intercept"]][1L])
  prog_tick(pr, label = paste("engine", gate))

  for (j in 1:2) {
    what <- c("contrast", "one point")[j]
    ks <- suppressWarnings(stats::ks.test(gen[, j], eng[, j]))
    rows[[length(rows) + 1L]] <- data.frame(
      gate = gate, quantity = what,
      sd_gen = sd(gen[, j]), sd_eng = sd(eng[, j]),
      ratio = sd(eng[, j]) / sd(gen[, j]),
      q05_gen = quantile(gen[, j], 0.05, names = FALSE),
      q05_eng = quantile(eng[, j], 0.05, names = FALSE),
      q95_gen = quantile(gen[, j], 0.95, names = FALSE),
      q95_eng = quantile(eng[, j], 0.95, names = FALSE),
      ks_D = unname(ks$statistic), ks_p = ks$p.value)
  }
  draws_kept[[gate]] <- list(gen = gen, eng = eng)
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
saveRDS(list(rows = res, draws = draws_kept, complete = TRUE),
        "_dev/sbc-prior-check.rds")
print(res, digits = 3, row.names = FALSE)
