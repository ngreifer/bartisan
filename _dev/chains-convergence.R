# Is the chain-count effect on coverage real, or a pattern in unconverged fits?
#
# WHAT IS BEING TESTED
#
# `_dev/chains-vs-length.R` found that splitting a fixed sweep budget across more
# chains widens the posterior and raises coverage: one chain of 4600 draws gives
# 0.949 coverage at width 0.538, four chains of 1000 give 0.975 at 0.636, at
# identical compute. Two readings, and they imply opposite things:
#
#   Artifact. Neither arm is converged -- split-R-hat was 1.13 on the log
#   likelihood and 1.27 on the splitting proportions, against `diagnose()`'s
#   threshold of 1.01 -- so the two arms are differently-wrong rather than
#   differently-informative, and two estimates of the same posterior disagreeing
#   is just what non-convergence looks like. If so, the gap must close as the
#   chains are run longer, and the whole effect is a reason to run longer rather
#   than a reason to use more chains.
#
#   Real. A single chain settles into one variable-selection state and reports
#   the spread within it, while several chains find several and report the spread
#   across them. If so the gap persists at any run length, because it is
#   multimodality in the sparsity prior rather than Monte Carlo error, and the
#   single-chain interval is understating a real posterior.
#
# The second reading names a mechanism, so it is testable directly and not only
# by running longer: the mode in question is the set of predictors the forest
# splits on. `sparsity = FALSE` removes the Dirichlet prior over that set. If the
# gap is sparsity multimodality it should shrink or vanish with sparsity off,
# whatever the run length. If it survives with sparsity off, it is not that.
#
# WHAT THE OUTCOME MEASUREMENT IS
#
# Friedman at n = 4000 with 25 of 30 predictors irrelevant. Coverage of the true
# regression function at 500 held-out points, nominal 95%, and the mean interval
# width. Four replicates.
#
# **The decisive column is the one-chain-minus-four-chain gap in interval width,
# not coverage.** Width is a mean over 500 points and moves on four replicates;
# coverage is a proportion and would need twenty. The question is whether that
# gap goes to zero, and width answers it with far less noise.
#
# The ladder varies total sweeps at a *constant* number of stored draws, by
# thinning: `num_draws = 800` throughout with `num_thin` of 1, 4, 16 and 64, so
# the arms run 800, 3200, 12800 and 51200 sweeps past warmup while every arm
# reports 800 draws per chain. That holds Monte Carlo error from the draw count
# fixed across the ladder, so a change in width is a change in the posterior
# being sampled rather than in how finely it was sampled -- and it keeps memory
# bounded, since 4 chains of 51200 retained draws at n = 4000 would be gigabytes
# of `eta`.
#
# WHAT EACH OUTCOME WOULD MEAN
#
#   Gap closes as sweeps rise, with sparsity on. The effect was non-convergence.
#   The recommendation becomes a draws default, not a chains default, and the
#   previous entry's coverage table is an artifact to be labeled as one.
#
#   Gap persists with sparsity on but closes with sparsity off. The effect is
#   multimodality in the variable-selection state. Then the single-chain interval
#   really is too narrow on sparse problems, several chains is the honest
#   default there, and this is documentation rather than a code change.
#
#   Gap persists with sparsity off too. Neither explanation holds and something
#   else is generating it; the next thing to look at would be the leaf scale,
#   which the propensity entry in `_dev/TASKS.md` already records as failing to
#   settle on problems the likelihood likes.
#
#   R-hat never reaches 1.01 anywhere on the ladder. That is a finding in its own
#   right and bears on the defaults rather than on chains: it would mean this
#   model cannot be run to the package's own convergence standard at this sample
#   size in any practical time, and `?bartisan_control` should say which rows of
#   `diagnose()` bind in practice instead of implying all of them can pass.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

SMOKE <- nzchar(Sys.getenv("CONV_SMOKE"))

N_TRAIN <- 4000L
N_TEST <- 500L
P <- 30L
LEVEL <- 0.95
DRAWS <- 800L
BURN <- 200L

# thin 1, 4, 16 at both chain counts and both sparsity settings, plus one deeper
# rung at a single chain to see whether R-hat ever passes.
ARMS <- rbind(
  expand.grid(thin = c(1L, 4L, 16L), chains = c(1L, 4L),
              sparsity = c(TRUE, FALSE)),
  data.frame(thin = 64L, chains = 1L, sparsity = TRUE))

REPS <- 4L
OUT <- "_dev/chains-convergence.rds"

if (SMOKE) {
  N_TRAIN <- 500L
  REPS <- 1L
  DRAWS <- 40L
  BURN <- 20L
  ARMS <- ARMS[ARMS$thin <= 4L, ]
  OUT <- "_dev/chains-convergence-smoke.rds"
}

sim <- function(rep) {
  set.seed(31000 + 131 * rep)
  n <- N_TRAIN + N_TEST
  X <- matrix(stats::runif(n * P), n, P,
              dimnames = list(NULL, paste0("x", seq_len(P))))
  truth <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] + 5 * X[, 5]
  dat <- cbind(data.frame(y = truth + stats::rnorm(n)), as.data.frame(X))
  train <- seq_len(N_TRAIN)

  list(train = dat[train, ], test = dat[-train, ], truth = truth[-train])
}

pick <- function(tb, what, col) {
  row <- match(what, tb[["quantity"]])
  if (is.na(row)) NA_real_ else as.numeric(tb[[col]][row])
}

run <- function(i_arm, rep) {
  d <- sim(rep)
  spec <- ARMS[i_arm, ]

  elapsed <- system.time(
    fit <- bartisan(y ~ ., data = d$train, family = gaussian(),
                    chains = spec$chains, num_burn = BURN,
                    num_draws = DRAWS, num_thin = spec$thin,
                    sparsity = spec$sparsity)
  )[["elapsed"]]

  eta <- predict(fit, newdata = d$test, type = "link", draws = TRUE)
  lo <- apply(eta, 2L, stats::quantile, probs = (1 - LEVEL) / 2)
  hi <- apply(eta, 2L, stats::quantile, probs = 1 - (1 - LEVEL) / 2)
  mid <- colMeans(eta)

  diag <- tryCatch({
    tb <- diagnose(fit)[["table"]]
    c(pick(tb, "loglik", "rhat"), pick(tb, "loglik", "ess_bulk"),
      pick(tb, "splits.eta", "rhat"), pick(tb, "splits.eta", "ess_bulk"),
      pick(tb, "eta.eta (worst 5% of observations)", "rhat"),
      pick(tb, "eta.eta (worst 5% of observations)", "ess_bulk"))
  }, error = function(e) rep(NA_real_, 6L))

  data.frame(
    rep = rep, chains = spec$chains, thin = spec$thin,
    sparsity = spec$sparsity,
    sweeps_post_warmup = DRAWS * spec$thin,
    total_sweeps = spec$chains * (BURN + DRAWS * spec$thin),
    pooled = nrow(eta),
    coverage = mean(d$truth >= lo & d$truth <= hi),
    width = mean(hi - lo),
    abs_bias = mean(abs(mid - d$truth)),
    post_sd = mean(apply(eta, 2L, stats::sd)),
    rmse = sqrt(mean((mid - d$truth)^2)),
    rhat_loglik = diag[1L], ess_loglik = diag[2L],
    rhat_splits = diag[3L], ess_splits = diag[4L],
    rhat_eta = diag[5L], ess_eta = diag[6L],
    seconds = elapsed)
}

cells <- expand.grid(i_arm = seq_len(nrow(ARMS)), rep = seq_len(REPS))
cells <- cells[order(cells$rep, cells$i_arm), ]

pr <- prog_init(total = nrow(cells), title = "Chain count against convergence",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()

for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]
  spec <- ARMS[cell$i_arm, ]
  label <- sprintf("rep %d: %d chains, thin %d, sparsity %s", cell$rep,
                   spec$chains, spec$thin, spec$sparsity)

  rows[[length(rows) + 1L]] <- prog_do(pr, i, function(.) {
    run(cell$i_arm, cell$rep)
  }, label = label)

  # After every fit, so a killed run leaves every finished one behind.
  saveRDS(list(res = do.call(rbind, rows), arms = ARMS, complete = FALSE,
               done = length(rows), total = nrow(cells)), OUT)
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, arms = ARMS, complete = TRUE, done = nrow(res),
             total = nrow(cells)), OUT)

on.exit()
prog_end(pr, "done")

print(res, row.names = FALSE)
