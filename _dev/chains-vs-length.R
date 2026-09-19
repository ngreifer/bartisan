# Many short chains against one long chain, scored on interval coverage.
#
# This is the part of Krantsevich, He and Hahn (2023) that survives the transient
# measurement in `_dev/TASKS.md`. Their warm-start BCF beats BCF on coverage
# (0.96 and 0.92 against 0.90 and 0.73 for ATE and CATE at n = 5000), and their
# own long run -- 20000 draws after 20000 burn-in -- still does not match it. So
# the gain is not chain length, and the transient measurement rules out burn-in
# length. What is left in their protocol is that they run s - b *separately
# initialized* chains in parallel and pool them, against one chain for BCF.
#
# If between-chain diversity is what widens their intervals, then splitting a
# fixed sampling budget across more chains should improve coverage here too, and
# `chains` already does that with no grow-from-root anywhere. If it does not,
# their gain comes from something else and there is nothing here to take.
#
# The low-noise reading is the posterior standard deviation, not coverage.
# Coverage is a proportion averaged over replicates and needs many of them to
# resolve a point or two; the posterior spread is a mean of 500 numbers per fit
# and moves immediately. If pooling 16 chains does not widen the posterior
# relative to one chain, the mechanism is absent and coverage cannot improve by
# this route whatever the coverage column happens to say.
#
# Arms. The first four hold total sweeps at 4800, so they cost the same and
# differ only in how the budget is cut up; `pooled` is what survives warmup and
# falls as chains rise, which is the real cost of many chains. Arm E is the
# interesting one: it matches B's pooled draws from 16 separate starts by
# spending less on each warmup, which the transient measurement licenses, since
# 50 sweeps is past the knee at this sample size. The last two ignore cost
# entirely and ask whether chain count does anything at all.
#
#   A   1 chain   200 warmup  4600 draws   4800 sweeps   4600 pooled
#   B   4 chains  200 warmup  1000 draws   4800 sweeps   4000 pooled
#   D  16 chains  200 warmup   100 draws   4800 sweeps   1600 pooled
#   E  16 chains   50 warmup   250 draws   4800 sweeps   4000 pooled
#   F   1 chain   200 warmup   400 draws    600 sweeps    400 pooled
#   G  16 chains  200 warmup   400 draws   9600 sweeps   6400 pooled
#
# Scored against the true regression function at 500 held-out points, nominal
# 95%, equal-tailed. The diagnostic columns are the ones in the coverage entry of
# `_dev/TASKS.md`, so the two tables can be read together.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

SMOKE <- nzchar(Sys.getenv("CHAINS_SMOKE"))

N_TRAIN <- 4000L
N_TEST <- 500L
P <- 30L
LEVEL <- 0.95

ARMS <- data.frame(
  arm    = c("A", "B", "D", "E", "F", "G"),
  chains = c(1L, 4L, 16L, 16L, 1L, 16L),
  burn   = c(200L, 200L, 200L, 50L, 200L, 200L),
  draws  = c(4600L, 1000L, 100L, 250L, 400L, 400L),
  stringsAsFactors = FALSE)

REPS <- 20L
OUT <- "_dev/chains-vs-length.rds"

if (SMOKE) {
  N_TRAIN <- 500L
  REPS <- 1L
  ARMS$draws <- pmax(ARMS$draws %/% 40L, 10L)
  ARMS$burn <- pmax(ARMS$burn %/% 10L, 5L)
  OUT <- "_dev/chains-vs-length-smoke.rds"
}

# Friedman, with `P - 5` irrelevant predictors after the five it uses. The
# training and test sets come from one draw so the test points lie in the same
# region, and the truth is kept separately from the noisy response.
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

run <- function(arm, rep) {
  d <- sim(rep)
  spec <- ARMS[ARMS$arm == arm, ]

  elapsed <- system.time(
    fit <- bartisan(y ~ ., data = d$train, family = gaussian(),
                    chains = spec$chains, num_burn = spec$burn,
                    num_draws = spec$draws, num_thin = 1L)
  )[["elapsed"]]

  # Draws by test points of the regression function. Gaussian with an identity
  # link, so the additive predictor is the thing being covered.
  eta <- predict(fit, newdata = d$test, type = "link", draws = TRUE)

  lo <- apply(eta, 2L, stats::quantile, probs = (1 - LEVEL) / 2)
  hi <- apply(eta, 2L, stats::quantile, probs = 1 - (1 - LEVEL) / 2)
  mid <- colMeans(eta)
  psd <- apply(eta, 2L, stats::sd)

  # R-hat and effective sample size, from `diagnose()` rather than from the fit,
  # which carries no `rhat` element. Three quantities: the log likelihood, the
  # splitting proportions, which the diagnostics entry records as the slowest
  # thing in the model and which this design stresses with 25 irrelevant
  # predictors, and the worst 5% of the predictor. R-hat is split-R-hat, so it
  # is defined for the single-chain arms too.
  pick <- function(tb, what, col) {
    row <- match(what, tb[["quantity"]])
    if (is.na(row)) NA_real_ else as.numeric(tb[[col]][row])
  }

  diag <- tryCatch({
    tb <- diagnose(fit)[["table"]]
    c(pick(tb, "loglik", "rhat"), pick(tb, "loglik", "ess_bulk"),
      pick(tb, "splits.eta", "rhat"), pick(tb, "splits.eta", "ess_bulk"),
      pick(tb, "eta.eta (worst 5% of observations)", "rhat"))
  }, error = function(e) rep(NA_real_, 5L))

  data.frame(
    arm = arm, rep = rep, chains = spec$chains, burn = spec$burn,
    draws_per_chain = spec$draws, pooled = nrow(eta),
    sweeps = spec$chains * (spec$burn + spec$draws),
    coverage = mean(d$truth >= lo & d$truth <= hi),
    width = mean(hi - lo),
    abs_bias = mean(abs(mid - d$truth)),
    post_sd = mean(psd),
    ratio = mean(abs(mid - d$truth)) / mean(psd),
    rmse = sqrt(mean((mid - d$truth)^2)),
    rhat_loglik = diag[1L], ess_loglik = diag[2L],
    rhat_splits = diag[3L], ess_splits = diag[4L],
    rhat_eta_worst = diag[5L],
    seconds = elapsed)
}

cells <- expand.grid(arm = ARMS$arm, rep = seq_len(REPS),
                     stringsAsFactors = FALSE)
cells <- cells[order(cells$rep, match(cells$arm, ARMS$arm)), ]

pr <- prog_init(total = nrow(cells), title = "Chains against chain length",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()

for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]
  label <- sprintf("rep %d arm %s (%d chains)", cell$rep, cell$arm,
                   ARMS$chains[match(cell$arm, ARMS$arm)])

  rows[[length(rows) + 1L]] <- prog_do(pr, i, function(.) {
    run(cell$arm, cell$rep)
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
