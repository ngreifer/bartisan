# How often does the soft-rule bandwidth actually need resampling?
#
# `bandwidth_every` is 1 by default: every tree attempts a bandwidth move on
# every sweep. That move is the most expensive thing in a soft-rule sweep -- it
# rebuilds every membership weight in the tree and evaluates the whole
# likelihood -- so raising it is the cheapest speedup on offer. The question is
# what it costs, and the answer has three parts that have to be measured
# together:
#
#   predictive   held-out RMSE of the mean function
#   inferential  coverage of 95% credible intervals for that function
#   mixing       effective sample size per second, on the worst quantity
#
# The third is the one that decides it. Fewer bandwidth moves make each sweep
# cheaper but each sweep a smaller step, so a setting can be faster per sweep
# and worse per second. Only ESS/sec settles that, and it is why this cannot be
# answered with a timing benchmark alone.
#
# `update_bandwidth = FALSE` is included as the limiting case, but it is a
# different model rather than a different sampler: the bandwidth is then fixed
# at `bandwidth` instead of learned, so any difference in fit is partly that.
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))
Sys.setenv(PROGRESS_ROOT = "/Users/NoahGreifer/Dropbox/Research/R/bartisan/.progress-runs")

suppressMessages(library(bartisan))

REPS <- 8L
N_TRAIN <- 500L
N_TEST <- 500L
CHAINS <- 4L
NUM_BURN <- 1000L
NUM_DRAWS <- 1000L
OUT <- file.path("_dev", "bandwidth-every.rds")

# Smooth: what soft rules are for, so this is where losing bandwidth moves
# should hurt most. Step: a hard boundary, where the bandwidth should want to
# be small and there is little to learn about it.
truth <- function(X, shape) {
  if (identical(shape, "smooth")) {
    3 * sin(pi * X[, 1L]) + 2 * X[, 2L] + 3 * (X[, 3L] - 0.5)^2
  }
  else {
    3 * (X[, 1L] > 0.5) + 2 * (X[, 2L] > 0.3) - 1.5 * (X[, 3L] > 0.7)
  }
}

make <- function(n, p, shape, seed) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  d <- as.data.frame(X)
  mu <- truth(X, shape)
  d$y <- mu + rnorm(n, 0, 0.7)
  list(data = d, mu = mu)
}

# The worst-mixing quantity, dropping the averages over observations: an
# average mixes about as well as the draw count whatever the sampler does, so
# keeping it would dilute the comparison.
worst_ess <- function(fit) {
  tab <- diagnose(fit)[["table"]]
  tab <- tab[!grepl("average over observations", tab[["quantity"]], fixed = TRUE), ]
  c(ess_min = min(tab[["ess_bulk"]], na.rm = TRUE),
    ess_med = stats::median(tab[["ess_bulk"]], na.rm = TRUE),
    rhat_max = max(tab[["rhat"]], na.rm = TRUE))
}

quiet <- function(expr) {
  invisible(utils::capture.output(out <- suppressMessages(suppressWarnings(expr))))
  out
}

settings <- list(
  list(label = "every = 1",   every = 1L,  update = TRUE),
  list(label = "every = 5",   every = 5L,  update = TRUE),
  list(label = "every = 10",  every = 10L, update = TRUE),
  list(label = "every = 25",  every = 25L, update = TRUE),
  list(label = "fixed",       every = 1L,  update = FALSE)
)

shapes <- c("smooth", "step")

pr <- prog_init(total = length(shapes) * length(settings) * REPS,
                title = "Bandwidth update frequency", unit = "fit",
                kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
k <- 0L

for (shape in shapes) {
  test <- make(N_TEST, 10L, shape, seed = 77L + match(shape, shapes))

  for (st in settings) {
    for (r in seq_len(REPS)) {
      train <- make(N_TRAIN, 10L, shape, seed = 5000L + 101L * r +
                      7L * match(shape, shapes))

      ctrl <- bartisan_control(num_trees = 50L, num_burn = NUM_BURN,
                               num_draws = NUM_DRAWS, chains = CHAINS,
                               bandwidth_every = st$every,
                               update_bandwidth = st$update)

      t0 <- Sys.time()
      fit <- quiet(bartisan(y ~ ., train$data, family = stats::gaussian(),
                            control = ctrl))
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      draws <- predict(fit, newdata = test$data, type = "link", draws = TRUE)
      post <- colMeans(draws)
      lo <- apply(draws, 2L, stats::quantile, 0.025)
      hi <- apply(draws, 2L, stats::quantile, 0.975)

      ess <- quiet(worst_ess(fit))
      k <- k + 1L

      rows[[k]] <- data.frame(
        shape = shape, setting = st$label, every = st$every,
        update = st$update, rep = r, secs = secs,
        rmse = sqrt(mean((post - test$mu)^2)),
        coverage = mean(lo <= test$mu & test$mu <= hi),
        width = mean(hi - lo),
        ess_min = ess[["ess_min"]], ess_med = ess[["ess_med"]],
        rhat_max = ess[["rhat_max"]],
        ess_per_sec = ess[["ess_min"]] / secs,
        bandwidth = mean(fit[["bandwidth"]]),
        stringsAsFactors = FALSE)

      prog_tick(pr, i = k, secs = secs,
                label = sprintf("%s %s rep %d", shape, st$label, r))
    }
  }
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, reps = REPS, n_train = N_TRAIN, n_test = N_TEST,
             chains = CHAINS, num_burn = NUM_BURN, num_draws = NUM_DRAWS), OUT)

on.exit()
prog_end(pr, "done")

agg <- aggregate(cbind(secs, rmse, coverage, width, ess_min, ess_per_sec,
                       rhat_max, bandwidth) ~ shape + setting,
                 data = res, FUN = mean)
print(agg[order(agg$shape, agg$setting), ], row.names = FALSE)
cat("\nwrote", OUT, "\n")
