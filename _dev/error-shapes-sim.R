# The measurements behind the error-shape table in vignette("families")
# § Numeric Responses: four families that make different assumptions about the
# error, over five shapes of it.
#
# Reconstructed from the vignette's own description of the design, because the
# original script was not kept -- `_dev/` is gitignored apart from an allowlist.
# The design is the one the vignette states: 1000 training and 1000 test
# observations, 50 trees, 500 draws after 500 warmup, three replicates, and every
# error centered so that E[Y | x] is the same function in every column. The
# numbers it produces therefore replace the whole table rather than one row of
# it, since the seeds cannot be recovered.
#
# Run with: Rscript _dev/error-shapes-sim.R
# Writes:   _dev/error-shapes.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

n_train <- 1000L
n_test <- 1000L
reps <- 3L
level <- 0.95

control <- bartisan_control(num_trees = 50L, num_burn = 500L,
                            num_draws = 500L, verbose = FALSE)

sim_x <- function(n) {
  data.frame(x1 = stats::runif(n), x2 = stats::runif(n), x3 = stats::runif(n))
}

mean_fun <- function(x) 1 + 0.8 * sin(pi * x$x1) + 0.6 * x$x2

# Each generator returns an error with mean zero, so E[Y | x] is `mu` in every
# column and the RMSE is scored against the same target throughout.
errors <- list(
  normal = function(n, x) stats::rnorm(n, sd = 0.3),

  `$t_3$` = function(n, x) {
    e <- stats::rt(n, df = 3)
    0.3 * e / sqrt(3)
  },

  skewed = function(n, x) {
    e <- stats::rgamma(n, shape = 1.5, rate = 1.5)
    0.3 * (e - 1)
  },

  bimodal = function(n, x) {
    side <- stats::rbinom(n, 1L, 0.5)
    0.3 * ifelse(side == 1L, 1, -1) + stats::rnorm(n, sd = 0.1)
  },

  heteroskedastic = function(n, x) {
    stats::rnorm(n, sd = 0.15 + 0.45 * x$x3)
  }
)

fits <- list(
  `gaussian()` = function(train, test) {
    fit <- bartisan(y ~ x1 + x2 + x3, data = train, family = gaussian(),
                    control = control)
    list(mean = drop(predict(fit, newdata = test, type = "response")),
         score = sum(predict(fit, newdata = test, type = "density",
                             log = TRUE)),
         fit = fit)
  },
  `dpm()` = function(train, test) {
    fit <- bartisan(y ~ x1 + x2 + x3, data = train, family = dpm(),
                    control = control)
    list(mean = drop(predict(fit, newdata = test, type = "response")),
         score = sum(predict(fit, newdata = test, type = "density",
                             log = TRUE)),
         fit = fit)
  },
  `gaussian_ls()` = function(train, test) {
    fit <- bartisan(y ~ x1 + x2 + x3, data = train, family = gaussian_ls(),
                    control = control)
    list(mean = drop(predict(fit, newdata = test, type = "response")),
         score = sum(predict(fit, newdata = test, type = "density",
                             log = TRUE)),
         fit = fit)
  },
  # Binned onto N_BINS quantiles, the same way `_dev/positive-sim.R` does it and
  # the way vignette("families") describes it: every value is replaced by the
  # mean of its bin, which leaves a numeric response whose sorted unique values
  # are the categories. `type = "mean"` then reads those values back, so the
  # prediction is already on the response's own scale and needs no mapping.
  #
  # Handing it the raw y instead makes every distinct value a category, which is
  # a thousand cutpoints at this sample size: an order of magnitude slower and a
  # different model from the one the vignette recommends.
  #
  # The cutpoints absorb the marginal distribution and the forest explains only
  # the ordering. The density is on the binned scale and so is not comparable
  # with the others'.
  `ordinal("probit")` = function(train, test) {
    edges <- unique(stats::quantile(train$y,
                                    seq(0, 1, length.out = N_BINS + 1L)))
    idx <- cut(train$y, edges, include.lowest = TRUE, labels = FALSE)

    binned <- train
    binned$y <- stats::ave(train$y, idx)

    fit <- bartisan(y ~ x1 + x2 + x3, data = binned,
                    family = ordinal("probit"), control = control)

    list(mean = drop(predict(fit, newdata = test, type = "mean")),
         score = NA_real_, fit = fit)
  }
)

OUT <- file.path("_dev", "error-shapes.rds")

N_BINS <- 25L

n_total <- length(errors) * reps * length(fits)
rows <- list()

pr <- prog_init(total = n_total, title = "Error shapes by family",
                unit = "fit", kind = "simulation", workers = 1L,
                command = "Rscript _dev/error-shapes-sim.R")

# A crash still closes the run rather than leaving it looking live.
on.exit(prog_end(pr, "failed", "aborted before the last cell"), add = TRUE)

checkpoint <- function(complete) {
  saveRDS(list(res = do.call(rbind, rows), complete = complete,
               done = length(rows), total = n_total, reps = reps,
               n_train = n_train, n_test = n_test, level = level,
               n_bins = N_BINS),
          OUT)
}

for (shape in names(errors)) {
  for (rep in seq_len(reps)) {
    set.seed(1000L * match(shape, names(errors)) + rep)

    x_train <- sim_x(n_train)
    x_test <- sim_x(n_test)
    mu_train <- mean_fun(x_train)
    mu_test <- mean_fun(x_test)

    train <- x_train
    train$y <- mu_train + errors[[shape]](n_train, x_train)
    test <- x_test
    test$y <- mu_test + errors[[shape]](n_test, x_test)

    for (family in names(fits)) {
      t0 <- proc.time()[["elapsed"]]
      got <- tryCatch(
        suppressMessages(suppressWarnings(fits[[family]](train, test))),
        error = function(e) e)
      secs <- proc.time()[["elapsed"]] - t0

      if (inherits(got, "error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          shape = shape, family = family, rep = rep, rmse = NA_real_,
          score = NA_real_, coverage = NA_real_, secs = secs)
        prog_tick(pr, i = length(rows), ok = FALSE, secs = secs,
                  label = sprintf("%s / %s rep %d", shape, family, rep),
                  msg = conditionMessage(got))
        next
      }

      # Coverage of the regression function, from the interval on the mean.
      draws <- suppressMessages(suppressWarnings(
        if (family == 'ordinal("probit")') {
          # Already on the response's scale, since the categories are the bin
          # means.
          predict(got$fit, newdata = test, type = "mean", draws = TRUE)
        } else {
          predict(got$fit, newdata = test, type = "response", draws = TRUE)
        }))
      q <- apply(draws, 2L, stats::quantile,
                 probs = c((1 - level) / 2, 1 - (1 - level) / 2))
      cover <- mean(mu_test >= q[1L, ] & mu_test <= q[2L, ])

      rows[[length(rows) + 1L]] <- data.frame(
        shape = shape, family = family, rep = rep,
        rmse = sqrt(mean((got$mean - mu_test)^2)),
        score = got$score, coverage = cover, secs = secs)

      prog_tick(pr, i = length(rows), secs = secs,
                label = sprintf("%s / %s rep %d", shape, family, rep))
    }

    # After every replicate, so a killed run keeps every finished fit.
    checkpoint(FALSE)
  }
}

checkpoint(TRUE)

on.exit()
prog_end(pr, "done", sprintf("%d fits over %d error shapes", length(rows),
                             length(errors)))
