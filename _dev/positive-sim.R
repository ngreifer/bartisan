# The measurements behind vignette("families") § Positive Continuous Responses.
#
# Every family here targets E[Y | x] on the *original* scale, so the RMSE column
# compares the same quantity and the log score column compares densities taken
# with respect to the same measure. That is the point of the design: an earlier
# version of this table compared `Gamma("log")` against `gaussian()` fitted to
# log(y), which is a different estimand (the mean of the log rather than the log
# of the mean) and a density on a different scale, so neither column meant what
# it appeared to.
#
# `gaussian("log")` is the like-for-like Gaussian comparison: the link is
# composed from R onto the identity-link family, so the forest is on the log
# mean while the error stays additive and of constant variance.
#
# `dpm()` has no link argument, so it appears only on the raw scale. Its
# additive predictor is defined to be the conditional mean, which is what makes
# the centered mixture identified, and a log link would be a different
# construction rather than a composition.
#
# Run with: Rscript _dev/positive-sim.R
# Writes:   _dev/positive-results.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

n_train <- 800L
n_test <- 800L
reps <- 4L
n_bins <- 25L

control <- bartisan_control(num_trees = 50L, num_burn = 500L,
                            num_draws = 500L)

# One mean function on the original scale, shared by every truth, so that the
# rows differ only in the shape of the error around it.
mean_fun <- function(x) {
  exp(1 + 0.8 * sin(pi * x$x1) + 0.6 * x$x2)
}

# Each generator returns a response whose conditional mean is exactly `mu`, so
# the RMSE column is scored against the same target in every row.
truths <- list(
  `gamma, constant dispersion` = function(mu, x) {
    shape <- 4
    stats::rgamma(length(mu), shape = shape, rate = shape / mu)
  },

  `gamma, varying dispersion` = function(mu, x) {
    shape <- 1 / exp(-1.5 + 1.8 * x$x3)
    stats::rgamma(length(mu), shape = shape, rate = shape / mu)
  },

  lognormal = function(mu, x) {
    sdlog <- 0.6
    mu * stats::rlnorm(length(mu), meanlog = -sdlog^2 / 2, sdlog = sdlog)
  },

  # A contaminated gamma: a tenth of the draws come from a component with a much
  # heavier right tail. Both components have mean one, so the product has mean
  # `mu` exactly rather than approximately.
  `heavy tail` = function(mu, x) {
    n <- length(mu)
    heavy <- stats::runif(n) < 0.1
    v <- numeric(n)
    v[!heavy] <- stats::rgamma(sum(!heavy), shape = 6, rate = 6)
    v[heavy] <- stats::rgamma(sum(heavy), shape = 0.5, rate = 0.5)
    mu * v
  }
)

sim_x <- function(n) {
  data.frame(x1 = stats::runif(n), x2 = stats::runif(n), x3 = stats::runif(n))
}

# `ordinal()` on a binned response: the cutpoints absorb the marginal
# distribution and `type = "mean"` reads the bin means back, so the prediction
# is on the response's own scale and comparable with the rest.
bin_response <- function(y, n_bins) {
  edges <- unique(stats::quantile(y, seq(0, 1, length.out = n_bins + 1L)))
  idx <- cut(y, edges, include.lowest = TRUE, labels = FALSE)
  stats::ave(y, idx)
}

fits <- list(
  `Gamma("log")` = function(train, test) {
    fit <- bartisan(y ~ x1 + x2 + x3, data = train, family = Gamma("log"),
                    control = control)
    list(mean = predict(fit, newdata = test, type = "response"),
         score = sum(predict(fit, newdata = test, type = "density",
                             log = TRUE)))
  },

  `Gamma_ls()` = function(train, test) {
    fit <- bartisan(y ~ x1 + x2 + x3, data = train, family = Gamma_ls(),
                    control = bartisan_control(num_trees = c(50L, 20L),
                                               num_burn = 500L,
                                               num_draws = 500L))
    list(mean = predict(fit, newdata = test, type = "response"),
         score = sum(predict(fit, newdata = test, type = "density",
                             log = TRUE)))
  },

  `gaussian("log")` = function(train, test) {
    fit <- bartisan(y ~ x1 + x2 + x3, data = train, family = gaussian("log"),
                    control = control)
    list(mean = predict(fit, newdata = test, type = "response"),
         score = sum(predict(fit, newdata = test, type = "density",
                             log = TRUE)))
  },

  `dpm()` = function(train, test) {
    fit <- bartisan(y ~ x1 + x2 + x3, data = train, family = dpm(),
                    control = control)
    list(mean = predict(fit, newdata = test, type = "response"),
         score = sum(predict(fit, newdata = test, type = "density",
                             log = TRUE)))
  },

  `ordinal("probit")` = function(train, test) {
    binned <- train
    binned$y <- bin_response(train$y, n_bins)
    fit <- bartisan(y ~ x1 + x2 + x3, data = binned,
                    family = ordinal("probit"), control = control)
    # The log score is on the binned scale, so it is not comparable with the
    # densities above and is deliberately not returned.
    list(mean = predict(fit, newdata = test, type = "mean"), score = NA_real_)
  }
)

OUT <- file.path("_dev", "positive-results.rds")

out <- list()
n_total <- length(truths) * reps * length(fits)

pr <- prog_init(total = n_total, title = "Positive continuous families",
                unit = "fit", kind = "simulation", workers = 1L,
                command = "Rscript _dev/positive-sim.R")

# A crash still closes the run, so the widget and `progress-status` report a
# failure rather than a job that looks as though it is still going.
on.exit(prog_end(pr, "failed", "aborted before the last cell"), add = TRUE)

checkpoint <- function(complete) {
  saveRDS(list(res = do.call(rbind, out), complete = complete,
               done = length(out), total = n_total, reps = reps,
               n_train = n_train, n_test = n_test, n_bins = n_bins),
          OUT)
}

for (truth in names(truths)) {
  for (rep in seq_len(reps)) {
    set.seed(1000L * match(truth, names(truths)) + rep)

    x_train <- sim_x(n_train)
    x_test <- sim_x(n_test)
    mu_train <- mean_fun(x_train)
    mu_test <- mean_fun(x_test)

    train <- x_train
    train$y <- truths[[truth]](mu_train, x_train)
    test <- x_test
    test$y <- truths[[truth]](mu_test, x_test)

    for (family in names(fits)) {
      started <- proc.time()[["elapsed"]]
      got <- fits[[family]](train, test)
      seconds <- proc.time()[["elapsed"]] - started

      out[[length(out) + 1L]] <- data.frame(
        truth = truth, family = family, rep = rep,
        rmse = sqrt(mean((got$mean - mu_test)^2)),
        score = got$score,
        seconds = seconds
      )

      prog_tick(pr, i = length(out), secs = seconds,
                label = sprintf("%s / %s rep %d", truth, family, rep))
    }

    # After every replicate rather than after the last one: a killed run leaves
    # every finished fit on disk, and `complete` stops the partial file being
    # read as a whole one.
    checkpoint(FALSE)
  }
}

checkpoint(TRUE)
results <- do.call(rbind, out)

on.exit()
prog_end(pr, "done", sprintf("%d fits over %d truths", length(out),
                             length(truths)))

# The medians rather than the means, because the heavy-tail row's log score is
# dominated by whichever test point landed furthest out.
agg <- aggregate(cbind(rmse, score, seconds) ~ truth + family, results,
                 function(z) stats::median(z, na.rm = TRUE), na.action = NULL)

print(agg[order(agg$truth, agg$family), ], digits = 4)
