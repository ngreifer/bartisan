# Where `dpm_aft()`'s over-coverage on a linear mean comes from.
#
# What is tested. On the survival part of `_dev/dpm-packages-coverage.R`
# (log T = 1 + 0.8 x1 - 0.6 x2 + W, n = 700, 30% censored), 90% intervals for
# S(t | x) from `dpm_aft()` with its baseline drawn (the default since
# 2026-11-03) covered .963 on a log-normal error and .947 on a bimodal one,
# with intervals .135 and .088 wide, against `anovaDDP()`'s .891 and .892 at
# .077 and .066. `anovaDDP()` assumes a linear mean, which is the truth here,
# and centers its priors on a log-normal `survreg()` fit to the same data;
# `dpm_aft()`'s mean is a sum of trees. The claim: the extra width is the
# trees' uncertainty about a mean they do not assume linear, not the mixture.
# Three predictions, each of which could fail: (i) with the true mean given as
# an offset, so that only the error is estimated, `dpm_aft()`'s coverage of
# S(t | x) falls to near nominal; (ii) `lognormal_aft()`, whose error model is
# exactly right on the log-normal truth and whose mean is the same sum of
# trees, over-covers there as much as `dpm_aft()` does; (iii) both families'
# 90% intervals for E[log T | x] over-cover the true linear mean.
#
# The measurement. Replicates 1 to 20 of the log-normal and bimodal truths of
# that study, the same seeds and draws. Arms: `full` (`dpm_aft()` at its
# defaults, 50 trees, 500 + 500), `oracle` (`dpm_aft()` with the true mean as
# an offset and a forest on x3, which the mean does not use), and `lognormal`
# (`lognormal_aft()` at its defaults). Per fit: `cover` and `width` of 90%
# intervals for S(t | x) on the same 100 test subjects by 20 times; for `full`
# and `lognormal`, `mean_cover` and `mean_width` of 90% intervals for
# E[log T | x] at the test subjects against the truth 1 + 0.8 x1 - 0.6 x2
# (the error has mean zero in both truths).
#
# What each outcome would mean.
#   - (i) to (iii) hold: the over-coverage is the trees' honest uncertainty
#     about the mean, and nominal coverage against a linear truth is bought by
#     assuming linearity, which a BART family does not do. Nothing in
#     `anovaDDP()` carries over, and the documentation should say that on a
#     mean that a parametric model gets right, BART's intervals are
#     conservative.
#   - (i) fails, the oracle still over-covering: the mixture's own intervals
#     are too wide, and the drawn baseline's priors come under suspicion.
#   - (ii) fails, `lognormal_aft()` near nominal: the width is the mixture's
#     contribution, not the trees'.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#           PKG_LIB=<library with spBayesSurv and the other packages> \
#             Rscript _dev/linear-mean-coverage.R [workers]
# Writes:   _dev/linear-mean-coverage.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

libs <- c(Sys.getenv("BARTISAN_LIB"), Sys.getenv("PKG_LIB"))
.libPaths(c(libs[nzchar(libs)], .libPaths()))

args <- commandArgs(TRUE)
workers <- if (length(args) >= 1L) as.integer(args[1L]) else 6L

# The data come from `_dev/dpm-packages-coverage.R`'s own definitions, so the
# draws are that study's.
lines <- readLines("_dev/dpm-packages-coverage.R")
eval(parse(text = lines[grep("^study <- local", lines):(grep("^dir.create", lines) - 1L)]))

units <- expand.grid(rep = 1:20, shape = c("lognormal", "bimodal"),
                     stringsAsFactors = FALSE)

# That study numbered its 80 density units first and then the survival units
# by shape (lognormal, bimodal, skewed), rep fastest, with seed 20261104 + k.
units$seed <- 20261104L + 80L +
  (match(units$shape, names(study$errors)) - 1L) * 20L + units$rep
units <- split(units, seq_len(nrow(units)))

one_unit <- function(unit) {
  suppressMessages(library(survival))
  errors <- study$errors
  level <- 0.90

  set.seed(unit$seed)
  n <- 700L
  mk <- function(n) {
    data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
               x3 = stats::runif(n))
  }
  m <- function(d) 1 + 0.8 * d$x1 - 0.6 * d$x2

  train <- mk(n)
  test <- mk(100L)
  log_t <- m(train) + errors[[unit$shape]]$draw(n)

  # The true mean as a column, since `bartisan()` evaluates `offset` in the
  # data (`predict()` takes the vector). Adding it draws nothing, so the data are unchanged.
  train$off <- m(train)
  test$off <- m(test)
  log_c <- stats::rnorm(n, mean(log_t) + 1.2, 1.5)
  train$time <- exp(pmin(log_t, log_c))
  train$status <- as.numeric(log_t <= log_c)

  times <- exp(stats::quantile(log_t, seq(0.05, 0.95, length.out = 20L),
                               names = FALSE))
  truth <- t(vapply(log(times), function(l) {
    errors[[unit$shape]]$surv(l - m(test))
  }, numeric(nrow(test))))

  ctrl <- bartisan::bartisan_control(num_trees = 50, num_burn = 500,
                                     num_draws = 500)
  q <- c((1 - level) / 2, 1 - (1 - level) / 2)

  arms <- list(
    full = list(formula = survival::Surv(time, status) ~ x1 + x2 + x3,
                family = bartisan::dpm_aft(), offset = FALSE),
    oracle = list(formula = survival::Surv(time, status) ~ x3,
                  family = bartisan::dpm_aft(), offset = TRUE),
    lognormal = list(formula = survival::Surv(time, status) ~ x1 + x2 + x3,
                     family = bartisan::lognormal_aft(), offset = FALSE))

  rows <- lapply(names(arms), function(arm) {
    a <- arms[[arm]]
    if (a$offset) {
      fit <- bartisan::bartisan(a$formula, data = train, family = a$family,
                                offset = off, control = ctrl, verbose = FALSE)
      drawn <- predict(fit, newdata = test, type = "survival", times = times,
                       draws = TRUE, offset = test$off)
    }
    else {
      fit <- bartisan::bartisan(a$formula, data = train, family = a$family,
                                control = ctrl, verbose = FALSE)
      drawn <- predict(fit, newdata = test, type = "survival", times = times,
                       draws = TRUE)
    }
    lo <- t(apply(drawn, c(2L, 3L), stats::quantile, q[1L]))
    hi <- t(apply(drawn, c(2L, 3L), stats::quantile, q[2L]))

    mean_cover <- NA_real_
    mean_width <- NA_real_

    if (!a$offset) {
      eta <- predict(fit, newdata = test, type = "link", draws = TRUE)
      elo <- apply(eta, 2L, stats::quantile, q[1L])
      ehi <- apply(eta, 2L, stats::quantile, q[2L])
      mean_cover <- mean(m(test) >= elo & m(test) <= ehi)
      mean_width <- mean(ehi - elo)
    }

    data.frame(truth = unit$shape, rep = unit$rep, arm = arm,
               cover = mean(truth >= lo & truth <= hi),
               width = mean(hi - lo), mean_cover = mean_cover,
               mean_width = mean_width)
  })

  do.call(rbind, rows)
}

pr <- prog_init(total = length(units), title = "Linear mean, coverage",
                unit = "data set", workers = workers, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

future::plan(future::multisession, workers = workers)

results <- prog_future_lapply(pr, units, one_unit,
                              label = function(u) paste(u$shape, u$rep))

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(main = results, complete = TRUE, when = Sys.time()),
        "_dev/linear-mean-coverage.rds")
on.exit()
prog_end(pr, "done", sprintf("%d data sets", length(units)))
cat("wrote _dev/linear-mean-coverage.rds\n")
