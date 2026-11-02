# Does a tighter prior on `dpm_ls()`'s scale forest close its gap to the
# paper's form when the scale is constant, without giving up what it gains when
# the scale varies?
#
# What is tested. `_dev/dpm-ls-forms.R` found `dpm_ls()` behind
# `dpm_ls_component()` on skewed errors with a constant scale (mean RMSE .0900
# against .0792, difference .011, SE .002), and `_dev/dpm-ls-homoskedastic.R`
# traced it to spurious x-variation in the scale forest: with that forest
# emptied the RMSE was .0758, and the fitted log_sd varied with sd .059 across
# test points where the truth is flat. Under a skewed W the mean sits a
# distance proportional to exp(g) from the bulk of the error, so that variation
# moves the mean. The scale forest's leaf standard deviation is drawn under a
# half-Cauchy prior whose median is 3 * .5 / (k sqrt(10)), .237 at the default
# k = 2, borrowed from `gaussian_ls()`. The claim is that a smaller median, or
# the Dirichlet sparsity prior on that forest alone, suppresses the spurious
# variation enough to close the gap, and that the half-Cauchy's tail leaves
# room for a real scale signal, so the heteroskedastic scenarios keep their
# gains. It could fail either way: the tighter prior may not reach the gap, or
# it may cost the scenarios where the scale really varies.
#
# What is measured. The main study's 120 data sets on its seeds (six scenarios
# by 20; see that file for the design), refit with `dpm_ls()` under three
# settings of the scale forest, everything else as there:
#
#   k4       k = c(mean = 2, log_sd = 4): the median halved, to .119.
#   k8       k = c(mean = 2, log_sd = 8): quartered, to .059.
#   sparse   sparsity = c(mean = FALSE, log_sd = TRUE), k at its default.
#
# Each is compared within data set to the main study's stored fits of
# `dpm_ls()` at the default (`whole`), which this script's `check` run refits
# and must reproduce exactly, and of `dpm_ls_component()` (`component`). The
# columns are the main study's (score, rmse, mean_cover, cover_lo, cover_hi,
# seconds) and g_spread, the sd over test points of the posterior mean of
# log_sd, which says whether the prior did what it was meant to.
#
# THE COMPARISONS TO READ: rmse on `homoskedastic` against `component`, and
# score on `normal`, `heavy`, `skewed` and `bimodal` against `whole`.
#
# What each outcome would mean.
#   - A setting qualifies to become `dpm_ls()`'s default if on `homoskedastic`
#     its rmse is no worse than `component`'s beyond two standard errors, and
#     on each of `normal`, `heavy`, `skewed` and `bimodal` its score is no worse
#     than `whole`'s beyond two standard errors (so it keeps the lead over
#     `component` there). If more than one qualifies, the one closest to the
#     present default is preferred (k4 over k8), and `sparse` is preferred to
#     either only if it is ahead of them on the homoskedastic rmse, since it is
#     the larger change to how the forest behaves.
#   - If none qualifies but one closes the homoskedastic gap at a cost on the
#     heteroskedastic scenarios, the trade is real: the default stays, and the
#     documentation says that a constant scale is `dpm()`'s case.
#   - If none moves the homoskedastic rmse at all, while g_spread does fall,
#     the account in `_dev/dpm-ls-homoskedastic.R` is incomplete and needs
#     revisiting before anything else is tried.
#   - `component` (the paper's model as truth) is reported but not part of the
#     rule; the whole-error form cannot represent it.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/dpm-ls-scale-prior.R [reps] [workers] [check]
#   where `check` (any third argument) also refits `whole` and stops unless it
#   matches the stored fit to the last digit.
# Writes:   _dev/dpm-ls-scale-prior/<scenario>-<rep>.rds per data set, and
#           _dev/dpm-ls-scale-prior.rds with all of them at the end. Reads
#           _dev/dpm-ls-forms/, which `_dev/dpm-ls-forms.R` writes.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

lib <- Sys.getenv("BARTISAN_LIB")

if (nzchar(lib)) {
  .libPaths(c(lib, .libPaths()))
}

library(bartisan)

args <- commandArgs(TRUE)
reps <- if (length(args) >= 1L) as.integer(args[1L]) else 20L
workers <- if (length(args) >= 2L) as.integer(args[2L]) else 6L
check <- length(args) >= 3L

study <- local({
  n_train <- 500L
  n_test <- 1000L
  units_dir <- "_dev/dpm-ls-scale-prior"
  stored_dir <- "_dev/dpm-ls-forms"

  # The generators of `_dev/dpm-ls-forms.R`, verbatim, so that a seed gives the
  # same data set here as there.
  sim_x <- function(n) {
    as.data.frame(matrix(stats::runif(5L * n), n, 5L,
                         dimnames = list(NULL, paste0("x", 1:5))))
  }

  mean_fun <- function(x) 1 + 0.8 * sin(pi * x$x1) + 0.6 * x$x2
  sd_fun <- function(x) exp(-1.2 + 1.2 * x$x3)

  std_gamma <- function(n) (stats::rgamma(n, shape = 1.5) - 1.5) / sqrt(1.5)

  scenarios <- list(
    normal = function(x) {
      s <- sd_fun(x)
      list(e = s * stats::rnorm(nrow(x)), sd = s)
    },
    skewed = function(x) {
      s <- sd_fun(x)
      list(e = s * std_gamma(nrow(x)), sd = s)
    },
    heavy = function(x) {
      s <- sd_fun(x)
      list(e = s * stats::rt(nrow(x), df = 3) / sqrt(3), sd = s)
    },
    bimodal = function(x) {
      s <- sd_fun(x)
      side <- ifelse(stats::rbinom(nrow(x), 1L, 0.5) == 1L, 0.9, -0.9)
      list(e = s * (side + stats::rnorm(nrow(x), sd = sqrt(0.19))), sd = s)
    },
    component = function(x) {
      s <- sd_fun(x)
      side <- ifelse(stats::rbinom(nrow(x), 1L, 0.5) == 1L, 0.6, -0.6)
      list(e = side + 0.6 * s * stats::rnorm(nrow(x)),
           sd = 0.6 * sqrt(1 + s^2))
    },
    homoskedastic = function(x) {
      s <- rep.int(0.6, nrow(x))
      list(e = s * std_gamma(nrow(x)), sd = s)
    }
  )

  predictors <- paste0("x", 1:5)
  formula <- list(stats::reformulate(predictors, "y"),
                  log_sd = stats::reformulate(predictors))

  settings <- list(
    whole = list(),
    k4 = list(k = c(mean = 2, log_sd = 4)),
    k8 = list(k = c(mean = 2, log_sd = 8)),
    sparse = list(sparsity = c(mean = FALSE, log_sd = TRUE))
  )

  fit_one <- function(setting, train) {
    # Namespaced, because the worker is sent a wrapper around this function
    # and *future* does not see past it to attach the package.
    args <- c(list(formula, data = train, family = bartisan::dpm_ls(),
                   num_trees = c(50L, 10L), num_burn = 500L, num_draws = 500L,
                   verbose = FALSE),
              settings[[setting]])
    do.call(bartisan::bartisan, args)
  }

  # The main study's score_fit(), with the spread of the fitted scale added.
  score_fit <- function(fit, test, truth) {
    mu <- stats::predict(fit, newdata = test, type = "response", draws = TRUE)
    mean_hat <- colMeans(mu)
    mean_lo <- apply(mu, 2L, stats::quantile, 0.05)
    mean_hi <- apply(mu, 2L, stats::quantile, 0.95)

    dens <- stats::predict(fit, newdata = test, type = "density", log = TRUE)

    reps <- bartisan:::posterior_predict.bartisan_fit(fit, newdata = test)
    lo <- apply(reps, 2L, stats::quantile, 0.05)
    hi <- apply(reps, 2L, stats::quantile, 0.95)
    inside <- test$y >= lo & test$y <= hi

    breaks <- unique(stats::quantile(truth$sd, c(0, 1 / 3, 2 / 3, 1)))
    third <- {
      if (length(breaks) < 4L) rep(NA_integer_, length(truth$sd))
      else cut(truth$sd, breaks, include.lowest = TRUE, labels = FALSE)
    }

    log_sd <- stats::predict(fit, newdata = test, type = "link")[, 2L]

    c(score = mean(dens),
      rmse = sqrt(mean((mean_hat - truth$mean)^2)),
      mean_cover = mean(truth$mean >= mean_lo & truth$mean <= mean_hi),
      cover_lo = mean(inside[which(third == 1L)]),
      cover_hi = mean(inside[which(third == 3L)]),
      g_spread = stats::sd(log_sd))
  }

  one_dataset <- function(unit, check) {
    stored <- readRDS(file.path(stored_dir, sprintf("%s-%02d.rds",
                                                    unit$scenario, unit$rep)))

    set.seed(unit$seed)
    train <- sim_x(n_train)
    test <- sim_x(n_test)
    draw_train <- scenarios[[unit$scenario]](train)
    draw_test <- scenarios[[unit$scenario]](test)
    train$y <- mean_fun(train) + draw_train$e
    test$y <- mean_fun(test) + draw_test$e
    truth <- list(mean = mean_fun(test), sd = draw_test$sd)

    # The fits draw from the stream the data left behind, as the main study's
    # first fit did, so only `whole` run first reproduces it.
    names <- c(if (check) "whole", "k4", "k8", "sparse")

    rows <- lapply(names, function(name) {
      started <- Sys.time()
      fit <- fit_one(name, train)
      seconds <- as.numeric(Sys.time() - started, units = "secs")
      data.frame(scenario = unit$scenario, rep = unit$rep, fit = name,
                 t(score_fit(fit, test, truth)), seconds = seconds)
    })

    out <- do.call(rbind, rows)

    if (check) {
      mine <- out[out$fit == "whole", c("score", "rmse", "mean_cover")]
      theirs <- stored[stored$fit == "whole", c("score", "rmse", "mean_cover")]

      if (!isTRUE(all.equal(unlist(mine), unlist(theirs), tolerance = 0))) {
        stop("`whole` does not reproduce the stored fit for ", unit$scenario,
             " ", unit$rep)
      }
    }

    # The stored fits join the new ones, with g_spread missing since the main
    # study did not record it.
    keep <- stored[stored$fit %in% c("whole", "component"), ]
    keep$g_spread <- NA_real_
    out <- rbind(keep[names(out)][!(keep$fit %in% out$fit), ], out)

    saveRDS(out, file.path(units_dir, sprintf("%s-%02d.rds", unit$scenario,
                                              unit$rep)))
    out
  }

  environment()
})

dir.create(study$units_dir, showWarnings = FALSE)

# The main study's units, rep fastest over its scenarios in its order, seeded
# from 20261101 + 1.
all_units <- expand.grid(rep = seq_len(20L),
                         scenario = c("normal", "skewed", "heavy", "bimodal",
                                      "component", "homoskedastic"),
                         stringsAsFactors = FALSE)
all_units$seed <- 20261101L + seq_len(nrow(all_units))
units <- all_units[all_units$rep <= reps, ]
units <- split(units, seq_len(nrow(units)))

pr <- prog_init(total = length(units), title = "dpm_ls() scale prior",
                unit = "data set", workers = workers, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

future::plan(future::multisession, workers = workers)

results <- prog_future_lapply(pr, units,
                              function(u) study$one_dataset(u, check),
                              label = function(u) paste(u$scenario, u$rep))

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(res = results, complete = TRUE), "_dev/dpm-ls-scale-prior.rds")

on.exit()
prog_end(pr, "done")

# Within data set: each setting against `whole` on score and against
# `component` on rmse, with standard errors over the data sets.
paired <- function(metric, against) {
  rows <- list()

  for (scenario in unique(results$scenario)) {
    d <- results[results$scenario == scenario, ]
    base <- d[d$fit == against, c("rep", metric)]

    for (setting in c("k4", "k8", "sparse")) {
      o <- merge(d[d$fit == setting, c("rep", metric)], base, by = "rep")
      diff <- o[[2L]] - o[[3L]]
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = scenario, setting = setting, n = nrow(o),
        diff = mean(diff), se = stats::sd(diff) / sqrt(nrow(o)))
    }
  }

  out <- do.call(rbind, rows)
  names(out)[4:5] <- paste0(metric, "_vs_", against, c("", "_se"))
  out
}

print(paired("score", "whole"), digits = 3)
print(paired("rmse", "component"), digits = 3)

means <- stats::aggregate(cbind(score, rmse, mean_cover, cover_lo, cover_hi,
                                g_spread, seconds) ~ scenario + fit, results,
                          mean, na.action = stats::na.pass)
print(means[order(means$scenario, means$fit), ], digits = 3)
