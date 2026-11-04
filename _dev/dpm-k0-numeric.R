# Drawing the Dirichlet process baseline's k_0 and lambda, and nu = 3, on
# uncensored responses: `dpm()` and `dpm_ls()`.
#
# What is tested. `_dev/dpm-k0-oracle.R` and `_dev/dpm-k0-reps.R` test, on
# survival times, whether drawing the baseline's k_0 and lambda (and setting
# nu = 3 in place of 10) makes the mixture's intervals cover. A default changes
# for every Dirichlet process family at once, so the same arms are needed on
# the uncensored families, where the quantities people read are the
# conditional mean, the predictive distribution, and the error density. The
# claim: on errors that are not normal in shape the drawn arms raise the
# coverage of `error_density()`'s intervals and the log predictive score, and on
# no error shape do they cost the conditional mean's coverage or accuracy, or
# the score, beyond two SE. It could fail by the looser baseline letting the
# mixture's location trade against the mean forest's level (mean coverage or
# RMSE), or by nu = 3 opening narrow components on noise (a lower score on the
# normal and heavy-tailed errors).
#
# A second question, added before this run at the maintainer's direction:
# `dpm()` stays the default for continuous responses as long as its coverage is
# at least `gaussian()`'s on every error shape, and where both fall short the
# shortfall is documented as a limitation. So each data set is also fitted by
# the parametric counterpart, `gaussian()` for `dpm()` and `gaussian_ls()` for
# `dpm_ls()`, and the claim is that the default arm's coverage of every
# quantity below is at least the parametric family's on every shape, within
# two SE. It fails wherever the parametric family covers better: on the normal
# error, where it is the true model, that is the price of the default; on the
# others it would be a reason to reconsider it.
#
# The measurement. The design of `_dev/dpm-ls-forms.R`: 500 training and 1000
# test observations, five uniform predictors, E[y | x] = 1 + 0.8 sin(pi x1) +
# 0.6 x2, and W of mean 0 and variance 1 in four shapes: normal, skewed (a
# standardized gamma(1.5)), heavy (a standardized t_3), and bimodal (+-0.9 +
# N(0, .19)). `dpm()` gets the error 0.6 W; `dpm_ls()` gets s(x) W with
# s(x) = exp(-1.2 + 1.2 x3). Twenty data sets per shape and family, fits at 500
# warmup + 500 draws, one chain, 50 trees on the mean and 10 on the scale. Arms:
# `fixed` (the defaults), `both` (k_0 and lambda drawn under gamma priors of
# shape 1.5 centered on the fixed values), `both_nu3` (the same, nu = 3), and
# `parametric` (`gaussian()` or `gaussian_ls()`).
# Per fit, on the test set: `score`, the mean log predictive density;
# `rmse` and `mean_cover`, the posterior mean's error and the 90% intervals'
# coverage of E[y | x]; `pred_cover`, the coverage of 90% posterior predictive
# intervals for y; and for `dpm()` and `gaussian()`, `dens_cover`, the share of
# 101 points over the middle 96% of the true error distribution whose 90%
# interval holds the true density (from `error_density()`, or from the normal
# density at each draw's sigma), and `exceed_cover`, the coverage of 90%
# intervals for P(Y > c | x) at c = E[y | x] + the true error's 10th, 25th,
# 50th, 75th and 90th percentiles, at every test point: the continuous
# response's counterpart of S(t | x), whose true values are .9, .75, .5, .25
# and .1. Arms are compared within data set.
# Read `dens_cover` and `score` first, then the paired changes in `mean_cover`
# and `rmse`.
#
# What each outcome would mean.
#   - The claim holds for an arm: that arm is consistent across the censored
#     and uncensored families and stays a candidate default.
#   - An arm costs the conditional mean (coverage or RMSE beyond two SE) on any
#     shape: it cannot be the default as it stands, whatever it does for the
#     error distribution, since the conditional mean is what most fits are read
#     for.
#   - `dens_cover` does not improve here though it did on survival times: the
#     censored result does not carry over, and the reason is looked for before
#     any default changes.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/dpm-k0-numeric.R [reps] [workers]
# Writes:   _dev/dpm-k0-numeric.rds, and one file per data set in
#           _dev/dpm-k0-numeric/ as it finishes.

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

study <- local({
  n_train <- 500L
  n_test <- 1000L
  units_dir <- "_dev/dpm-k0-numeric"
  level <- 0.90

  sim_x <- function(n) {
    as.data.frame(matrix(stats::runif(5L * n), n, 5L,
                         dimnames = list(NULL, paste0("x", 1:5))))
  }

  mean_fun <- function(x) 1 + 0.8 * sin(pi * x$x1) + 0.6 * x$x2
  sd_fun <- function(x) exp(-1.2 + 1.2 * x$x3)

  # Each shape as a generator of W and its density.
  shapes <- list(
    normal = list(draw = function(n) stats::rnorm(n),
                  dens = function(w) stats::dnorm(w)),
    skewed = list(
      draw = function(n) (stats::rgamma(n, shape = 1.5) - 1.5) / sqrt(1.5),
      dens = function(w) stats::dgamma(1.5 + sqrt(1.5) * w, 1.5) * sqrt(1.5)),
    heavy = list(draw = function(n) stats::rt(n, df = 3) / sqrt(3),
                 dens = function(w) stats::dt(w * sqrt(3), 3) * sqrt(3)),
    bimodal = list(
      draw = function(n) {
        ifelse(stats::rbinom(n, 1L, 0.5) == 1L, 0.9, -0.9) +
          stats::rnorm(n, sd = sqrt(0.19))
      },
      dens = function(w) {
        0.5 * stats::dnorm(w, 0.9, sqrt(0.19)) +
          0.5 * stats::dnorm(w, -0.9, sqrt(0.19))
      })
  )

  # The middle 96% of each W, from a large sample, for the density grid.
  ranges <- lapply(shapes, function(sh) {
    set.seed(1)
    stats::quantile(sh$draw(1e6), c(0.02, 0.98), names = FALSE)
  })

  # The arms name `update_baseline = FALSE` explicitly where they mean the fixed
  # baseline, which was the default when this was run and is not since
  # 2026-11-03.
  arms <- list(fixed = list(update_baseline = FALSE),
               both = list(k0_shape = 1.5, lambda_shape = 1.5),
               both_nu3 = list(k0_shape = 1.5, lambda_shape = 1.5, nu = 3),
               parametric = NULL)

  # The true error's percentiles, for the exceedance thresholds.
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  percentiles <- lapply(shapes, function(sh) {
    set.seed(2)
    stats::quantile(sh$draw(1e6), probs, names = FALSE)
  })

  predictors <- paste0("x", 1:5)
  mean_formula <- stats::reformulate(predictors, "y")
  scale_formula <- stats::reformulate(predictors)

  fit_one <- function(family, arm, train) {
    two <- family == "dpm_ls"

    if (arm == "parametric") {
      fam <- if (two) gaussian_ls() else stats::gaussian()
    }
    else {
      fam <- if (two) dpm_ls() else dpm()
      fam[names(arms[[arm]])] <- arms[[arm]]
    }

    formula <- if (two) list(mean_formula, log_sd = scale_formula) else mean_formula
    bartisan(formula, data = train, family = fam,
             num_trees = if (two) c(50L, 10L) else 50L, num_burn = 500L,
             num_draws = 500L, verbose = FALSE)
  }

  score_fit <- function(fit, test, truth, unit) {
    mu <- predict(fit, newdata = test, type = "response", draws = TRUE)
    mean_hat <- colMeans(mu)
    mean_lo <- apply(mu, 2L, stats::quantile, (1 - level) / 2)
    mean_hi <- apply(mu, 2L, stats::quantile, 1 - (1 - level) / 2)

    dens <- predict(fit, newdata = test, type = "density", log = TRUE)

    reps <- bartisan:::posterior_predict.bartisan_fit(fit, newdata = test)
    lo <- apply(reps, 2L, stats::quantile, (1 - level) / 2)
    hi <- apply(reps, 2L, stats::quantile, 1 - (1 - level) / 2)

    dens_cover <- NA_real_
    exceed_cover <- NA_real_
    draws <- seq_len(nrow(mu))
    gaussian_fit <- identical(fit[["family"]][["family"]], "gaussian")

    if (unit$family == "dpm") {
      sh <- shapes[[unit$shape]]
      at <- 0.6 * seq(ranges[[unit$shape]][1L], ranges[[unit$shape]][2L],
                      length.out = 101L)
      true_dens <- sh$dens(at / 0.6) / 0.6

      if (gaussian_fit) {
        sigma <- fit[["aux"]][, "sigma"]
        dd <- t(vapply(sigma, function(s) stats::dnorm(at, 0, s),
                       numeric(length(at))))
        ed <- list(lower = apply(dd, 2L, stats::quantile, (1 - level) / 2),
                   upper = apply(dd, 2L, stats::quantile, 1 - (1 - level) / 2))
      }
      else {
        ed <- error_density(fit, at = at, level = level)
      }

      dens_cover <- mean(true_dens >= ed$lower & true_dens <= ed$upper)

      # P(Y > c | x) at each draw, for every threshold at every test point.
      cuts <- outer(0.6 * percentiles[[unit$shape]], truth$mean, "+")
      exceed <- vapply(draws, function(s) {
        gap <- sweep(cuts, 2L, mu[s, ])
        if (gaussian_fit) {
          stats::pnorm(gap, 0, fit[["aux"]][s, "sigma"], lower.tail = FALSE)
        }
        else {
          bartisan:::dpm_survival(fit, s, gap)
        }
      }, cuts)
      lo_e <- apply(exceed, c(1L, 2L), stats::quantile, (1 - level) / 2)
      hi_e <- apply(exceed, c(1L, 2L), stats::quantile, 1 - (1 - level) / 2)
      p_true <- matrix(1 - probs, nrow(cuts), ncol(cuts))
      exceed_cover <- mean(p_true >= lo_e & p_true <= hi_e)
    }

    mixture <- !is.null(fit[["aux"]]) && "clusters" %in% colnames(fit[["aux"]])

    c(score = mean(dens),
      rmse = sqrt(mean((mean_hat - truth$mean)^2)),
      mean_cover = mean(truth$mean >= mean_lo & truth$mean <= mean_hi),
      pred_cover = mean(test$y >= lo & test$y <= hi),
      dens_cover = dens_cover,
      exceed_cover = exceed_cover,
      k_0 = if (mixture) mean(bartisan:::dpm_hyper(fit, draws, "k_0")) else NA,
      lambda = if (mixture) mean(bartisan:::dpm_hyper(fit, draws, "lambda"))
               else NA,
      clusters = if (mixture) mean(fit[["aux"]][, "clusters"]) else NA)
  }

  one_dataset <- function(unit) {
    set.seed(unit$seed)
    train <- sim_x(n_train)
    test <- sim_x(n_test)
    s_train <- if (unit$family == "dpm") 0.6 else sd_fun(train)
    s_test <- if (unit$family == "dpm") 0.6 else sd_fun(test)
    train$y <- mean_fun(train) + s_train * shapes[[unit$shape]]$draw(n_train)
    test$y <- mean_fun(test) + s_test * shapes[[unit$shape]]$draw(n_test)
    truth <- list(mean = mean_fun(test))

    rows <- lapply(names(arms), function(arm) {
      started <- Sys.time()
      fit <- fit_one(unit$family, arm, train)
      seconds <- as.numeric(Sys.time() - started, units = "secs")
      data.frame(family = unit$family, shape = unit$shape, rep = unit$rep,
                 arm = arm, t(score_fit(fit, test, truth, unit)),
                 seconds = seconds)
    })

    out <- do.call(rbind, rows)
    saveRDS(out, file.path(units_dir, sprintf("%s-%s-%02d.rds", unit$family,
                                              unit$shape, unit$rep)))
    out
  }

  environment()
})

dir.create(study$units_dir, showWarnings = FALSE)

units <- expand.grid(rep = seq_len(reps), shape = names(study$shapes),
                     family = c("dpm", "dpm_ls"), stringsAsFactors = FALSE)
units$seed <- 20261103L + seq_len(nrow(units))
units <- split(units, seq_len(nrow(units)))

pr <- prog_init(total = length(units), title = "Drawing k_0 and lambda, uncensored",
                unit = "data set", workers = workers, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

future::plan(future::multisession, workers = workers)

results <- prog_future_lapply(pr, units, study$one_dataset,
                              label = function(u) {
                                paste(u$family, u$shape, u$rep)
                              })

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(main = results, complete = TRUE, when = Sys.time()),
        "_dev/dpm-k0-numeric.rds")
on.exit()
prog_end(pr, "done", sprintf("%d data sets", length(units)))
cat("wrote _dev/dpm-k0-numeric.rds\n")
