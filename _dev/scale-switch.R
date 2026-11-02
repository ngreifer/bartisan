# Can the data switch `dpm_ls()`'s scale forest off when the spread is constant?
#
# What is tested. On a skewed error of constant spread, `dpm_ls()` estimates the
# conditional mean worse than `dpm()` or the component form (`_dev/dpm-ls-forms.R`,
# RMSE .090 against .079), and `_dev/nftbart-compare.R` found nftbart, which
# fits the same model, worse by about as much, so the cost is the model's and
# not this implementation's. The cost comes from the scale forest's chance
# variation; with the forest emptied (`log_sd = ~ 1`) the RMSE was .076
# (`_dev/dpm-ls-homoskedastic.R`). The remedy proposed is a two-point prior on
# the whole forest, "constant" against "varies", so that the fit averages over
# the two models by their posterior probabilities. Built into the sampler, an
# indicator on a whole forest mixes badly (conditional on a forest fitted to
# the data, switching it off always looks costly), so this study tests the
# premise first, at the level of whole fits: the full fit and the flat fit,
# averaged with weights from PSIS-LOO, which approximate the posterior
# probabilities of the two models by their estimated predictive performance.
# The claim, which could be false, is that those weights go to the flat model
# when the spread is constant and to the full one when it varies, so that the
# average gets the flat fit's RMSE in the first case and keeps the full fit's
# score and RMSE in the second. The one earlier hint is in its favor: under
# constant spread the flat fit's held-out score was the better one, by .0045
# per observation (SE .0012), so held-out prediction can tell the two apart.
#
# What is measured. The 120 data sets of `_dev/dpm-ls-forms.R`, reproduced
# exactly from its seeds (six error scenarios, 500 training and 1000 test
# points). Fits, all with hard rules and 500 + 500 iterations for speed: `full`,
# `dpm_ls()` with the scale forest on all five predictors; `flat`, `dpm_ls()`
# with `log_sd = ~ 1`; `component`, `dpm_ls_component()` for reference; and two
# averages of `full` and `flat` from PSIS-LOO, `stack` with stacking weights and
# `bma` with pseudo-BMA+ weights; and the same two from 5-fold cross-validation
# on folds shared by both fits, `kstack` and `kbma`. The K-fold pair was added
# after a smoke test: on the first constant-spread data set the flat fit's
# held-out score was better by .017 per observation, yet PSIS-LOO put the full
# fit ahead by 9.8 in elpd and gave the flat one 2% of the weight, as if LOO
# were optimistic about a fit with a second forest. An average's mean is the
# weighted mean of the two posterior means, its density the weighted mixture,
# and its draws for intervals are taken from the two fits in proportion to the
# weights. Per fit:
#
#   rmse        RMSE of the posterior mean of E[y | x] against the truth. THE
#               COLUMN TO READ on `homoskedastic`.
#   score       mean held-out log predictive density per test observation. THE
#               COLUMN TO READ on the other five scenarios.
#   mean_cover  share of test points whose 90% interval for E[y | x] holds the
#               truth (nominal .90).
#   cover_lo, cover_hi  coverage of 90% predictive intervals in the lowest and
#               highest third of the true conditional SD.
#   w_flat      the weight an average gives the flat fit (for `stack`, `bma`).
#
# What each outcome would mean.
#   - The switch works if `stack` (or `bma`) has an RMSE on `homoskedastic`
#     within two SE of `flat`'s and below `full`'s, and on the five other
#     scenarios a score and an RMSE within two SE of `full`'s. Then a two-point
#     prior is worth building into the sampler, or the average is worth
#     offering as it stands, since it needs nothing new in the engine.
#   - If the LOO weights stay low under constant spread but the K-fold ones go
#     to the flat fit, held-out data can tell the two apart and PSIS-LOO cannot:
#     the switch has to be decided by cross-validation, and an indicator inside
#     the sampler, which sees only the fitted forest, would share LOO's bias.
#   - If both stay low, the data cannot tell the two models apart well enough
#     at this sample size, a prior indicator would not switch the forest off
#     either, and option 2 (capping how much the scale forest can vary) is next.
#   - If `w_flat` is high where the spread varies but only a little (`normal`
#     and `heavy` have the same spread function as the rest, so this would be
#     about the error's shape), the switch costs score there, and that cost is
#     what has to be weighed against the RMSE it saves.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/scale-switch.R [reps] [workers]
# Writes:   _dev/scale-switch/<scenario>-<rep>.rds per data set, and
#           _dev/scale-switch.rds with all of them at the end.

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

# Everything a worker needs travels in this environment with `one_dataset()`.
study <- local({
  n_train <- 500L
  n_test <- 1000L
  units_dir <- "_dev/scale-switch"

  # The generators and their draw order are `_dev/dpm-ls-forms.R`'s, so each
  # seed gives the same data set there and here.
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
  mean_formula <- stats::reformulate(predictors, "y")
  scale_formula <- stats::reformulate(predictors)

  # The formulas take this frame as their environment, where `kfold()` looks
  # for `train` when it refits.
  fit_one <- function(family, scale, train) {
    mean_f <- mean_formula
    environment(mean_f) <- environment()
    environment(scale) <- environment()
    bartisan::bartisan(list(mean_f, log_sd = scale), data = train,
                       family = family, num_trees = c(50L, 10L),
                       gate = "hard", num_burn = 500L, num_draws = 500L,
                       verbose = FALSE)
  }

  # What scoring needs from a fit, kept so the averages can be formed from it.
  parts <- function(fit, test) {
    list(mu = stats::predict(fit, newdata = test, type = "response",
                             draws = TRUE),
         dens = colMeans(stats::predict(fit, newdata = test,
                                        type = "density", draws = TRUE)),
         reps = bartisan:::posterior_predict.bartisan_fit(fit, newdata = test),
         loo = suppressWarnings(loo::loo(fit)),
         fit = fit)
  }

  # Draws from two fits in proportion to a weight, spaced evenly through each
  # chain.
  pool <- function(a, b, w) {
    s <- nrow(a)
    k <- round(w * s)
    take <- function(m, j) m[unique(round(seq(1, s, length.out = j))), ,
                             drop = FALSE]
    rbind(if (k > 0L) take(a, k), if (k < s) take(b, s - k))
  }

  average <- function(full, flat, w_flat) {
    list(mu = pool(flat$mu, full$mu, w_flat),
         mean = w_flat * colMeans(flat$mu) + (1 - w_flat) * colMeans(full$mu),
         dens = w_flat * flat$dens + (1 - w_flat) * full$dens,
         reps = pool(flat$reps, full$reps, w_flat))
  }

  summarize <- function(p, test, truth) {
    mean_hat <- if (is.null(p$mean)) colMeans(p$mu) else p$mean
    mean_lo <- apply(p$mu, 2L, stats::quantile, 0.05)
    mean_hi <- apply(p$mu, 2L, stats::quantile, 0.95)
    lo <- apply(p$reps, 2L, stats::quantile, 0.05)
    hi <- apply(p$reps, 2L, stats::quantile, 0.95)
    inside <- test$y >= lo & test$y <= hi
    breaks <- unique(stats::quantile(truth$sd, c(0, 1 / 3, 2 / 3, 1)))
    third <- {
      if (length(breaks) < 4L) rep(NA_integer_, length(truth$sd))
      else cut(truth$sd, breaks, include.lowest = TRUE, labels = FALSE)
    }

    c(score = mean(log(p$dens)),
      rmse = sqrt(mean((mean_hat - truth$mean)^2)),
      mean_cover = mean(truth$mean >= mean_lo & truth$mean <= mean_hi),
      cover_lo = mean(inside[which(third == 1L)]),
      cover_hi = mean(inside[which(third == 3L)]))
  }

  one_dataset <- function(unit) {
    set.seed(unit$seed)
    train <- sim_x(n_train)
    test <- sim_x(n_test)
    draw_train <- scenarios[[unit$scenario]](train)
    draw_test <- scenarios[[unit$scenario]](test)
    train$y <- mean_fun(train) + draw_train$e
    test$y <- mean_fun(test) + draw_test$e
    truth <- list(mean = mean_fun(test), sd = draw_test$sd)

    seconds <- c()
    timed <- function(name, expr) {
      started <- Sys.time()
      out <- expr
      seconds[[name]] <<- as.numeric(Sys.time() - started, units = "secs")
      out
    }

    full <- timed("full", parts(fit_one(bartisan::dpm_ls(), scale_formula,
                                        train), test))
    flat <- timed("flat", parts(fit_one(bartisan::dpm_ls(), ~ 1, train), test))
    component <- timed("component",
                       parts(fit_one(bartisan:::dpm_ls_component(),
                                     scale_formula, train), test))

    lpd_loo <- cbind(full$loo$pointwise[, "elpd_loo"],
                     flat$loo$pointwise[, "elpd_loo"])
    w_stack <- loo::stacking_weights(lpd_loo)[[2L]]
    w_bma <- loo::pseudobma_weights(lpd_loo)[[2L]]

    # The same folds for both, so the comparison is paired by observation.
    folds <- loo::kfold_split_random(K = 5L, N = n_train)
    lpd_kfold <- timed("kfold", cbind(
      loo::kfold(full$fit, folds = folds)$pointwise[, "elpd_kfold"],
      loo::kfold(flat$fit, folds = folds)$pointwise[, "elpd_kfold"]
    ))
    w_kstack <- loo::stacking_weights(lpd_kfold)[[2L]]
    w_kbma <- loo::pseudobma_weights(lpd_kfold)[[2L]]

    elpd_diff <- sum(lpd_loo[, 1L] - lpd_loo[, 2L])
    kfold_diff <- sum(lpd_kfold[, 1L] - lpd_kfold[, 2L])

    results <- list(full = full, flat = flat, component = component,
                    stack = average(full, flat, w_stack),
                    bma = average(full, flat, w_bma),
                    kstack = average(full, flat, w_kstack),
                    kbma = average(full, flat, w_kbma))
    w_flat <- c(full = 0, flat = 1, component = NA, stack = w_stack,
                bma = w_bma, kstack = w_kstack, kbma = w_kbma)

    out <- do.call(rbind, lapply(names(results), function(name) {
      data.frame(scenario = unit$scenario, rep = unit$rep, fit = name,
                 t(summarize(results[[name]], test, truth)),
                 w_flat = w_flat[[name]], elpd_diff = elpd_diff,
                 kfold_diff = kfold_diff,
                 seconds = if (name %in% names(seconds)) seconds[[name]]
                           else NA_real_)
    }))

    saveRDS(out, file.path(units_dir, sprintf("%s-%02d.rds", unit$scenario,
                                              unit$rep)))
    out
  }

  environment()
})

dir.create(study$units_dir, showWarnings = FALSE)

# The forms study's unit grid and seeds, in its order.
units <- expand.grid(rep = seq_len(20L), scenario = names(study$scenarios),
                     stringsAsFactors = FALSE)
units$seed <- 20261101L + seq_len(nrow(units))
units <- units[units$rep <= reps, ]
units <- split(units, seq_len(nrow(units)))

pr <- prog_init(total = length(units), title = "dpm_ls() scale switch",
                unit = "data set", workers = workers, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

future::plan(future::multisession, workers = workers)

results <- prog_future_lapply(pr, units, study$one_dataset,
                              label = function(u) paste(u$scenario, u$rep))

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(res = results, complete = TRUE), "_dev/scale-switch.rds")

on.exit()
prog_end(pr, "done")

means <- stats::aggregate(cbind(score, rmse, mean_cover, cover_lo, cover_hi,
                                w_flat) ~ scenario + fit,
                          results, mean, na.action = stats::na.pass)
print(means, digits = 3)

# Paired differences against the full fit, within data set.
for (scenario in unique(results$scenario)) {
  d <- results[results$scenario == scenario, ]
  base <- d[d$fit == "full", c("rep", "score", "rmse")]

  for (other in c("flat", "component", "stack", "bma", "kstack", "kbma")) {
    o <- merge(base, d[d$fit == other, c("rep", "score", "rmse")], by = "rep")
    ds <- o$score.y - o$score.x
    dr <- o$rmse.y - o$rmse.x
    cat(sprintf("%s: %s - full score %+.4f (SE %.4f), rmse %+.4f (SE %.4f)\n",
                scenario, other, mean(ds), stats::sd(ds) / sqrt(nrow(o)),
                mean(dr), stats::sd(dr) / sqrt(nrow(o))))
  }
}
