# `dpm_ls()` and `dpm_ls_component()` against nftbart, the reference
# implementation of Sparapani et al. (2023), with hard decision rules in all
# three.
#
# What is tested. The paper writes its error model (eq. 6, and eq. 7 for the
# censored draws) as e_i ~ N(mu_i, sigma_i^2 s(x_i)^2): only the component
# variances scale. `dpm_ls_component()` was written from that. But nftbart's
# code fits the other model. Its sampler (`cnft2.h` in version 1.5, from the
# paper's supplement) draws each censored time from N(f + s mu_c, s^2 sigma_c^2),
# subtracts s mu_c before the mean forest's update, and hands the mixture
# (y - f) / s; and the installed 2.3's `predict.nft2()` builds the density from
# components N(f + s mu_h, s^2 sigma_h^2). That is y = f + s W with W a DP
# mixture, the model `dpm_ls()` fits. Two claims follow, and either could be
# false:
#   (1) nftbart behaves like `dpm_ls()`, not like `dpm_ls_component()`, on the
#       scenarios where the two forms differ. If the reading of the code is
#       wrong, nftbart would instead lose to `dpm_ls()` on `bimodal` and
#       `skewed` and beat it on `component`, as the component form did in
#       `_dev/dpm-ls-forms.R`.
#   (2) `dpm_ls()` is no worse than nftbart, the published implementation of
#       the same model, in any scenario. The two differ in their priors (the
#       mean forest's leaf scale, k = 5 in nftbart and 2 here; nftbart's LIO
#       hyperpriors on the DP baseline, drawn, against a fixed conjugate
#       baseline here; nftbart's product-of-trees scale forest against a sum
#       on the log scale here), in nftbart's constraint step, which recenters
#       and rescales the mixture inside the chain, and in its scale update,
#       which holds the s mu_c term in the mean at the previous s.
#
# What is measured. The 120 data sets of `_dev/dpm-ls-forms.R`, reproduced
# exactly from its seeds: six error scenarios, 500 training and 1000 test
# points, five uniform predictors, scale driven by x3. Three fits per data set,
# all with hard rules, 50 trees on the mean and 10 on the scale, and the same
# chain length: nftbart at its defaults (1000 adaptation, 1000 burn-in, 2000
# kept), and `dpm_ls()` and `dpm_ls_component()` with `gate = "hard"`, 2000
# burn-in and 2000 kept. nftbart is given exp(y) as event times, all observed,
# so its log time is y. Per fit:
#
#   score       mean held-out log predictive density per test observation,
#               the log of the posterior mean density. THE COLUMN TO READ. For
#               nftbart it is computed from its own predictive formula, the
#               formula of `predict.nft2()`.
#   rmse        RMSE of the posterior mean of E[y | x] against the truth. For
#               nftbart, E[y | x] is the mean of its predictive mixture,
#               f + s sum_h w_h mu_h; `rmse_f` uses f alone, which the paper's
#               constraint means to be the same thing.
#   mean_cover  share of test points whose 90% interval for E[y | x] holds the
#               truth (nominal .90).
#   cover_lo, cover_hi  coverage of 90% posterior predictive intervals in the
#               lowest and highest third of the true conditional SD.
#   seconds     fitting time.
#
# Differences are taken within data set, `dpm_ls()` minus each other fit.
#
# What each outcome would mean.
#   - Claim (1) holds if nftbart minus `dpm_ls_component()` has the sign
#     `dpm_ls()` minus `dpm_ls_component()` had: higher score on `skewed` and
#     `bimodal`, lower on `component`. Then the paper's equations and its code
#     disagree, `dpm_ls()` is the model the authors ship, and the component
#     form is the one without an implementation behind it. If nftbart instead
#     tracks the component form, the reading of the code is wrong and has to
#     be redone before anything is said about it.
#   - Claim (2) holds if `dpm_ls()`'s score is within two SE of nftbart's or
#     above it in every scenario. If nftbart is ahead somewhere beyond two SE,
#     one of the differences listed above is worth having, and the next step
#     is to find which, starting with the mean forest's k and the DP baseline's
#     hyperpriors, which are the largest.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/nftbart-compare.R [reps] [workers]
# Writes:   _dev/nftbart-compare/<scenario>-<rep>.rds per data set, and
#           _dev/nftbart-compare.rds with all of them at the end.

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
  units_dir <- "_dev/nftbart-compare"
  num_burn <- 2000L
  num_draws <- 2000L

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

  summarize <- function(mu, mu_f, dens, reps, test, truth) {
    lo <- apply(reps, 2L, stats::quantile, 0.05)
    hi <- apply(reps, 2L, stats::quantile, 0.95)
    inside <- test$y >= lo & test$y <= hi
    breaks <- unique(stats::quantile(truth$sd, c(0, 1 / 3, 2 / 3, 1)))
    third <- {
      if (length(breaks) < 4L) rep(NA_integer_, length(truth$sd))
      else cut(truth$sd, breaks, include.lowest = TRUE, labels = FALSE)
    }
    mean_lo <- apply(mu, 2L, stats::quantile, 0.05)
    mean_hi <- apply(mu, 2L, stats::quantile, 0.95)

    c(score = mean(log(dens)),
      rmse = sqrt(mean((colMeans(mu) - truth$mean)^2)),
      rmse_f = sqrt(mean((colMeans(mu_f) - truth$mean)^2)),
      mean_cover = mean(truth$mean >= mean_lo & truth$mean <= mean_hi),
      cover_lo = mean(inside[which(third == 1L)]),
      cover_hi = mean(inside[which(third == 3L)]))
  }

  fit_bartisan <- function(family, train, test, truth) {
    fit <- bartisan::bartisan(list(mean_formula, log_sd = scale_formula), data = train,
                    family = family, num_trees = c(50L, 10L), gate = "hard",
                    num_burn = num_burn, num_draws = num_draws,
                    verbose = FALSE)
    mu <- stats::predict(fit, newdata = test, type = "response", draws = TRUE)
    dens <- stats::predict(fit, newdata = test, type = "density")
    reps <- bartisan:::posterior_predict.bartisan_fit(fit, newdata = test)
    summarize(mu, mu, dens, reps, test, truth)
  }

  # nftbart's predictive mixture at each test point and draw, as
  # `predict.nft2()` forms it: component h is N(f + s mu_h, (s sigma_h)^2) with
  # weight w_h, where f and s are its own predictions at the test points.
  fit_nftbart <- function(train, test, truth) {
    x <- as.matrix(train[predictors])
    xt <- as.matrix(test[predictors])
    utils::capture.output({
      fit <- nftbart::nft(x, exp(train$y), rep(1, nrow(x)), x.test = xt,
                          tc = 1L, printevery = 1e9)
      # `K = 0` falls back to the fit's grid of 100 event times and computes
      # survival curves at every test point; `K = NULL` gives f and s only.
      pred <- stats::predict(fit, xt, K = NULL, tc = 1L)
    })

    f <- pred$f.test
    s <- pred$s.test
    w <- fit$dpwt.
    m <- fit$dpmu.
    sd <- fit$dpsd.
    empty <- is.na(w) | w == 0
    w[empty] <- 0
    m[empty] <- 0
    sd[empty] <- 1

    shift <- rowSums(w * m)
    mu <- f + s * shift
    dens <- matrix(0, nrow(f), ncol(f))

    for (h in seq_len(ncol(w))) {
      dens <- dens + w[, h] * stats::dnorm(rep(test$y, each = nrow(f)),
                                           f + s * m[, h], s * sd[, h])
    }

    # One draw from the predictive at each test point and draw.
    pick <- t(apply(w, 1L, function(p) {
      sample.int(length(p), ncol(f), replace = TRUE, prob = p)
    }))
    rows <- matrix(seq_len(nrow(f)), nrow(f), ncol(f))
    reps <- f + s * (m[cbind(c(rows), c(pick))] +
                       sd[cbind(c(rows), c(pick))] *
                       stats::rnorm(length(f)))

    summarize(mu, f, colMeans(dens), reps, test, truth)
  }

  fits <- list(
    nftbart = function(train, test, truth) fit_nftbart(train, test, truth),
    whole = function(train, test, truth) {
      fit_bartisan(bartisan::dpm_ls(), train, test, truth)
    },
    component = function(train, test, truth) {
      fit_bartisan(bartisan:::dpm_ls_component(), train, test, truth)
    }
  )

  one_dataset <- function(unit) {
    set.seed(unit$seed)
    train <- sim_x(n_train)
    test <- sim_x(n_test)
    draw_train <- scenarios[[unit$scenario]](train)
    draw_test <- scenarios[[unit$scenario]](test)
    train$y <- mean_fun(train) + draw_train$e
    test$y <- mean_fun(test) + draw_test$e
    truth <- list(mean = mean_fun(test), sd = draw_test$sd)

    rows <- lapply(names(fits), function(name) {
      started <- Sys.time()
      res <- fits[[name]](train, test, truth)
      seconds <- as.numeric(Sys.time() - started, units = "secs")
      data.frame(scenario = unit$scenario, rep = unit$rep, fit = name,
                 t(res), seconds = seconds)
    })

    out <- do.call(rbind, rows)
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

pr <- prog_init(total = length(units), title = "dpm_ls() against nftbart",
                unit = "data set", workers = workers, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

future::plan(future::multisession, workers = workers)

results <- prog_future_lapply(pr, units, study$one_dataset,
                              label = function(u) paste(u$scenario, u$rep))

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(res = results, complete = TRUE), "_dev/nftbart-compare.rds")

on.exit()
prog_end(pr, "done")

means <- stats::aggregate(cbind(score, rmse, rmse_f, mean_cover, cover_lo,
                                cover_hi, seconds) ~ scenario + fit,
                          results, mean, na.action = stats::na.pass)
print(means, digits = 3)

for (scenario in unique(results$scenario)) {
  d <- results[results$scenario == scenario, ]
  base <- d[d$fit == "whole", c("rep", "score", "rmse")]

  for (other in c("nftbart", "component")) {
    o <- merge(base, d[d$fit == other, c("rep", "score", "rmse")], by = "rep")
    ds <- o$score.x - o$score.y
    dr <- o$rmse.x - o$rmse.y
    cat(sprintf("%s: whole - %s score %.4f (SE %.4f), rmse %.4f (SE %.4f)\n",
                scenario, other, mean(ds), stats::sd(ds) / sqrt(nrow(o)),
                mean(dr), stats::sd(dr) / sqrt(nrow(o))))
  }

  o <- merge(d[d$fit == "nftbart", c("rep", "score")],
             d[d$fit == "component", c("rep", "score")], by = "rep")
  ds <- o$score.x - o$score.y
  cat(sprintf("%s: nftbart - component score %.4f (SE %.4f)\n", scenario,
              mean(ds), stats::sd(ds) / sqrt(nrow(o))))
}
