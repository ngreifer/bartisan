# Does `dpm_ls()`, which scales the whole error, do at least as well as the
# heteroskedastic mixture of Sparapani et al. (2023), which scales only each
# component's standard deviation?
#
# What is tested. `dpm_ls()` fits y = f(x) + exp(g(x)) W with W a Dirichlet
# process mixture, so the error keeps one shape and only its size moves with x.
# The paper's form, kept unexported as `dpm_ls_component()`, fits
# e ~ N(mu_c, sigma_c^2 exp(2 g(x))): the component means do not move with x.
# The claim is that the first does no worse than the second anywhere and better
# where the shape is not normal and the size varies. It could be false in two
# ways. The paper's form can represent a shape that changes with x, which the
# whole-error form cannot, so it should win when the truth is its own model
# (scenario `component`). And the whole-error form reports its conditional
# mean as f + exp(g) * center, which carries the noise in g, and the one smoke
# run before this (skewed errors, n = 400) had its mean RMSE at 0.14 against the
# paper form's 0.11; if that holds up it is a cost on the quantity most people
# read.
#
# What is measured. Six error scenarios on the same design: 500 training and
# 1000 test observations, five uniform predictors, E[y | x] = 1 + 0.8 sin(pi x1)
# + 0.6 x2, and a conditional standard deviation s(x) = exp(-1.2 + 1.2 x3) that
# depends on a predictor the mean does not, so each forest has its own signal.
# W has mean 0 and variance 1 in every scenario but `component`, so s(x) is the
# conditional SD and f the conditional mean.
#
#   normal       W normal: both forms are right, and both should tie.
#   skewed       W a standardized gamma(1.5): shape fixed, size varies.
#   heavy        W a standardized t_3.
#   bimodal      W = +-0.9 + N(0, 0.19): two modes that move apart as s grows.
#   component    the paper's model as truth: e = +-0.6 + 0.6 s(x) z, so the
#                modes stay at +-0.6 and only their widths grow. Its
#                conditional SD is 0.6 sqrt(1 + s(x)^2).
#   homoskedastic  W skewed as above with s(x) fixed at 0.6: `dpm()`'s home
#                ground, where the scale forest has nothing to find.
#
# Four fits per data set, all 500 warmup + 500 draws, one chain, 50 trees on the
# mean and 10 on the scale: `dpm_ls()`, `dpm_ls_component()`, and `gaussian_ls()`
# and `dpm()` for reference. Twenty data sets per scenario, paired across fits.
# Per fit, on the test set:
#
#   score      mean log predictive density per test observation. THE COLUMN TO
#              READ: it scores the whole conditional distribution, which is what
#              the two forms differ in.
#   rmse       RMSE of the posterior mean of E[y | x] against the truth.
#   mean_cover share of test points whose 90% interval for E[y | x] holds the
#              truth (nominal .90).
#   cover_lo, cover_hi  coverage of 90% posterior predictive intervals among the
#              test points in the lowest and highest third of the true
#              conditional SD (nominal .90 each). A form that compresses the
#              spread over-covers the first and under-covers the second.
#   seconds    fitting time.
#
# Differences are taken within data set, whole minus component, and summarized
# by their mean and standard error over the twenty.
#
# What each outcome would mean.
#   - If the whole-error form's score is at or above the component form's in
#     every scenario but `component`, within two standard errors, and its rmse
#     is not worse by more than two standard errors anywhere but `component`,
#     the claim holds: `dpm_ls()` ships as written and the paper's form stays
#     internal, for reference.
#   - If it loses only on `component`, that is the one model it cannot
#     represent, and it is recorded as the whole-error form's limit rather than
#     as a reason to change it; the claim is then reported as false there.
#   - If its rmse is worse beyond two standard errors on scenarios where it
#     scores better, the mean it reports is paying for the noise in g, and the
#     reported mean has to be revisited before shipping, whatever the score
#     says.
#   - If it loses on `normal` or `homoskedastic`, where it is as right as the
#     other form, the cost is in the sampler (the scale forest is on the
#     general path), and that is the thing to look at.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/dpm-ls-forms.R [reps] [workers]
# Writes:   _dev/dpm-ls-forms/<scenario>-<rep>.rds per data set, as each
#           finishes, and _dev/dpm-ls-forms.rds with all of them at the end.

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

# Everything a worker needs lives in this one environment, which travels with
# `one_dataset()` when it is sent to a worker. Defined at the top level instead,
# the scenario generators' own helpers were not found there, since *future*
# looks for globals in the function it sends and not in the closures it calls.
study <- local({
  n_train <- 500L
  n_test <- 1000L
  units_dir <- "_dev/dpm-ls-forms"

  sim_x <- function(n) {
    as.data.frame(matrix(stats::runif(5L * n), n, 5L,
                         dimnames = list(NULL, paste0("x", 1:5))))
  }

  mean_fun <- function(x) 1 + 0.8 * sin(pi * x$x1) + 0.6 * x$x2
  sd_fun <- function(x) exp(-1.2 + 1.2 * x$x3)

  std_gamma <- function(n) (stats::rgamma(n, shape = 1.5) - 1.5) / sqrt(1.5)

  # Each returns the error and the true conditional SD at each x.
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

  families <- list(
    whole = function() dpm_ls(),
    component = function() bartisan:::dpm_ls_component(),
    gaussian_ls = function() gaussian_ls(),
    dpm = function() dpm()
  )

  fit_one <- function(name, train) {
    two <- name != "dpm"
    formula <- if (two) list(mean_formula, log_sd = scale_formula) else mean_formula
    bartisan(formula, data = train, family = families[[name]](),
             num_trees = if (two) c(50L, 10L) else 50L, num_burn = 500L,
             num_draws = 500L, verbose = FALSE)
  }

  score_fit <- function(fit, test, truth) {
    mu <- predict(fit, newdata = test, type = "response", draws = TRUE)
    mean_hat <- colMeans(mu)
    mean_lo <- apply(mu, 2L, stats::quantile, 0.05)
    mean_hi <- apply(mu, 2L, stats::quantile, 0.95)

    dens <- predict(fit, newdata = test, type = "density", log = TRUE)

    reps <- bartisan:::posterior_predict.bartisan_fit(fit, newdata = test)
    lo <- apply(reps, 2L, stats::quantile, 0.05)
    hi <- apply(reps, 2L, stats::quantile, 0.95)
    inside <- test$y >= lo & test$y <= hi
    # Thirds of a constant SD are not defined, and `homoskedastic` reports no
    # coverage by spread.
    breaks <- unique(stats::quantile(truth$sd, c(0, 1 / 3, 2 / 3, 1)))
    third <- {
      if (length(breaks) < 4L) rep(NA_integer_, length(truth$sd))
      else cut(truth$sd, breaks, include.lowest = TRUE, labels = FALSE)
    }

    c(score = mean(dens),
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

    rows <- lapply(names(families), function(name) {
      started <- Sys.time()
      fit <- fit_one(name, train)
      seconds <- as.numeric(Sys.time() - started, units = "secs")
      data.frame(scenario = unit$scenario, rep = unit$rep, fit = name,
                 t(score_fit(fit, test, truth)), seconds = seconds)
    })

    out <- do.call(rbind, rows)
    saveRDS(out, file.path(units_dir, sprintf("%s-%02d.rds", unit$scenario,
                                              unit$rep)))
    out
  }

  environment()
})

dir.create(study$units_dir, showWarnings = FALSE)

units <- expand.grid(rep = seq_len(reps), scenario = names(study$scenarios),
                     stringsAsFactors = FALSE)
units$seed <- 20261101L + seq_len(nrow(units))
units <- split(units, seq_len(nrow(units)))

pr <- prog_init(total = length(units), title = "dpm_ls() forms",
                unit = "data set", workers = workers, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

future::plan(future::multisession, workers = workers)

results <- prog_future_lapply(pr, units, study$one_dataset,
                              label = function(u) paste(u$scenario, u$rep))

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(res = results, complete = TRUE), "_dev/dpm-ls-forms.rds")

on.exit()
prog_end(pr, "done")

# The summary the header asks for: whole minus component within data set.
wide <- reshape(results[results$fit %in% c("whole", "component"),
                        c("scenario", "rep", "fit", "score", "rmse")],
                idvar = c("scenario", "rep"), timevar = "fit",
                direction = "wide")
paired <- do.call(rbind, lapply(split(wide, wide$scenario), function(d) {
  ds <- d$score.whole - d$score.component
  dr <- d$rmse.whole - d$rmse.component
  data.frame(scenario = d$scenario[1L], n = nrow(d),
             score_diff = mean(ds), score_se = stats::sd(ds) / sqrt(nrow(d)),
             rmse_diff = mean(dr), rmse_se = stats::sd(dr) / sqrt(nrow(d)))
}))
print(paired, digits = 3)

means <- stats::aggregate(cbind(score, rmse, mean_cover, cover_lo, cover_hi,
                                seconds) ~ scenario + fit, results, mean,
                          na.action = stats::na.pass)
print(means[order(means$scenario, means$fit), ], digits = 3)
