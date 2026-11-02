# Why does `dpm_ls()` lose to the paper's form when the scale is constant?
#
# What is tested. `_dev/dpm-ls-forms.R` found the whole-error form (`dpm_ls()`)
# behind `dpm_ls_component()` on its `homoskedastic` scenario (skewed W, s(x)
# fixed at .6): mean RMSE .090 against .080 (difference .010, SE .002) and
# `dpm()` at .079. The hypothesis is that the reported mean, f + exp(g) c,
# carries the posterior wiggle of g through c, the raw mixture center. With a
# skewed W the center is away from zero, and with a constant true scale the
# line along which the center trades against the level of f is nearly flat, so
# the center can sit far out. The component form's mean is f + c, which g does
# not touch. The claim could be false: the cost could be in the whole form's
# sampler generally (its scale forest is on the general path), in which case it
# would show with the scale forest emptied too.
#
# What is measured. The 20 homoskedastic data sets of the main study, on the
# same seeds, refit with: `dpm_ls()` and its full scale forest (`whole`), a
# repeat of the main study's fit; `dpm_ls()` with `log_sd = ~ 1` (`whole_flat`),
# where g cannot vary with x; `dpm_ls_component()` (`component`); and `dpm()`.
# Then 20 more data sets with the same design but normal W (`normal_homo`),
# fit with `whole` and `component`. Per fit:
#
#   rmse        RMSE of the posterior mean of E[y | x]. THE COLUMN TO READ.
#   score       mean held-out log predictive density.
#   center      posterior mean of |center| / spread, the raw center of W in its
#               own standard deviations (whole-error fits only).
#   g_spread    sd over the test points of the posterior mean of log_sd.
#
# What each outcome would mean.
#   - If `whole_flat` matches `dpm()` on rmse, `whole` stays behind both, and
#     under normal W the two forms tie (with `center` near zero), the cost is
#     the coupling of g with an off-center c, and the fix to look for keeps the
#     raw center near zero in the sampler.
#   - If `whole_flat` is behind `dpm()` as well, or normal W shows the gap too,
#     the cost is not that coupling, and the sampler or the prior of the
#     whole-error form is the thing to look at.
#
# Run with: BARTISAN_LIB=<library holding the working tree> \
#             Rscript _dev/dpm-ls-homoskedastic.R [reps] [workers]
# Writes:   _dev/dpm-ls-homoskedastic/<case>-<rep>.rds per data set, and
#           _dev/dpm-ls-homoskedastic.rds with all of them at the end.

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
  units_dir <- "_dev/dpm-ls-homoskedastic"

  sim_x <- function(n) {
    as.data.frame(matrix(stats::runif(5L * n), n, 5L,
                         dimnames = list(NULL, paste0("x", 1:5))))
  }

  mean_fun <- function(x) 1 + 0.8 * sin(pi * x$x1) + 0.6 * x$x2
  std_gamma <- function(n) (stats::rgamma(n, shape = 1.5) - 1.5) / sqrt(1.5)

  # The draw order matches `_dev/dpm-ls-forms.R`'s, so `skewed_homo` on its
  # seeds reproduces that study's `homoskedastic` data sets exactly.
  errors <- list(
    skewed_homo = function(x) 0.6 * std_gamma(nrow(x)),
    normal_homo = function(x) 0.6 * stats::rnorm(nrow(x))
  )

  predictors <- paste0("x", 1:5)
  mean_formula <- stats::reformulate(predictors, "y")
  scale_formula <- stats::reformulate(predictors)

  fits <- list(
    whole = list(family = function() dpm_ls(), scale = scale_formula),
    whole_flat = list(family = function() dpm_ls(), scale = ~ 1),
    component = list(family = function() bartisan:::dpm_ls_component(),
                     scale = scale_formula),
    dpm = list(family = function() dpm(), scale = NULL)
  )

  fit_one <- function(spec, train) {
    two <- !is.null(spec$scale)
    formula <- if (two) list(mean_formula, log_sd = spec$scale) else mean_formula
    bartisan(formula, data = train, family = spec$family(),
             num_trees = if (two) c(50L, 10L) else 50L, num_burn = 500L,
             num_draws = 500L, verbose = FALSE)
  }

  one_dataset <- function(unit) {
    set.seed(unit$seed)
    train <- sim_x(n_train)
    test <- sim_x(n_test)
    train$y <- mean_fun(train) + errors[[unit$case]](train)
    test$y <- mean_fun(test) + errors[[unit$case]](test)
    truth <- mean_fun(test)

    names <- if (unit$case == "skewed_homo") names(fits) else c("whole", "component")

    rows <- lapply(names, function(name) {
      fit <- fit_one(fits[[name]], train)
      mean_hat <- stats::predict(fit, newdata = test, type = "response")
      score <- mean(stats::predict(fit, newdata = test, type = "density",
                                   log = TRUE))
      aux <- fit[["aux"]]
      whole <- name %in% c("whole", "whole_flat")
      g_spread <- {
        if (name == "dpm") NA_real_
        else stats::sd(stats::predict(fit, newdata = test, type = "link")[, 2L])
      }

      data.frame(case = unit$case, rep = unit$rep, fit = name,
                 rmse = sqrt(mean((mean_hat - truth)^2)), score = score,
                 center = if (whole) mean(abs(aux[, "center"]) / aux[, "spread"])
                          else NA_real_,
                 g_spread = g_spread)
    })

    out <- do.call(rbind, rows)
    saveRDS(out, file.path(units_dir, sprintf("%s-%02d.rds", unit$case,
                                              unit$rep)))
    out
  }

  environment()
})

dir.create(study$units_dir, showWarnings = FALSE)

# The main study's units were numbered with `rep` fastest over its six
# scenarios in order, `homoskedastic` last, from seed 20261101 + 1; its data
# sets are therefore seeds 20261101 + 100 + rep. The normal ones take fresh
# seeds after the main study's 120.
units <- rbind(
  data.frame(case = "skewed_homo", rep = seq_len(reps),
             seed = 20261101L + 100L + seq_len(reps)),
  data.frame(case = "normal_homo", rep = seq_len(reps),
             seed = 20261101L + 200L + seq_len(reps))
)
units <- split(units, seq_len(nrow(units)))

pr <- prog_init(total = length(units), title = "dpm_ls() homoskedastic",
                unit = "data set", workers = workers, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

future::plan(future::multisession, workers = workers)

results <- prog_future_lapply(pr, units, study$one_dataset,
                              label = function(u) paste(u$case, u$rep))

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(res = results, complete = TRUE), "_dev/dpm-ls-homoskedastic.rds")

on.exit()
prog_end(pr, "done")

means <- stats::aggregate(cbind(rmse, score, center, g_spread) ~ case + fit,
                          results, mean, na.action = stats::na.pass)
print(means, digits = 3)

# Paired differences against the whole-error fit, within data set.
for (case in unique(results$case)) {
  d <- results[results$case == case, ]
  base <- d[d$fit == "whole", c("rep", "rmse")]

  for (other in setdiff(unique(d$fit), "whole")) {
    o <- merge(base, d[d$fit == other, c("rep", "rmse")], by = "rep")
    diff <- o$rmse.x - o$rmse.y
    cat(sprintf("%s: whole - %s rmse = %.4f (SE %.4f)\n", case, other,
                mean(diff), stats::sd(diff) / sqrt(length(diff))))
  }
}
