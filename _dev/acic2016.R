# The 2016 Atlantic Causal Inference Competition, as a bench for configuration
# choices (Dorie, Hill, Shalit, Scott and Cervone 2019).
#
# Not run as part of anything. Invoke it deliberately:
#
#   Rscript _dev/acic2016.R smoke                 # two settings, one sim each
#   Rscript _dev/acic2016.R run 8 5               # 8 settings, 5 sims each
#   Rscript _dev/acic2016.R run all 10            # every setting, 10 sims each
#
# Writes _dev/acic2016-results.rds and prints the table.
#
# ---- what this bench is -----------------------------------------------------
#
# 4802 real covariates from the Collaborative Perinatal Project, 58 of them, and
# a simulated treatment and outcome. 77 parameter settings crossed from six
# factors, each with 100 replications, so 7700 datasets in all:
#
#   model.trt    linear | polynomial | step        how treatment depends on x
#   root.trt     0.35 | 0.65                       proportion treated
#   overlap.trt  full | one-term                   whether every unit could be treated
#   model.rsp    linear | exponential | step       how the outcome depends on x
#   alignment    0 | 0.25 | 0.75                   how far the two share predictors
#   te.hetero    none | med | high                 how much the effect varies
#
# `dgp_2016()` returns both potential outcomes and the true propensity score, so
# every estimand is known exactly and a fit can be scored against the sample it
# was given rather than against a population.
#
# The competition's estimand is the **sample** average effect on the treated,
# `mean((y.1 - y.0)[z == 1])`, and it scored bias, root mean squared error,
# interval coverage and interval length. Those are what this reports.
#
# ---- installing the data ----------------------------------------------------
#
# `remotes::install_github("vdorie/aciccomp/2016")` is the documented route and
# needs the GitHub API. Where that is blocked, the tarball works:
#
#   curl -sSL -o acic.tgz \
#     https://codeload.github.com/vdorie/aciccomp/tar.gz/refs/heads/master
#   tar xzf acic.tgz
#   R CMD INSTALL aciccomp-master/2016
#
# See _dev/SIMULATION.md for the design this shares with the other benches.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0L) args[1L] else "none"

if (identical(mode, "none")) {
  cat("Nothing run. See the header of this file for how to invoke it.\n")
  quit(save = "no")
}

if (!requireNamespace("aciccomp2016", quietly = TRUE)) {
  stop("aciccomp2016 is not installed; see the header of this file")
}

# ---- the configurations under test ------------------------------------------
#
# Anchors first, for the reason in _dev/SIMULATION.md: `oracle` is handed the
# true propensity score and `none` gets none at all, and between them is the
# range any choice of propensity settings is competing over.

PS <- list(num_burn = 200L, num_draws = 400L, chains = 2L)

# `both` is what the mixing work recommends and `fixed_scale` half of it; the
# question here is whether either costs anything in bias where a propensity
# score is actually load-bearing, which is what poor overlap and low alignment
# are for.

configurations <- list(
  none        = list(propensity = FALSE),
  oracle      = list(propensity = "true"),
  default     = list(ps = list()),
  fixed_scale = list(ps = list(update_sigma_mu = FALSE)),
  both        = list(ps = list(update_sigma_mu = FALSE, sparsity = FALSE))
)

OUT <- bartisan_control(num_burn = 200L, num_draws = 400L, chains = 2L)

fit_one <- function(spec, d, covs) {
  form <- stats::reformulate(covs, response = "y")

  propensity <- {
    if (isFALSE(spec[["propensity"]])) FALSE
    else if (identical(spec[["propensity"]], "true")) d[["e"]]
    else TRUE
  }

  ps_args <- {
    if (is.null(spec[["ps"]])) list()
    else list(control = do.call(bartisan_control, c(PS, spec[["ps"]])))
  }

  bcf(form, treatment = ~ z, data = d, family = gaussian(),
      propensity = propensity, propensity_args = ps_args, control = OUT)
}

# The effect forest carries the conditional effect directly, so the estimand is
# read off it rather than recomputed by differencing predictions. It has to come
# through `coef(draws = TRUE)` and not from the forest's own values: a
# varying-coefficient model splits the level between the control function and
# the coefficient, and where that split falls moves from draw to draw, which
# inflates an interval built on the raw forest without touching the mean.
score <- function(fit, d) {
  slopes <- stats::coef(fit, draws = TRUE)[[1L]]
  tau_hat <- colMeans(slopes)

  treated <- d[["z"]] == 1L
  satt_draws <- rowMeans(slopes[, treated, drop = FALSE])
  truth <- mean((d[["y.1"]] - d[["y.0"]])[treated])
  ci <- stats::quantile(satt_draws, c(0.025, 0.975), names = FALSE)

  c(satt = mean(satt_draws),
    truth = truth,
    bias = mean(satt_draws) - truth,
    covered = as.numeric(ci[1L] <= truth && truth <= ci[2L]),
    width = ci[2L] - ci[1L],
    pehe = sqrt(mean((tau_hat - (d[["mu.1"]] - d[["mu.0"]]))^2)))
}

# ---- the run ----------------------------------------------------------------

x <- aciccomp2016::input_2016
covs <- names(x)
grid <- aciccomp2016::parameters_2016

spec <- {
  if (identical(mode, "smoke")) list(settings = c(1L, 40L), sims = 1L)
  else {
    n_set <- if (length(args) > 1L) args[2L] else "8"
    sims <- if (length(args) > 2L) as.integer(args[3L]) else 5L

    settings <- {
      if (identical(n_set, "all")) seq_len(nrow(grid))
      else if (grepl(",", n_set) || grepl("^[0-9]+$", n_set) &&
                 as.integer(n_set) > nrow(grid)) {
        # Named rows, for aiming at the factors a question turns on rather than
        # sampling the grid blind.
        as.integer(strsplit(n_set, ",", fixed = TRUE)[[1L]])
      }
      else {
        # A spread across the grid rather than the first few, so the six factors
        # are all represented even at a small budget.
        k <- as.integer(n_set)
        unique(as.integer(round(seq(1, nrow(grid), length.out = k))))
      }
    }

    list(settings = settings, sims = sims)
  }
}

n_total <- length(spec$settings) * spec$sims * length(configurations)

cat(sprintf("%d settings x %d sims x %d configurations = %d fits\n\n",
            length(spec$settings), spec$sims, length(configurations),
            n_total))

tag <- if (identical(mode, "smoke")) {
  "smoke"
} else {
  paste(range(spec$settings), collapse = "-")
}

OUT <- sprintf("_dev/acic2016-results-%s.rds", tag)

out <- list()
done <- 0L

pr <- prog_init(total = n_total, title = "ACIC 2016 configurations",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the last dataset"), add = TRUE)

# Written after every dataset rather than after the last one, so a run that is
# killed leaves its finished datasets readable. `complete` is what tells a
# reader which of the two it is looking at.
checkpoint <- function(complete) {
  saveRDS(list(res = do.call(rbind, out), complete = complete, done = done,
               total = n_total, settings = spec$settings, sims = spec$sims),
          OUT)
}

for (p in spec$settings) {
  for (s in seq_len(spec$sims)) {
    d <- aciccomp2016::dgp_2016(x, p, s)
    d <- cbind(x, d)

    for (nm in names(configurations)) {
      # Paired: the same dataset and the same stream for every configuration.
      set.seed(97)
      elapsed <- system.time(
        fit <- fit_one(configurations[[nm]], d, covs))[["elapsed"]]

      out[[length(out) + 1L]] <- data.frame(
        as.list(score(fit, d)), configuration = nm, parameter = p, sim = s,
        seconds = elapsed, grid[p, , drop = FALSE], row.names = NULL)

      done <- done + 1L
      prog_tick(pr, i = done, secs = elapsed,
                label = sprintf("setting %d sim %d / %s", p, s, nm))
    }

    checkpoint(FALSE)
  }
}

checkpoint(TRUE)
results <- do.call(rbind, out)

on.exit()
prog_end(pr, "done", sprintf("%d fits", done))

# ---- the report -------------------------------------------------------------
#
# Scaled by the outcome's own spread the way the competition reported it, so a
# number is comparable across settings whose effects are of different size.

report <- function(results) {
  cat(sprintf("\n=== over %d settings, %d datasets ===\n",
              length(unique(results$parameter)), nrow(results) / length(unique(results$configuration))))
  cat(sprintf("%-12s %9s %9s %8s %9s %9s %8s\n", "config", "bias", "RMSE",
              "cover", "width", "PEHE", "seconds"))

  for (nm in unique(results$configuration)) {
    r <- results[results$configuration == nm, ]
    cat(sprintf("%-12s %+9.4f %9.4f %8.2f %9.3f %9.4f %8.1f\n", nm,
                mean(r$bias), sqrt(mean(r$bias^2)), mean(r$covered),
                mean(r$width), mean(r$pehe), mean(r$seconds)))
  }

  # Where a configuration wins and loses is the point of the 77 settings, so the
  # factors get their own breakdown rather than being averaged away.
  for (f in c("overlap.trt", "te.hetero", "model.rsp")) {
    cat(sprintf("\n  RMSE by %s:\n", f))
    lv <- unique(results[[f]])
    cat(sprintf("  %-12s %s\n", "config",
                paste(sprintf("%10s", as.character(lv)), collapse = "")))

    for (nm in unique(results$configuration)) {
      r <- results[results$configuration == nm, ]
      cells <- vapply(lv, function(l)
        sqrt(mean(r$bias[r[[f]] == l]^2)), numeric(1L))
      cat(sprintf("  %-12s %s\n", nm,
                  paste(sprintf("%10.4f", cells), collapse = "")))
    }
  }
}

report(results)
