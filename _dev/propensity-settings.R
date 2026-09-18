# Which settings a propensity model should be fitted with, measured on the
# estimand rather than on the propensity score.
#
# The propensity score in `bcf()` is a nuisance: it goes into the control
# function so that the control function can absorb selection. So the question is
# not which settings predict treatment best, it is which settings leave the
# least confounding in the effect. Those are different questions, and the
# measurement has to be of the second.
#
# Run with: Rscript _dev/propensity-settings.R [reps]
# Writes:   _dev/propensity-settings.rds
#
# See _dev/SIMULATION.md for what makes this design work and how to reuse it.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

# `Rscript _dev/propensity-settings.R <reps> [design]`. Naming a design runs
# only that row of the grid below and writes its own file, so the four can go to
# the queue at once. Their `seconds` columns are then not comparable across
# jobs, which is why nothing is concluded from them; see `_dev/SIMULATION.md`.
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) > 0L) as.integer(args[1L]) else 25L
only <- if (length(args) > 1L) as.integer(args[2L]) else NA_integer_

# ---- the data ---------------------------------------------------------------
#
# Real covariates, simulated treatment and outcome. Real covariates keep the
# correlations, the skew and the discreteness a synthetic design smooths away;
# simulating the treatment and the outcome is what supplies a truth to score
# against. Two datasets of different shape, so a result that only holds at one
# aspect ratio shows itself.

covariates <- function(which) {
  if (identical(which, "rhc")) {
    data("rhc", package = "bartisan", envir = environment())
    keep <- c("age","sex","race","edu","aps","meanbp","resp","hema","pafi",
              "paco2","crea","surv2m","card")
    return(rhc[, keep])
  }

  data("lalonde", package = "cobalt", envir = environment())
  lalonde[, c("age","educ","race","married","nodegree","re74","re75")]
}

# A numeric design matrix, standardized, for building the true surfaces.
design <- function(x) {
  m <- stats::model.matrix(~ . - 1, data = x)
  m <- m[, apply(m, 2L, stats::sd) > 0, drop = FALSE]
  scale(m)
}

# The truth. `selection` says how the treatment depends on the covariates and
# `surface` how the outcome does; `hetero` scales the moderation. The treatment
# is centered to a target prevalence and the signals are scaled to a target
# strength, so that "linear" and "nonlinear" differ in shape and not in how much
# there is to find.
simulate <- function(x, selection, surface, hetero, prevalence = 0.4,
                     confounding = "aligned", noise = 0L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (noise > 0L) {
    extra <- matrix(stats::rnorm(nrow(x) * noise), nrow(x), noise)
    colnames(extra) <- sprintf("n%d", seq_len(noise))
    x <- cbind(x, as.data.frame(extra))
  }

  m <- design(x)
  p <- ncol(m)
  strong <- seq_len(min(6L, p))

  lin <- function(cols, w) as.vector(m[, cols, drop = FALSE] %*% w)

  w_t <- stats::rnorm(length(strong))
  w_y <- stats::rnorm(length(strong))

  # Under `targeted`, the covariates that drive treatment are ones the outcome
  # barely depends on. That is the case a propensity score exists for: the
  # outcome model shrinks a weak predictor away, the selection it carried goes
  # with it, and the bias that leaves is what the score in the control function
  # restores. With both surfaces on the same strong covariates, which is
  # `aligned`, the outcome model absorbs selection by itself and the score has
  # nothing to add.
  drivers <- {
    if (identical(confounding, "aligned")) strong
    else utils::tail(seq_len(p), min(3L, p))
  }

  w_t <- if (identical(confounding, "aligned")) w_t else stats::rnorm(length(drivers))

  et <- {
    if (identical(selection, "linear")) lin(drivers, w_t)
    else lin(drivers, w_t) + 1.5 * m[, drivers[1L]] * m[, drivers[2L]] -
           1.5 * abs(m[, drivers[min(3L, length(drivers))]])
  }
  et <- 1.2 * as.vector(scale(et))

  # Solved rather than guessed, so prevalence is what was asked for whatever the
  # shape of the selection.
  shift <- stats::uniroot(function(a) mean(stats::plogis(a + et)) - prevalence,
                          c(-20, 20))$root
  e <- stats::plogis(shift + et)
  z <- stats::rbinom(nrow(m), 1L, e)

  mu <- {
    if (identical(surface, "linear")) lin(strong, w_y)
    else lin(strong, w_y) + 2 * sin(pi * m[, strong[1L]]) +
           1.5 * m[, strong[2L]]^2
  }
  mu <- 2 * as.vector(scale(mu))

  # The drivers keep a weak hand in the outcome, so they are genuine confounders
  # rather than instruments; an instrument would bias the effect rather than
  # help it, and the propensity score would have nothing to restore.
  if (!identical(confounding, "aligned")) {
    mu <- mu + 0.35 * as.vector(scale(lin(drivers, rep.int(1, length(drivers)))))
  }

  tau <- 1 + hetero * as.vector(scale(m[, strong[2L]] + 0.5 * m[, strong[4L]]))
  y <- mu + z * tau + stats::rnorm(nrow(m))

  data.frame(x, z = z, y = y, e = e, tau = tau)
}

# ---- the settings under test ------------------------------------------------
#
# Two anchors bracket the comparison. `none` fits no propensity model at all and
# says how much the score is worth; `oracle` hands the fit the true score and
# says how much better any estimate of it could be. Anything between them is the
# range the settings are competing over, and a difference that is small next to
# that range is small whatever its standard error.

settings <- list(
  none        = list(propensity = FALSE),
  oracle      = list(propensity = "true"),
  default     = list(args = list()),
  no_sparsity = list(args = list(sparsity = FALSE)),
  fixed_scale = list(args = list(update_sigma_mu = FALSE)),
  both        = list(args = list(update_sigma_mu = FALSE, sparsity = FALSE)),
  undersmooth = list(args = list(k = 1, sparsity = FALSE)),
  under_fixed = list(args = list(k = 1, sparsity = FALSE,
                                 update_sigma_mu = FALSE))
)

PS_CONTROL <- list(num_burn = 200L, num_draws = 400L, chains = 2L)
OUT_CONTROL <- bartisan_control(num_burn = 200L, num_draws = 400L, chains = 2L)

fit_one <- function(spec, d, covs) {
  form <- stats::reformulate(covs, response = "y")

  propensity <- {
    if (isFALSE(spec[["propensity"]])) FALSE
    else if (identical(spec[["propensity"]], "true")) d[["e"]]
    else TRUE
  }

  args <- {
    if (is.null(spec[["args"]])) list()
    else list(control = do.call(bartisan_control,
                                c(PS_CONTROL, spec[["args"]])))
  }

  bcf(form, treatment = ~ z, data = d, family = gaussian(),
      propensity = propensity, propensity_args = args, control = OUT_CONTROL)
}

# ---- scoring ----------------------------------------------------------------
#
# The estimand is the sample average effect on the treated, which is what the
# effect forest reports per observation, so it is read off `coef()` rather than
# recomputed by prediction. The interval is over the posterior of that average.

score <- function(fit, d) {
  # `coef(draws = TRUE)` is the identified coefficient per draw. The forest's
  # own values are not: a varying-coefficient model splits the level between the
  # control function and the coefficient, and where that split falls moves from
  # draw to draw. Averaging the raw forest carries that movement into the
  # interval and made coverage 0.998 at a width five times the sampling spread
  # of the estimate. `vc_recenter()`, which `coef()` applies and the raw forest
  # does not, is the difference.
  slopes <- stats::coef(fit, draws = TRUE)[[1L]]
  tau_hat <- colMeans(slopes)

  treated <- d[["z"]] == 1L
  att_draws <- rowMeans(slopes[, treated, drop = FALSE])
  truth <- mean(d[["tau"]][treated])
  ci <- stats::quantile(att_draws, c(0.025, 0.975), names = FALSE)

  c(att = mean(att_draws),
    bias = mean(att_draws) - truth,
    covered = as.numeric(ci[1L] <= truth && truth <= ci[2L]),
    width = ci[2L] - ci[1L],
    cate_rmse = sqrt(mean((tau_hat - d[["tau"]])^2)))
}

# ---- the run ----------------------------------------------------------------

designs <- rbind(
  expand.grid(data = c("rhc", "lalonde"), selection = c("linear", "nonlinear"),
              surface = "nonlinear", hetero = 0.5, confounding = "aligned",
              noise = 0L, stringsAsFactors = FALSE),
  expand.grid(data = c("rhc", "lalonde"), selection = "nonlinear",
              surface = "nonlinear", hetero = 0.5, confounding = "targeted",
              noise = 20L, stringsAsFactors = FALSE)
)

out <- list()

rows <- if (is.na(only)) seq_len(nrow(designs)) else only
OUT <- if (is.na(only)) {
  "_dev/propensity-settings.rds"
} else {
  sprintf("_dev/propensity-settings-%d.rds", only)
}

n_total <- length(rows) * reps * length(settings)
done <- 0L

title <- if (is.na(only)) {
  "Propensity settings"
} else {
  sprintf("Propensity settings: design %d", only)
}

pr <- prog_init(total = n_total, unit = "fit", kind = "simulation",
                title = title)
on.exit(prog_end(pr, "failed", "aborted before the last replicate"),
        add = TRUE)

# Written after every replicate rather than after the last one, so a run that
# is killed leaves its finished replicates readable. `complete` is what tells a
# reader which of the two it is looking at.
checkpoint <- function(complete) {
  saveRDS(list(res = do.call(rbind, out), complete = complete, done = done,
               total = n_total, reps = reps, designs = designs[rows, ]), OUT)
}

for (di in rows) {
  dg <- designs[di, ]
  x <- covariates(dg[["data"]])
  covs <- names(x)

  for (r in seq_len(reps)) {
    d <- simulate(x, dg[["selection"]], dg[["surface"]], dg[["hetero"]],
                  confounding = dg[["confounding"]], noise = dg[["noise"]],
                  seed = 10000 * di + r)
    covs <- setdiff(names(d), c("z", "y", "e", "tau"))

    for (nm in names(settings)) {
      # The same data and the same stream for every setting, so the comparison
      # is paired and the differences carry far less noise than the levels.
      set.seed(97)
      elapsed <- system.time(fit <- fit_one(settings[[nm]], d, covs))[["elapsed"]]

      out[[length(out) + 1L]] <- data.frame(
        as.list(score(fit, d)), setting = nm, rep = r, seconds = elapsed,
        data = dg[["data"]], selection = dg[["selection"]],
        surface = dg[["surface"]], hetero = dg[["hetero"]],
        confounding = dg[["confounding"]]
      )

      done <- done + 1L
      prog_tick(pr, i = done, secs = elapsed,
                label = sprintf("%s / %s / %s / %s rep %d", dg[["data"]],
                                dg[["selection"]], dg[["confounding"]], nm, r))
    }

    checkpoint(FALSE)
  }
}

checkpoint(TRUE)
results <- do.call(rbind, out)

on.exit()
prog_end(pr, "done", sprintf("%d fits", done))

# ---- the report -------------------------------------------------------------

report <- function(results) {
  results[["design"]] <- paste(results[["data"]], results[["selection"]],
                               results[["confounding"]])

  for (dn in unique(results[["design"]])) {
    r <- results[results[["design"]] == dn, ]
    if (nrow(r) == 0L) next

    cat(sprintf("\n=== %s (%d replicates) ===\n", dn,
                length(unique(r[["rep"]]))))
    cat(sprintf("%-12s %8s %8s %8s %8s %9s %8s\n", "setting", "bias", "RMSE",
                "cover", "width", "cateRMSE", "seconds"))

    base <- r[r[["setting"]] == "default", ]
    base <- base[order(base[["rep"]]), "bias"]

    for (nm in unique(r[["setting"]])) {
      s <- r[r[["setting"]] == nm, ]
      s <- s[order(s[["rep"]]), ]
      cat(sprintf("%-12s %+8.4f %8.4f %8.2f %8.3f %9.4f %8.2f\n", nm,
                  mean(s[["bias"]]), sqrt(mean(s[["bias"]]^2)),
                  mean(s[["covered"]]), mean(s[["width"]]),
                  mean(s[["cate_rmse"]]), mean(s[["seconds"]])))
    }
  }
}

report(results)
