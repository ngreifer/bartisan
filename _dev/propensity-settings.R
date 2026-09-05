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

library(bartisan)

reps <- {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) > 0L) as.integer(a[1L]) else 25L
}

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
                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  m <- design(x)
  p <- ncol(m)
  strong <- seq_len(min(6L, p))

  lin <- function(cols, w) as.vector(m[, cols, drop = FALSE] %*% w)

  w_t <- stats::rnorm(length(strong))
  w_y <- stats::rnorm(length(strong))

  et <- {
    if (identical(selection, "linear")) lin(strong, w_t)
    else lin(strong, w_t) + 1.5 * m[, strong[1L]] * m[, strong[2L]] -
           1.5 * abs(m[, strong[3L]])
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
  tau_hat <- as.vector(stats::coef(fit)[, 1L])
  draws <- fit[["eta"]][[2L]]
  draws <- draws - mean(draws) + mean(tau_hat)

  treated <- d[["z"]] == 1L
  att_draws <- rowMeans(draws[, treated, drop = FALSE])
  truth <- mean(d[["tau"]][treated])
  ci <- stats::quantile(att_draws, c(0.025, 0.975), names = FALSE)

  c(att = mean(att_draws),
    bias = mean(att_draws) - truth,
    covered = as.numeric(ci[1L] <= truth && truth <= ci[2L]),
    width = ci[2L] - ci[1L],
    cate_rmse = sqrt(mean((tau_hat - d[["tau"]])^2)))
}

# ---- the run ----------------------------------------------------------------

designs <- expand.grid(
  data = c("rhc", "lalonde"),
  selection = c("linear", "nonlinear"),
  surface = c("nonlinear"),
  hetero = c(0.5),
  stringsAsFactors = FALSE
)

out <- list()

for (di in seq_len(nrow(designs))) {
  dg <- designs[di, ]
  x <- covariates(dg[["data"]])
  covs <- names(x)

  for (r in seq_len(reps)) {
    d <- simulate(x, dg[["selection"]], dg[["surface"]], dg[["hetero"]],
                  seed = 10000 * di + r)

    for (nm in names(settings)) {
      # The same data and the same stream for every setting, so the comparison
      # is paired and the differences carry far less noise than the levels.
      set.seed(97)
      elapsed <- system.time(fit <- fit_one(settings[[nm]], d, covs))[["elapsed"]]

      out[[length(out) + 1L]] <- data.frame(
        as.list(score(fit, d)), setting = nm, rep = r, seconds = elapsed,
        data = dg[["data"]], selection = dg[["selection"]],
        surface = dg[["surface"]], hetero = dg[["hetero"]]
      )
    }

    cat(sprintf("  %s / %s: rep %d of %d\n", dg[["data"]], dg[["selection"]],
                r, reps))
    utils::flush.console()
  }
}

results <- do.call(rbind, out)
saveRDS(results, "_dev/propensity-settings.rds")

# ---- the report -------------------------------------------------------------

report <- function(results) {
  for (dn in unique(results[["data"]])) for (sn in unique(results[["selection"]])) {
    r <- results[results[["data"]] == dn & results[["selection"]] == sn, ]
    if (nrow(r) == 0L) next

    cat(sprintf("\n=== %s, %s selection (%d replicates) ===\n", dn, sn,
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
