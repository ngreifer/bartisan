# BART as the outcome nuisance for difference-in-differences, in the
# Sant'Anna & Zhao (2020) / Callaway & Sant'Anna (2021) framework rather than
# the DiD-BCF one.
#
# With the IPW weights set to 1 this is the outcome-regression estimator of
# Heckman, Ichimura & Todd (1997), which `did` calls est_method = "reg". The
# nuisance is E[dY | X, D = 0]: a regression of the CHANGE in outcome on
# baseline covariates, fitted on the comparison group only. Nothing about it
# has to be linear, so BART goes straight in.
#
# Measured on `did::mpdta`, all with xformla = ~ lpop:
#
#   2x2 (cohort 2007 vs never, 2006 -> 2007)
#     BART OR, conditional on treated units   -0.0277 [-0.0454, -0.0111] w 0.034
#     BART OR, + Bayesian bootstrap           -0.0279 [-0.0597, +0.0039] w 0.064
#     DRDID::reg_did_panel (linear nuisance)  -0.0288 [-0.0605, +0.0029] w 0.063
#
#   staggered, via did::att_gt(est_method = )
#     est_method = "reg"                      -0.0420  se 0.0114
#     est_method = "dr"                       -0.0418  se 0.0115
#     est_method = bart_or                    -0.0421  se 0.0100
#
# Point estimates agree with the linear versions to the third decimal. The
# standard errors do not, and that is the part to be careful about: dropping
# the propensity term costs Neyman orthogonality, so error in the fitted
# nuisance enters the estimator's first-order behavior. The plug-in influence
# function below omits that term and reads about 13% small; the Bayesian
# posterior conditional on the treated units is narrower still, because it
# holds their realized outcome changes fixed, and a Bayesian bootstrap over the
# treated units restores the population-ATT width to within 0.3%.
#
# To keep the orthogonality, supply both nuisances -- a `binomial()` BART for
# the propensity score alongside this one -- and cross-fit. That is the
# double/debiased ML route and `did`'s est_method hook takes it.

suppressPackageStartupMessages(library(bartisan))
data("mpdta", package = "did"); d <- as.data.frame(mpdta)

# An `est_method` for did::att_gt(). Signature fixed by `did`; it is the same
# one DRDID::reg_did_panel has. BART estimates the outcome nuisance
# E[dY | X, D = 0]; the IPW weights are 1, so this is the outcome-regression
# estimator, not the doubly robust one.
bart_or <- function(y1, y0, D, covariates, i.weights = NULL, inffunc = FALSE, ...) {
  n  <- length(y1)
  dY <- y1 - y0
  w  <- if (is.null(i.weights)) rep(1, n) else i.weights

  X <- as.data.frame(covariates)
  X <- X[, vapply(X, function(z) length(unique(z)) > 1L, logical(1L)), drop = FALSE]

  if (ncol(X) == 0L) {                       # no covariates: unconditional DiD
    m0 <- rep(stats::weighted.mean(dY[D == 0], w[D == 0]), n)
  } else {
    names(X) <- paste0("v", seq_len(ncol(X)))
    dat <- cbind(dY = dY, X)
    fit <- suppressMessages(bartisan(
      dY ~ ., data = dat[D == 0, , drop = FALSE], weights = w[D == 0],
      family = stats::gaussian(),
      control = bartisan_control(num_burn = 200L, num_draws = 400L,
                                 verbose = FALSE)))
    m0 <- as.vector(stats::predict(fit, newdata = dat))
  }

  pbar <- stats::weighted.mean(D, w)
  att  <- stats::weighted.mean((dY - m0)[D == 1], w[D == 1])

  # Plug-in influence function: it treats the fitted nuisance as known, so it
  # omits the term from having estimated it. Without the propensity part the
  # moment is not Neyman-orthogonal, so that term is first order and this
  # understates. Use it for the point estimate; see the note on inference.
  inf <- (w * D / pbar) * (dY - m0 - att)

  list(ATT = att, att.inf.func = if (inffunc) inf else NULL)
}

fit_one <- function(lab, est) {
  a <- did::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                   gname = "first.treat", xformla = ~ lpop, data = d,
                   control_group = "nevertreated", est_method = est,
                   bstrap = FALSE)
  s <- did::aggte(a, type = "simple", na.rm = TRUE)
  cat(sprintf("%-28s ATT %+.4f  se %.4f  [%+.4f, %+.4f]\n", lab,
              s$overall.att, s$overall.se,
              s$overall.att - 1.96 * s$overall.se,
              s$overall.att + 1.96 * s$overall.se))
  invisible(s)
}
set.seed(4)
fit_one("did, est_method = 'reg'",  "reg")
fit_one("did, est_method = 'dr'",   "dr")
fit_one("did, est_method = bart_or", bart_or)
