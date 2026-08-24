#' Counterfactual estimands with marginaleffects
#'
#' @description
#' A `genbart` fit works with the \pkg{marginaleffects} package, so that
#' predictions, comparisons and slopes -- and hypothesis tests on any of them --
#' can be computed without extracting draws by hand. Load \pkg{marginaleffects}
#' and call its functions on the fit directly; there is nothing to set up.
#'
#' @details
#' # How the uncertainty is computed
#'
#' A forest has no coefficient vector and no variance-covariance matrix, so the
#' delta method \pkg{marginaleffects} uses for a frequentist model has nothing
#' to work with. It has something better here: the posterior draws. Every
#' estimand is computed by pushing all of the draws through the same
#' transformation and summarizing at the end, so an interval is a posterior
#' quantile rather than a normal approximation, and a nonlinear estimand needs no
#' approximation at all. This is the same path \pkg{marginaleffects} takes for
#' \pkg{brms} and \pkg{rstanarm} fits.
#'
#' # Slopes need a linear predictor transform
#'
#' A slope is a numerical derivative, and taking one requires the fitted function
#' to be differentiable in the predictor *as the caller supplies it*. The default
#' `x_transform = "quantile"` maps each predictor through its empirical
#' distribution function before any rule sees it, and an empirical distribution
#' function is a step function -- so the fit is a step function of the original
#' predictor whatever the decision rules are, and its difference quotient grows
#' without bound as the step shrinks. Measured on a smooth surface where the
#' average derivative is zero:
#'
#' | step | `x_transform = "quantile"` | `x_transform = "range"` |
#' |---|---|---|
#' | 1e-4 | -4.79 | -0.28 |
#' | 1e-2 | -0.40 | -0.30 |
#' | 5e-2 | -0.29 | -0.25 |
#'
#' So **use `x_transform = "range"` if slopes are the estimand**, which maps each
#' predictor linearly and leaves a soft-rule fit differentiable. Hard rules give a
#' piecewise-constant fit under either transform, and a derivative of one is not
#' a meaningful quantity however it is computed.
#'
#' None of this affects [marginaleffects::predictions()] or
#' [marginaleffects::comparisons()], which evaluate the fit at two points a
#' substantive distance apart rather than dividing by a vanishing one. Those are
#' the estimands to reach for with the default transform.
#'
#' # What is not covered
#'
#' `type = "link"` is refused for a family with more than one additive predictor
#' -- `location_scale()`, the zero-inflated families, `multinomial()` -- because
#' there is no single link there to be talking about. Those families work on the
#' response scale, which is one number per observation whatever the family, and
#' `multinomial()` and `ordinal()` work on the probability scale, which gives one
#' group per category. To reach a *particular* predictor of a multi-predictor
#' family, call [predict.genbart()] directly.
#'
#' Extrapolation is worth keeping in mind for `comparisons()`: a forest is
#' constant outside the range of the predictor it was fitted on, so a contrast
#' that steps a predictor beyond that range reports the boundary value rather
#' than an extrapolated one.
#'
#' @param model,x,object,formula a fitted `genbart` object. The last is named
#'   that only because [stats::model.frame()] names its first argument that way.
#' @param newdata data at which to evaluate the fit. Defaults to the data the
#'   model was fitted to, which is retained in the fit for this purpose.
#' @param type `"response"` for the fitted mean, `"link"` for the additive
#'   predictor, or `"prob"` for the per-category probabilities of a categorical
#'   family. A categorical family has no mean, so `"response"` gives `"prob"`
#'   there.
#' @param coefs ignored; a forest has no coefficient vector.
#' @param ... further arguments, passed on to [predict.genbart()] where it makes
#'   sense to and ignored otherwise.
#'
#' @returns
#' `get_predict()` returns a data frame with columns `rowid`, `group` and
#' `estimate`, carrying the draws in a `"posterior_draws"` attribute of
#' observations by draws, which is the interface \pkg{marginaleffects}
#' documents. The other methods exist to satisfy the generic and return what
#' their names suggest.
#'
#' @examplesIf requireNamespace("marginaleffects", quietly = TRUE)
#' set.seed(1)
#' n <- 200
#' d <- data.frame(x1 = runif(n), x2 = runif(n))
#' d$y <- 2 * sin(pi * d$x1) - d$x2 + rnorm(n)
#'
#' fit <- genbart(y ~ x1 + x2, data = d,
#'                control = genbart_control(num_trees = 10, num_burn = 50,
#'                                          num_save = 50))
#'
#' marginaleffects::avg_predictions(fit)
#' marginaleffects::avg_comparisons(fit)
#'
#' @name genbart-marginaleffects
#' @importFrom stats family formula terms model.frame nobs
NULL

# Support for marginaleffects, and for the insight package it reaches through.
#
# The point of this file is that a fitted forest has no coefficient vector and no
# variance-covariance matrix, so the delta method marginaleffects uses for
# frequentist models has nothing to work with. What it does have is the posterior
# draws, and marginaleffects handles that case natively: a `get_predict()` method
# that attaches the draws lets every one of its estimands -- predictions,
# comparisons, slopes, and any hypothesis test on them -- be computed by pushing
# the draws through the same transformation and summarizing at the end. That is
# the correct thing for this model rather than a workaround: the interval comes
# from the posterior rather than from a normal approximation to it.
#
# The methods are registered at load time only if the packages are present, so
# they stay in Suggests.

# Base generics first. insight falls back on these for any class it does not know
# specifically, so defining them is most of what makes the fit legible to it.

#' @rdname genbart-marginaleffects
#' @export
formula.genbart <- function(x, ...) {
  x[["formula"]]
}

#' @rdname genbart-marginaleffects
#' @export
terms.genbart <- function(x, ...) {
  x[["terms"]]
}

#' @rdname genbart-marginaleffects
#' @export
model.frame.genbart <- function(formula, ...) {
  formula[["model"]]
}

#' @rdname genbart-marginaleffects
#' @export
nobs.genbart <- function(object, ...) {
  object[["n"]]
}

#' @rdname genbart-marginaleffects
#' @export
family.genbart <- function(object, ...) {
  object[["family"]]
}

# The number of retained draws, which is also the number of rows every posterior
# matrix in the fit has.
num_draws <- function(object) {
  nrow(object[["sigma_mu"]])
}

# marginaleffects asks for a type by name. Its vocabulary is not quite ours, so
# the mapping is explicit rather than passed through: "link" is the additive
# predictor, "response" the fitted mean, and "prob" the per-category
# probabilities that a categorical family reports instead of a mean.
me_type <- function(object, type) {
  categorical <- object[["family"]][["family"]] %in%
    c("binomial", "ordinal", "multinomial")

  if (is_null(type) || identical(type, "response")) {
    return(if (categorical) "prob" else "response")
  }

  if (identical(type, "link")) {
    return("link")
  }

  if (identical(type, "prob")) {
    return("prob")
  }

  arg::err("{.arg type} must be {.val response}, {.val link} or {.val prob}")
}

# Draws of whatever quantity was asked for, as a matrix of draws by rows -- or,
# for a categorical family on the probability scale, a three-dimensional array of
# draws by rows by categories.
me_draws <- function(model, newdata, type) {
  out <- stats::predict(model, newdata = newdata, type = type, draws = TRUE)

  # A family with several additive predictors returns one matrix per predictor on
  # the link scale, and there is no single quantity for marginaleffects to be
  # talking about. The response and probability scales are a different matter:
  # those are one number per observation whatever the family, so they go through.
  if (is.list(out)) {
    if (length(out) > 1L) {
      arg::err("this family has {length(out)} additive predictors, so
                {.val link} is not a single quantity. Ask for {.val response},
                or use {.fn predict} directly to reach a particular predictor")
    }

    out <- out[[1L]]
  }

  out
}

#' @rdname genbart-marginaleffects
#' @exportS3Method marginaleffects::get_predict
get_predict.genbart <- function(model, newdata = NULL, type = NULL, ...) {
  if (is_null(newdata)) {
    newdata <- model[["model"]]
  }

  resolved <- me_type(model, type)
  draws <- me_draws(model, newdata, resolved)

  if (length(dim(draws)) == 3L) {
    # Draws by rows by categories. marginaleffects wants the rows stacked
    # category by category, and the draws matrix stacked the same way, so that
    # its row order and the draws' row order agree.
    levs <- dimnames(draws)[[3L]]

    if (is_null(levs)) {
      levs <- as.character(seq_len(dim(draws)[3L]))
    }

    estimate <- apply(draws, c(2L, 3L), mean)

    out <- data.frame(
      rowid = rep(seq_len(nrow(estimate)), times = ncol(estimate)),
      group = rep(levs, each = nrow(estimate)),
      estimate = as.vector(estimate))

    flat <- do.call(cbind, lapply(seq_len(dim(draws)[3L]),
                                  function(k) draws[, , k]))
    attr(out, "posterior_draws") <- t(flat)

    return(out)
  }

  out <- data.frame(rowid = seq_len(ncol(draws)),
                    group = "main_marginaleffect",
                    estimate = colMeans(draws))
  rownames(out) <- NULL
  attr(out, "posterior_draws") <- t(draws)

  out
}

#' @rdname genbart-marginaleffects
#' @exportS3Method marginaleffects::get_group_names
get_group_names.genbart <- function(model, ...) {
  levs <- model[["levels"]]

  if (is_null(levs)) {
    return("main_marginaleffect")
  }

  levs
}

# A forest has no coefficient vector. Returning an empty one rather than erroring
# is what lets marginaleffects reach the draws path: it checks for coefficients,
# finds none to perturb, and falls back on the posterior.
#' @rdname genbart-marginaleffects
#' @exportS3Method marginaleffects::get_coef
get_coef.genbart <- function(model, ...) {
  stats::setNames(numeric(0), character(0))
}

#' @rdname genbart-marginaleffects
#' @exportS3Method marginaleffects::set_coef
set_coef.genbart <- function(model, coefs, ...) {
  model
}

# No variance-covariance matrix, and this is not a gap to be filled: the
# uncertainty is in the draws, and a normal approximation built from a
# pseudo-vcov would be a worse answer than the posterior marginaleffects already
# has.
#' @rdname genbart-marginaleffects
#' @exportS3Method marginaleffects::get_vcov
get_vcov.genbart <- function(model, ...) {
  NULL
}

#' @rdname genbart-marginaleffects
#' @exportS3Method insight::get_data
get_data.genbart <- function(x, ...) {
  x[["model"]]
}

.onLoad <- function(libname, pkgname) {
  # marginaleffects checks a fitted object's class against a list of the ones it
  # knows how to handle, before any of the methods above are reached, and exposes
  # this option as the way a package outside it registers its own class. Setting
  # an option on load is not something to do lightly, but this is the documented
  # mechanism and the alternative is that the methods are never called.
  classes <- getOption("marginaleffects_model_classes", default = NULL)

  if (!"genbart" %in% classes) {
    options(marginaleffects_model_classes = c(classes, "genbart"))
  }

  invisible()
}
