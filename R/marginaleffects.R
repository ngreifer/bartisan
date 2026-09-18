#' Counterfactual estimands with marginaleffects
#'
#' @description
#' Methods that let a `<bartisan_fit>` object be used by the
#' \pkg{marginaleffects} package, which computes predictions, comparisons and
#' slopes (and hypothesis tests on any of them) from the posterior draws. There
#' is nothing to set up: load \pkg{marginaleffects} and call its functions on
#' the fit.
#'
#' @param model,x,object,formula a `<bartisan_fit>` object; the output of a call
#'   to [bartisan()]. The last is named `formula` only because
#'   [stats::model.frame()] names its first argument that way.
#' @param newdata optional; a data frame at which to evaluate the fit. Default
#'   is `NULL` to use the data the model was fitted to, which is retained in the
#'   fit for this purpose.
#' @param type string; the scale to work on. Allowable options include
#'   `"response"` (the default), `"link"`, `"prob"`, `"mean"`, `"stdlv"`, and
#'   `"survival"`. `"response"` is the fitted mean, `"link"` the additive
#'   predictor, `"prob"` the per-category probabilities of a categorical family,
#'   `"mean"` the mean of a categorical response with its labels read as
#'   numbers, `"stdlv"` the standardized latent variable of an ordinal fit, and
#'   `"survival"` the survival function at the times given in `times`; all but
#'   the first two are described under [predict.bartisan_fit()]. An `ordinal()`
#'   or `multinomial()` response has no single mean, so `"response"` gives
#'   `"prob"` for those families, while a binomial fit reports the probability
#'   of its second level, as [stats::glm()] does, and `"prob"` is what asks for
#'   both of its columns. `"probs"`, `"lp"`, `"lv"` and `"surv"` are accepted as
#'   aliases for `"prob"`, `"link"`, `"link"` and `"survival"`, since those are
#'   the names the same quantities go by for other ordinal fits in
#'   \pkg{marginaleffects}.
#' @param coefs ignored; a forest has no coefficient vector.
#' @param ... further arguments. `values`, `iterations`, `offset`, `weights`,
#'   `log` and `times` are passed on to [predict.bartisan_fit()]; anything else
#'   is ignored, since \pkg{marginaleffects} puts arguments of its own here too.
#'   Note that \pkg{marginaleffects} warns that it does not recognize `values`,
#'   which is expected (it is this package's argument, not one of its own), and
#'   the value is used regardless.
#'
#' @returns
#' `get_predict()` returns a data frame with columns `rowid`, `group` and
#' `estimate`, carrying the draws in a `"posterior_draws"` attribute of
#' observations by draws, which is the interface \pkg{marginaleffects}
#' documents. The other methods exist to satisfy the generic and return what
#' their names suggest.
#'
#' @details
#' ## Uncertainty
#'
#' A forest has no coefficient vector and no variance-covariance matrix, so the
#' delta method \pkg{marginaleffects} uses for a frequentist model has nothing to
#' work with. It has the posterior draws instead: every estimand is computed by
#' pushing all of the draws through the same transformation and summarizing at
#' the end, so an interval is a posterior quantile rather than a normal
#' approximation, and a nonlinear estimand needs no approximation at all. This is
#' the same path \pkg{marginaleffects} takes for \pkg{brms} and \pkg{rstanarm}
#' fits.
#'
#' One consequence is worth knowing: \pkg{marginaleffects} centers a posterior at
#' its median, where [predict.bartisan_fit()] reports its mean. The two summarize
#' the same draws, so a difference between them is the skewness of the posterior
#' rather than a disagreement.
#'
#' ## A Contrast of Exactly Zero Is Usually Real
#'
#' `avg_comparisons()` reporting an estimate of exactly `0` is the most common
#' surprise here, and two facts meet to produce it.
#'
#' The posterior of a contrast has an atom at exactly zero. In any draw where no
#' tree splits on the variable being contrasted, the fit does not depend on that
#' variable at all, so the two counterfactual predictions are identical to the
#' last bit and their difference is exactly zero. The Dirichlet sparsity prior,
#' `sparsity` in [bartisan_control()], is what makes those draws common, since
#' dropping a weak predictor from every tree is what a variable selection prior
#' is for. And \pkg{marginaleffects} centers a posterior at its median, so once
#' the atom holds more than half the mass the reported estimate is exactly zero
#' however large the rest of the posterior is: the median can land on the atom
#' while the posterior mean and the upper limit of the interval are both far from
#' zero.
#'
#' Four things are worth doing about it, in the order given. Look at
#' `prop_used` in `summary()`, the posterior probability that the predictor
#' appears anywhere in the forest, which is what the zero reports. Ask for the
#' mean instead, with `options(marginaleffects_posterior_center = mean)`, which
#' is the summary [predict.bartisan_fit()] reports and the one that behaves
#' sensibly against an atom. Reconsider the sparsity prior, since
#' `sparsity = FALSE` removes the atom almost entirely where a larger `num_trees`
#' leaves it in place; that is a modeling choice rather than a fix, and it is the
#' right one when a contrast on a particular predictor is the estimand. And run
#' several chains and compare them, because the variable selection state mixes
#' slowly and a single chain can look much more settled than the posterior is.
#' `vignette("effects")` works the choice through.
#'
#' ## Slopes and the Predictor Transform
#'
#' A slope is a numerical derivative, and taking one requires the fitted function
#' to be differentiable in the predictor as the caller supplies it. Whether it is
#' depends on `x_transform` in [bartisan_control()], because the fit is a smooth
#' function of the transformed predictor rather than of the original one.
#'
#' `x_transform = "range"` is what slopes want, and refitting with it is
#' worthwhile when a derivative is the quantity being reported rather than a
#' prediction. It is the only affine map of the three, so it is the only one
#' through which a slope of the fit is a slope of the original predictor. The
#' default, `"smoothcdf"`, carries an estimated density along with it and lands
#' close; `"quantile"` maps each predictor through a step function and so has no
#' derivative to take, and the number returned there is a property of the step
#' size rather than of the fit. `vignette("effects")` shows what each does as the
#' step shrinks.
#'
#' None of this affects \pkgfun{marginaleffects}{predictions} or
#' \pkgfun{marginaleffects}{comparisons}, which evaluate the fit at two points a
#' substantive distance apart rather than dividing by a vanishing one.
#'
#' ## The Usual Survival Estimand
#'
#' For a survival family, the estimand is usually a contrast in survival at a
#' horizon rather than in the predictor. `type = "survival"` with `times` gives
#' it:
#'
#' ```r
#' # The difference in one-year survival between treated and untreated
#' avg_comparisons(fit, variables = "trt", type = "survival", times = 1)
#' ```
#'
#' One time per call. \pkg{marginaleffects} checks the dots against a whitelist
#' of its own, hardcoded per model class, so it warns that it does not recognize
#' `times` while passing it through, which is what the warning says. There is no
#' hook for registering an argument with it, so the warning is expected and the
#' result is correct.
#'
#' ## Prediction Types
#'
#' `type = "link"` needs a family with a single additive predictor, since there
#' is no one link to be talking about in `gaussian_ls()`, the zero-inflated
#' families or `multinomial()`. Those families work on the response scale, which
#' is one number per observation whatever the family, and `multinomial()` and
#' `ordinal()` work on the probability scale, which gives one group per category.
#' [predict.bartisan_fit()] reaches a particular predictor of a multi-predictor
#' family, and computes `type = "class"` and `type = "density"`, which are not
#' one number per observation that an average could be taken of.
#'
#' Extrapolation is worth keeping in mind for `comparisons()`: a forest is
#' constant outside the range of the predictor it was fitted on, so a contrast
#' that steps a predictor beyond that range reports the boundary value rather
#' than an extrapolated one.
#'
#' @seealso
#' [predict.bartisan_fit()] for the prediction scales these estimands are
#' computed on; [bartisan_control()] for `sparsity` and `x_transform`;
#' [bartisan-interop] for the methods that let other packages assess the fit
#'
#' @examplesIf rlang::is_installed("marginaleffects")
#' data("rhc")
#' set.seed(123)
#'
#' # Whether a patient died, with every other variable a candidate predictor.
#' # The sparsity prior is turned off because a contrast on one predictor is
#' # the estimand rather than variable selection
#' fit <- bartisan(death ~ . - days, data = rhc, sparsity = FALSE,
#'                 num_trees = 10, num_burn = 50, num_draws = 50)
#'
#' # The effect of catheterization on the probability of death, as an average
#' # contrast of counterfactual predictions
#' marginaleffects::avg_comparisons(fit, variables = "rhc")
#'
#' # The same contrast within each sex, and a test that the two are equal
#' marginaleffects::avg_comparisons(fit, variables = "rhc", by = "sex",
#'                                  hypothesis = ~pairwise)
#'
#' # Centering each posterior at its mean rather than its median, which is the
#' # summary `predict()` reports
#' op <- options(marginaleffects_posterior_center = mean)
#' marginaleffects::avg_comparisons(fit, variables = "rhc")
#' options(op)
#'
#' @name bartisan-marginaleffects
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

#' @rdname bartisan-marginaleffects
#' @export
formula.bartisan_fit <- function(x, ...) {
  x[["formula"]]
}

#' @rdname bartisan-marginaleffects
#' @export
terms.bartisan_fit <- function(x, ...) {
  x[["terms"]]
}

#' @rdname bartisan-marginaleffects
#' @export
model.frame.bartisan_fit <- function(formula, ...) {
  formula[["model"]]
}

#' @rdname bartisan-marginaleffects
#' @export
nobs.bartisan_fit <- function(object, ...) {
  object[["n"]]
}

#' @rdname bartisan-marginaleffects
#' @export
family.bartisan_fit <- function(object, ...) {
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
# probabilities that a categorical family reports instead of a mean. The aliases
# are the names the same quantities go by elsewhere in marginaleffects -- "probs"
# and "lp" are what it uses for the WeightIt classes -- so that a call written
# for one ordinal fit does not have to be rewritten for this one.
#
# The types that do not appear here are the ones that are not a number per
# observation and so have nothing marginaleffects could average or contrast:
# "class" is a factor, and "density" needs the outcome, which a counterfactual
# grid does not have.
me_type <- function(object, type) {
  # A binomial fit is deliberately not in this list. Its two probabilities sum to
  # one, so reporting both gives every estimand twice, as mirror images, where
  # `glm()` gives one row; and because a binary outcome carries no `levels`, the
  # two rows were labelled by the outcome's own values rather than by anything
  # meaningful. Reporting the probability of the second level, which is what
  # `predict(type = "response")` returns, matches `glm()` and keeps
  # `hypothesis = ~pairwise` usable across subgroups. `type = "prob"` still asks
  # for both columns.
  categorical <- object[["family"]][["family"]] %in%
    c("ordinal", "multinomial", "mnp")

  if (is_null(type) || identical(type, "response")) {
    return(if (categorical) "prob" else "response")
  }

  known <- c(link = "link", lp = "link", lv = "link",
             prob = "prob", probs = "prob",
             mean = "mean", stdlv = "stdlv",
             # The survival function at named times, which is what makes the
             # usual survival estimand -- a contrast in t-year survival --
             # reachable through the estimand functions. The times come through
             # the dots and become the groups.
             survival = "survival", surv = "survival")

  if (!type %in% names(known)) {
    arg::err("{.arg type} must be {.or {.val {c('response', unique(names(known)))}}}")
  }

  unname(known[[type]])
}

# Draws of whatever quantity was asked for, as a matrix of draws by rows -- or,
# for a categorical family on the probability scale, a three-dimensional array of
# draws by rows by categories.
me_draws <- function(model, newdata, type, extra = list()) {
  args <- c(list(model, newdata = newdata, type = type, draws = TRUE), extra)
  out <- do.call(stats::predict, args)

  # One time is one quantity, not a group of one: dropping the trailing dimension
  # sends it down the ungrouped path, which is what makes a contrast in t-year
  # survival read like any other contrast.
  if (length(dim(out)) == 3L && dim(out)[3L] == 1L) {
    out <- out[, , 1L, drop = TRUE]
  }

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

#' @rdname bartisan-marginaleffects
#' @exportS3Method marginaleffects::get_predict
get_predict.bartisan_fit <- function(model, newdata = NULL, type = NULL, ...) {
  if (is_null(newdata)) {
    newdata <- model[["model"]]
  }

  resolved <- me_type(model, type)
  draws <- me_draws(model, newdata, resolved, predict_args(...))

  # marginaleffects prepends rows of its own, marked `rowid = -1`, and drops
  # them again by that marker once the predictions come back. Regenerating the
  # column as `seq_len(n)` destroys the marker, and those rows then leak into the
  # result: `predictions(fit, newdata = d[1:3, ])` returned five rows, two of
  # them marginaleffects' own. Carry the column through when it is there.
  row_id <- newdata[["rowid"]] %or% seq_len(nrow(newdata))

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
      rowid = rep.int(row_id, ncol(estimate)),
      group = rep(levs, each = nrow(estimate)),
      estimate = as.vector(estimate))

    flat <- lapply(seq_len(dim(draws)[3L]),
                   function(k) draws[, , k]) |>
      do_cbind()

    attr(out, "posterior_draws") <- t(flat)

    return(out)
  }

  out <- data.frame(rowid = row_id,
                    group = "main_marginaleffect",
                    estimate = colMeans(draws)) |>
    unrowname()

  attr(out, "posterior_draws") <- t(draws)

  out
}

# marginaleffects calls get_predict() with arguments of its own in the dots --
# `mfx`, and whatever the user passed to the estimand function -- so the dots
# cannot simply be forwarded to predict(); a stray name there would match one of
# predict()'s arguments positionally or partially. Only the arguments predict()
# actually has a use for are taken, by exact name.
predict_args <- function(...) {
  dots <- list(...)

  keep <- intersect(names(dots),
                    c("values", "iterations", "offset", "weights", "log",
                      "times"))

  dots[keep]
}

#' @rdname bartisan-marginaleffects
#' @exportS3Method marginaleffects::get_group_names
get_group_names.bartisan_fit <- function(model, ...) {
  model[["levels"]] %or% "main_marginaleffect"
}

# A forest has no coefficient vector. Returning an empty one rather than erroring
# is what lets marginaleffects reach the draws path: it checks for coefficients,
# finds none to perturb, and falls back on the posterior.
#' @rdname bartisan-marginaleffects
#' @exportS3Method marginaleffects::get_coef
get_coef.bartisan_fit <- function(model, ...) {
  setNames(numeric(), character())
}

#' @rdname bartisan-marginaleffects
#' @exportS3Method marginaleffects::set_coef
set_coef.bartisan_fit <- function(model, coefs, ...) {
  model
}

# No variance-covariance matrix, and this is not a gap to be filled: the
# uncertainty is in the draws, and a normal approximation built from a
# pseudo-vcov would be a worse answer than the posterior marginaleffects already
# has.
#' @rdname bartisan-marginaleffects
#' @exportS3Method marginaleffects::get_vcov
get_vcov.bartisan_fit <- function(model, ...) {
  NULL
}

#' @rdname bartisan-marginaleffects
#' @exportS3Method insight::get_data
get_data.bartisan_fit <- function(x, ...) {
  flatten_matrix_columns(x[["model"]])
}

# A survival response arrives as a two-column matrix -- `Surv(time, status)`, or
# `cbind(time, status)` -- and a model frame keeps it as a single matrix column.
# marginaleffects converts what it is given to a data.table, which cannot hold
# one: it reads the 2n cells as a column of length 2n and fails on the mismatch
# with every other column. That took out `avg_comparisons()` and its siblings for
# every survival family, whatever `type` was asked for, with an error that named
# data.table rather than the cause.
#
# The estimands need the predictors, not the response, so the matrix is split
# into ordinary columns. Their names are not meant to be used; they exist so the
# frame is rectangular.
flatten_matrix_columns <- function(data) {
  wide <- vapply(data, function(column) {
    is.matrix(column) && ncol(column) > 1L
  }, logical(1))

  if (!any(wide)) {
    return(data)
  }

  pieces <- lapply(names(data), function(nm) {
    column <- data[[nm]]

    if (!(is.matrix(column) && ncol(column) > 1L)) {
      return(setNames(list(column), nm))
    }

    lapply(seq_len(ncol(column)), function(j) column[, j]) |>
      setNames(paste0(nm, ".", seq_len(ncol(column))))
  })

  out <- do.call(c, pieces)
  structure(as.data.frame(out, stringsAsFactors = FALSE,
                          check.names = FALSE),
            terms = attr(data, "terms"))
}

.onLoad <- function(libname, pkgname) {
  # marginaleffects checks a fitted object's class against a list of the ones it
  # knows how to handle, before any of the methods above are reached, and exposes
  # this option as the way a package outside it registers its own class. Setting
  # an option on load is not something to do lightly, but this is the documented
  # mechanism and the alternative is that the methods are never called.
  classes <- getOption("marginaleffects_model_classes", default = NULL)

  if (!"bartisan_fit" %in% classes) {
    options(marginaleffects_model_classes = c(classes, "bartisan_fit"))
  }

  invisible()
}
