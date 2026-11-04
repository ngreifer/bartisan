#' Generalized BART as a SuperLearner library member
#'
#' @description
#' `SL.bartisan()` fits a [bartisan()] model in the form \CRANpkg{SuperLearner}
#' expects of a candidate learner, so that `"SL.bartisan"` can be named in the
#' `SL.library` argument of \pkgfun{SuperLearner}{SuperLearner} alongside any
#' other learner. `predict.SL.bartisan()` predicts from the fit it returns.
#'
#' @param Y the outcome, supplied by \pkg{SuperLearner}.
#' @param X,newX data frames of predictors at which to fit the model and at which
#'   to predict from it, supplied by \pkg{SuperLearner}.
#' @param family the `family` object passed to \pkgfun{SuperLearner}{SuperLearner},
#'   either `gaussian()` or `binomial()`. It is used as the model's family unless
#'   `bartisan_family` is supplied.
#' @param obsWeights observation weights, supplied by \pkg{SuperLearner}. They
#'   are passed to [bartisan()] as prior weights unless every one of them is 1,
#'   in which case none are passed, so that a family that does not take weights,
#'   such as [dpm()], can be requested through `bartisan_family`.
#' @param id cluster identifiers, supplied by \pkg{SuperLearner}; ignored, since
#'   they matter only to how \pkg{SuperLearner} forms its cross-validation
#'   folds.
#' @param bartisan_family optional; the family to fit in place of `family`, in
#'   any of the forms the `family` argument of [bartisan()] accepts. Default is
#'   `NULL` to use `family`. It has to be named differently from `family`
#'   because \pkg{SuperLearner} supplies that argument itself.
#' @param ... further arguments passed to [bartisan()], typically settings of
#'   [bartisan_control()] such as `num_trees` or `num_draws`. See Details for how
#'   to supply them.
#' @param object an `SL.bartisan` object; the `fit` component of the output of
#'   `SL.bartisan()`.
#' @param newdata a data frame of predictors at which to predict.
#'
#' @returns
#' `SL.bartisan()` returns a list with components `pred`, the posterior mean
#' predictions at `newX`, and `fit`, an `SL.bartisan` object holding the
#' `<bartisan_fit>` as its `object` component. `predict.SL.bartisan()` returns a
#' vector of posterior mean predictions, one per row of `newdata`.
#'
#' @details
#' The predictions \pkg{SuperLearner} combines are the posterior means of the
#' response (i.e., [predict.bartisan_fit()] with `type = "response"`), which is
#' the mean of the outcome for `family = gaussian()` and the probability that it
#' is 1 for `family = binomial()`. The ensemble weights are chosen to minimize
#' the cross-validated risk of those means, so the rest of the posterior plays
#' no part in the ensemble; the `<bartisan_fit>` stored in `fit$object` still
#' carries it for anything done with the learner afterwards.
#'
#' \pkg{SuperLearner} calls each learner with a fixed set of arguments, so
#' settings of the model are supplied by making a new learner. Either
#' \pkgfun{SuperLearner}{create.Learner} or a short wrapper function does this
#' (e.g., `SL.bartisan_100 <- function(...) SL.bartisan(..., num_trees = 100)`);
#' see Examples. \pkg{SuperLearner} finds a learner by its name, looking in the
#' environment given by its `env` argument (by default, the one it was called
#' from) and then the search path, so a wrapper must be defined there, and
#' \pkg{bartisan} must be attached with [library()] rather than only loaded for
#' `"SL.bartisan"` itself to be found.
#'
#' Every candidate is fit once per cross-validation fold and once more to the
#' full sample, so the cost of a \pkg{SuperLearner} fit is roughly the number of
#' folds plus one times that of a single [bartisan()] fit. Reducing `num_draws`
#' cuts that cost almost in proportion and changes only how precisely the
#' posterior means are estimated.
#'
#' A continuous outcome is fit with `gaussian()` by default because that is the
#' family \pkg{SuperLearner} names; the [bartisan()] default for an unnamed
#' family on a continuous outcome is [dpm()], which can be requested by setting
#' `bartisan_family = dpm()`. Either family targets the same conditional mean.
#'
#' @seealso
#' \pkgfun{SuperLearner}{SuperLearner} for fitting the ensemble; [bartisan()] for
#' the model and the arguments that can be passed through `...`;
#' [`bartisan-parsnip`] and [mice.impute.bartisan()] for the other prediction
#' workflows \pkg{bartisan} fits into; `vignette("workflows")` for a worked
#' example of all three
#'
#' @examplesIf rlang::is_installed("SuperLearner")
#' library("SuperLearner")
#' data("rhc")
#' set.seed(123)
#'
#' X <- rhc[c("age", "sex", "aps", "meanbp", "pafi", "crea")]
#' Y <- rhc$death
#'
#' # A smaller, faster version of the learner, made by
#' # passing settings through `...`
#' SL.bartisan_small <- function(...) {
#'   SL.bartisan(..., num_trees = 10, num_burn = 50,
#'               num_draws = 50)
#' }
#'
#' sl <- SuperLearner(Y, X, family = binomial(),
#'                    SL.library = c("SL.mean", "SL.glm",
#'                                   "SL.bartisan_small"),
#'                    cvControl = list(V = 3))
#' sl
#'
#' # The bartisan fit itself, for anything beyond the
#' # ensemble's predictions
#' fit <- sl$fitLibrary$SL.bartisan_small_All$object
#' head(predict(fit, type = "response"))
#'
#' @export
SL.bartisan <- function(Y, X, newX, family, obsWeights, id,
                        bartisan_family = NULL, ...) {
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)

  # The response gets a name that cannot collide with a predictor, since `X`
  # may well contain a column called `Y`.
  response <- make.unique(c(names(X), "Y"))[ncol(X) + 1L]

  data <- X
  data[[response]] <- Y

  # SuperLearner passes weights of 1 when the caller gave none, and passing them
  # on would make `bartisan_family = dpm()` an error for no reason.
  weights <- {
    if (is_null(obsWeights) || all(obsWeights == 1)) NULL
    else obsWeights
  }

  formula <- stats::reformulate(".", response = response)

  fit <- bartisan(formula,
                  data = data,
                  family = bartisan_family %or% family,
                  weights = weights,
                  ...)

  list(pred = stats::predict(fit, newdata = newX, type = "response"),
       fit = structure(list(object = fit), class = "SL.bartisan"))
}

#' @rdname SL.bartisan
#' @export
predict.SL.bartisan <- function(object, newdata, ...) {
  stats::predict(object[["object"]], newdata = as.data.frame(newdata),
                 type = "response")
}

