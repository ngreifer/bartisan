#' Predictions from a generalized BART model
#'
#' @description
#' Evaluates the stored posterior draws of the forests, either at the data used
#' to fit the model or at new data. Because every draw of every tree is kept,
#' predictions carry full posterior uncertainty rather than being a single point
#' estimate.
#'
#' @param object a fitted model from [bartisan()].
#' @param newdata optional data frame at which to predict. Omit it to use the
#'   data the model was fit to. Missing predictor values are allowed in the
#'   columns that had them when the model was fit, since only those columns'
#'   splitting rules carry an answer for one; a missing value anywhere else is an
#'   error. See the Missing predictor values section of [bartisan()].
#' @param type the scale of the prediction:
#'   \describe{
#'     \item{`"link"`}{the additive predictor, one column per predictor for
#'       families that have more than one.}
#'     \item{`"response"`}{the mean of the response; the median survival time
#'       for the accelerated failure time families; and, for a response with
#'       categories, the category probabilities, since there is no single mean
#'       to report.}
#'     \item{`"prob"`}{category probabilities, for the binomial, ordinal and
#'       multinomial families.}
#'     \item{`"class"`}{the most probable category, as a factor, for the same
#'       families.}
#'     \item{`"mean"`}{the mean of the response with the category labels read as
#'       numbers, for the same families. `"4"` counts as four. This is the
#'       summary an ordinal outcome with numeric labels usually wants, and it
#'       needs no assumption at the modeling stage: the model is still ordinal
#'       and only the reporting treats the categories as numbers. Use `values` to
#'       say what the categories are worth when the labels are not numbers, or
#'       are not the numbers you mean.}
#'     \item{`"stdlv"`}{the additive predictor divided by the standard deviation
#'       of the latent variable it indexes, for the ordinal and binomial
#'       families. Either response can be written as a threshold crossing of a
#'       continuous `y* = eta + e`, and the link fixes the distribution of `e`
#'       and so its variance; dividing by the standard deviation of `y*` puts
#'       fits with different links, or different amounts of signal, on one
#'       scale, which is what a standardized effect size on such an outcome
#'       needs. Available for the probit, logit and complementary log-log links,
#'       which are the ones with a latent distribution to name. See Details.}
#'     \item{`"density"`}{the conditional density of the outcome given the
#'       predictors, evaluated at the observed outcome. This requires the
#'       outcome, so `newdata` must contain it; omit `newdata` to use the data
#'       the model was fit to. The value is the likelihood contribution of the
#'       observation, so it is a density for a continuous response, a
#'       probability for a discrete one, and a survival probability for a
#'       censored survival time. Useful for held-out log scores and for
#'       posterior predictive checks. **The measure differs across the survival
#'       families**: the accelerated failure time families, `dpm_aft()` included,
#'       report the density of \eqn{\log T}, while `ph()` reports the density of
#'       \eqn{T}. The two differ by \eqn{\sum \log t}, so log
#'       scores are comparable within each group and not across them;
#'       `type = "survival"` is comparable throughout.}
#'     \item{`"survival"`}{the survival function \eqn{S(t \mid x)} at the times
#'       given in `times`, for the accelerated failure time families and `ph()`.
#'       This is the predictive distribution of a survival response, and the
#'       analogue of `"prob"` for a categorical one: `"response"` reports only the
#'       median. Returns one column per time, or a draws by rows by times array
#'       when `draws = TRUE`. It is also what makes the usual survival estimand
#'       reachable through \pkg{marginaleffects} -- a contrast in \eqn{t}-year
#'       survival; see [bartisan-marginaleffects].}
#'   }
#' @param draws return every posterior draw rather than the posterior mean. The
#'   result then gains a leading dimension indexing draws.
#' @param iterations optional integer vector selecting which stored draws to
#'   use. Defaults to all of them.
#' @param offset an offset for `newdata`, on the link scale. Required when the
#'   model was fit with an observation-level offset, because the offset is not a
#'   function of the predictors and so cannot be reconstructed.
#' @param values for `type = "mean"`, a numeric vector named for every response
#'   level, giving what each category is worth. Defaults to the level labels read
#'   as numbers, which fails with an error rather than a guess when they are not
#'   numbers.
#' @param weights prior weights for `newdata`, used only by
#'   `type = "density"`. For a binomial response given as proportions these are
#'   the numbers of trials. Ignored otherwise.
#' @param log for `type = "density"`, return the log of the value instead.
#' @param times for `type = "survival"`, the times at which to report the
#'   survival function. Required, because the horizon is a choice rather than a
#'   property of the fit.
#'   Summing across observations then gives a log score. With `draws = FALSE`
#'   the density is averaged over draws before the log is taken, so the result
#'   is the pointwise predictive density rather than the average log density.
#' @param ... ignored, present for compatibility with the generic.
#'
#' @returns
#' With `draws = FALSE`, a vector for a single-predictor family on the `"link"`
#' or `"response"` scale, a matrix of observations by categories for `"prob"`, a
#' factor for `"class"`, and a matrix of observations by predictors otherwise.
#' With `draws = TRUE`, a matrix of draws by observations, or a list of such
#' matrices when the family has several additive predictors, or an array of
#' draws by observations by categories for `"prob"`.
#'
#' @details
#' # Multinomial probit probabilities are simulated
#'
#' The likelihood of a `multinomial("probit")` fit has no closed form: the
#' probability of a category is the chance that the largest of several correlated
#' Gaussian variables is the one belonging to it, which is a multivariate orthant
#' probability. So `"prob"`, `"class"`, `"response"` and `"density"` are all
#' simulated, using the number of replicates the family was given. Fresh draws
#' are taken on every call, so two calls differ by Monte Carlo error; that error
#' is per posterior draw and averages down over them, which makes `draws = FALSE`
#' much more accurate than any single row of `draws = TRUE`.
#'
#' # The standardized latent variable
#'
#' `type = "stdlv"` reports `(eta - E[e]) / sd(y*)` for the latent
#' `y* = eta + e`, following
#' \pkgfun{WeightIt}{predict.ordinal_weightit}. Three parts of that need saying.
#'
#' The **scale** is `sd(y*) = sqrt(var(eta) + var(e))`, where `var(eta)` is taken
#' over the sample the model was fitted to, per draw, so it is a property of the
#' model rather than of whatever is being predicted; the same divisor is used when
#' predicting new data. `var(e)` is whatever the link implies: 1 for the probit
#' link, `pi^2 / 3` for the logit, `pi^2 / 6` for the complementary log-log.
#'
#' The **location** subtracts the latent error's mean, which shifts `y*` so that
#' its error is centered. That is invisible for the logit and probit links, whose
#' errors are already centered, and is the whole of the difference for the
#' complementary log-log link, whose error is a smallest extreme value variate.
#'
#' The **sign of that shift differs between the two families**, and only for the
#' complementary log-log link. A normal or logistic error is symmetric, so it does
#' not matter whether `e` or `-e` is the thing added to the index. A smallest
#' extreme value error is not symmetric, and the two families add it with opposite
#' signs: an ordinal model has `P(Y <= k) = G(c_k - eta)`, which is
#' `P(eta + e <= c_k)`, so its error has mean `-gamma`; a binomial model has
#' `P(Y = 1) = G(eta)`, which is `P(e <= eta)`, so its latent is `eta - e` and the
#' error has mean `+gamma`. The two are different models rather than the same one
#' written twice, which is also why a two-category ordinal complementary log-log
#' fit is not the same as a binomial one.
#'
#' Such a model is identified only up to a common shift of its thresholds and its
#' predictor, so the location of this quantity is a convention rather than a fact,
#' and the one used here is the same one the cutpoints use -- a predictor centered
#' over the fitted sample. Against \pkg{WeightIt}, which identifies by dropping the
#' intercept column instead, the two agree on the scale and differ by a constant;
#' measured on a linear truth, the standard deviations agree to under 1% and the
#' difference is constant to three decimals. Differences on this scale, which is
#' what a standardized quantity is for, are unaffected.
#'
#' @seealso [bartisan()]
#'
#' @examples
#' set.seed(1)
#'
#' n <- 150
#' d <- data.frame(x1 = runif(n), x2 = runif(n))
#' d$y <- rpois(n, exp(0.5 + d$x1))
#'
#' fit <- bartisan(y ~ x1 + x2, data = d, family = poisson(),
#'                control = bartisan_control(num_trees = 10, num_burn = 50,
#'                                          num_save = 50, verbose = FALSE))
#'
#' head(predict(fit))
#'
#' nd <- data.frame(x1 = c(0.1, 0.9), x2 = c(0.5, 0.5))
#' predict(fit, newdata = nd)
#'
#' # Posterior uncertainty for the two new points
#' apply(predict(fit, newdata = nd, draws = TRUE), 2,
#'       quantile, c(0.025, 0.975))
#'
#' # Held-out log score: the conditional density needs the outcome, so it must
#' # be present in `newdata`.
#' held_out <- data.frame(x1 = runif(20), x2 = runif(20))
#' held_out$y <- rpois(20, exp(0.5 + held_out$x1))
#' sum(log(predict(fit, newdata = held_out, type = "density")))
#'
#' @export
predict.bartisan_fit <- function(object, newdata = NULL, type = "response",
                            draws = FALSE, iterations = NULL, offset = NULL,
                            weights = NULL, values = NULL, log = FALSE,
                            times = NULL, ...) {

  type <- arg::match_arg(type, c("link", "response", "prob", "class",
                                 "mean", "stdlv", "density", "survival"))
  arg::arg_flag(draws)
  arg::arg_flag(log)

  family <- object[["family"]][["family"]]
  categorical <- family %in% c("binomial", "ordinal", "multinomial", "mnp")
  survival <- family %in% c("aft", "ph", "dpm_aft")

  if (identical(type, "survival")) {
    if (!survival) {
      arg::err("{.arg type} {.val survival} is available only for the survival
                families, {.fn weibull_aft}, {.fn loglogistic_aft},
                {.fn lognormal_aft} and {.fn ph}")
    }

    if (is_null(times)) {
      arg::err("{.arg times} says at which times to report survival, and has no
                default because the horizon is a choice rather than a property of
                the fit",
               i = "for example {.code times = c(1, 5)}")
    }

    arg::arg_numeric(times)

    if (any(!is.finite(times)) || any(times <= 0)) {
      arg::err("{.arg times} must be finite and strictly positive")
    }
  }
  else if (!is_null(times)) {
    arg::wrn("{.arg times} is ignored unless {.code type = \"survival\"}")
  }

  if (type %in% c("prob", "class", "mean") && !categorical) {
    arg::err("{.arg type} {.val {type}} is available only for the
              {.val binomial}, {.val ordinal} and {.val multinomial} families")
  }

  if (identical(type, "stdlv") && !family %in% c("ordinal", "binomial")) {
    arg::err("{.arg type} {.val stdlv} is available only for the
              {.val ordinal} and {.val binomial} families, whose response is a
              threshold crossing of a latent variable")
  }

  if (!identical(type, "mean") && !is_null(values)) {
    arg::wrn("{.arg values} is ignored unless {.code type = \"mean\"}")
  }

  parts <- predict_parts(object, newdata, offset, iterations)
  eta <- parts[["eta"]]
  aux <- parts[["aux"]]

  if (identical(type, "link")) {
    return(shape_link(eta, draws, object))
  }

  if (identical(type, "stdlv")) {
    return(shape_link(standardized_latent(object, eta), draws, object))
  }

  if (identical(type, "response")) {
    # A custom family supplies a log density, not a mean, so there is no
    # response scale to map onto and the predictor is what comes back.
    if (identical(family, "custom")) {
      return(shape_link(eta, draws, object))
    }

    return(response_scale(object, eta, aux, draws))
  }

  if (identical(type, "density")) {
    return(conditional_density(object, newdata, eta, aux, weights, draws, log,
                               parts[["iterations"]]))
  }

  if (identical(type, "survival")) {
    return(survival_scale(object, eta, aux, times, draws))
  }

  probs <- category_probs(object, eta, aux)

  if (identical(type, "prob")) {
    if (draws) {
      return(probs)
    }
    out <- apply(probs, c(2L, 3L), mean)
    dimnames(out) <- list(NULL, dimnames(probs)[[3L]])
    return(out)
  }

  if (identical(type, "mean")) {
    return(category_mean(object, probs, values, draws))
  }

  mean_probs <- apply(probs, c(2L, 3L), mean)
  levels <- dimnames(probs)[[3L]]

  # An ordinal response came in as an ordered factor, so send it back as one.
  factor(levels[apply(mean_probs, 1L, which.max)], levels = levels,
         ordered = identical(family, "ordinal"))
}

# The additive predictor on the scale of the standard deviation of the latent
# variable it is the index of.
#
# An ordinal or binary response can be written as a threshold crossing of a
# latent continuous response, `y* = eta + e`; the link fixes the distribution of
# `e` and therefore its variance, which is what makes the latent scale identified
# at all. Dividing by the standard deviation of `y*` puts models with different
# links, or with different amounts of signal, on one scale -- which is what a
# standardized effect size on such an outcome needs. This follows
# `WeightIt::predict.ordinal_weightit(type = "stdlv")`.
#
# Two details matter. The variance of `eta` is taken over the **fitted** sample,
# so it is a property of the model rather than of whatever is being predicted;
# and `e` is shifted to have mean zero, which moves its mean into the index. That
# second point is invisible for the logit and probit links, where the error is
# already centered, and is the whole of the difference for the complementary
# log-log link.
standardized_latent <- function(object, eta) {
  error <- latent_error(object[["family"]][["family"]],
                        object[["family"]][["link"]])

  # The predictor's own variance, per draw, over the observations the model was
  # fitted to.
  scale <- sqrt(apply(object[["eta"]][[1L]], 1L, stats::var) + error[2L])

  lapply(eta, function(m) (m - error[1L]) / scale[seq_len(nrow(m))])
}

# Mean and variance of the latent error the family thresholds.
#
# The mean is where the two families part company. A normal or logistic error is
# symmetric, so whether `e` or `-e` is written into `y*` makes no difference. A
# complementary log-log error is a smallest extreme value variate, which is not
# symmetric, and the two families put it in with opposite signs: the ordinal
# model has `P(Y <= k) = G(c_k - eta)`, which is `P(eta + e <= c_k)`, while the
# binomial model has `P(Y = 1) = G(eta)`, which is `P(e <= eta)` and so
# `y* = eta - e`. The variance is the same either way; the mean flips.
latent_error <- function(family, link) {
  # digamma(1) is -gamma, the mean of the smallest extreme value distribution.
  extreme_value_mean <- {
    if (identical(family, "binomial")) -digamma(1)
    else digamma(1)
  }

  switch(link,
         probit = c(0, 1),
         logit = c(0, pi^2 / 3),
         cloglog = c(extreme_value_mean, pi^2 / 6),
         arg::err("{.arg type} {.val stdlv} needs a link with a known latent
                   distribution, which {.val {link}} is not"))
}

# The mean of the response with the category labels read as numbers, which is
# what makes an ordinal outcome summarizable by a single number without
# pretending it is continuous at the modeling stage.
category_mean <- function(object, probs, values, draws) {
  levels <- dimnames(probs)[[3L]]

  arg::when_not_null(
    values,
    arg::arg_numeric
  )

  if (is_null(values)) {
    numeric_levels <- suppressWarnings(as.numeric(levels))

    if (anyNA(numeric_levels)) {
      arg::err(c("{.code type = \"mean\"} reads the response levels as numbers,
                  and {.val {levels[is.na(numeric_levels)]}} cannot be read that
                  way.",
                 i = "Give {.arg values} a named numeric vector, one element per
                      level, to say what each category is worth."))
    }

    values <- setNames(numeric_levels, levels)
  }
  else if (is_null(names(values)) || !all(levels %in% names(values))) {
    arg::err("{.arg values} must be a numeric vector named for every response
                level: {.val {levels}}")
  }

  weightsv <- values[levels]

  # One weighted sum per draw per observation.
  out <- apply(probs, c(1L, 2L), function(p) sum(p * weightsv))

  if (draws) {
    return(out)
  }

  colMeans(out)
}

# The additive predictors and the nuisance parameters at the requested draws,
# which is the front matter of every prediction and of the posterior predictive
# sampler alike.
predict_parts <- function(object, newdata = NULL, offset = NULL,
                          iterations = NULL) {

  num_save <- nrow(object[["sigma_mu"]])
  iterations <- resolve_iterations(iterations, num_save)

  eta <- {
    if (is_null(newdata)) {
      lapply(object[["eta"]], function(m) m[iterations, , drop = FALSE])
    }
    else {
      predict_eta(object, newdata, offset, iterations)
    }
  }

  aux <- {
    if (is_null(object[["aux"]])) NULL
    else object[["aux"]][iterations, , drop = FALSE]
  }

  list(eta = eta, aux = aux, iterations = iterations)
}

resolve_iterations <- function(iterations, num_save) {
  arg::when_not_null(
    iterations,
    arg::arg_whole_numeric
  )

  if (is_null(iterations)) {
    return(seq_len(num_save))
  }

  if (any(iterations < 1L) || any(iterations > num_save)) {
    arg::err("{.arg iterations} must lie between 1 and {num_save}, the number of
              stored draws")
  }

  as.integer(iterations)
}

# Rebuild the design matrix for new data exactly as it was built for fitting,
# then evaluate the stored forests.
predict_eta <- function(object, newdata, offset, iterations) {

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  tt <- stats::delete.response(object[["terms"]])
  mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                           xlev = object[["xlevels"]])

  x <- stats::model.matrix(tt, mf, contrasts.arg = object[["contrasts"]])
  assign <- attr(x, "assign")
  x <- x[, assign != 0L, drop = FALSE]

  expected <- names(object[["unit_maps"]])
  missing <- setdiff(expected, colnames(x))

  if (!is_null(missing)) {
    arg::err("{.arg newdata} does not reproduce the predictor column{?s}
              {.val {missing}}")
  }

  x <- x[, expected, drop = FALSE]

  # A missing value is only meaningful for a column that had missing values when
  # the model was fit, because only then does any splitting rule carry an answer
  # for it. Elsewhere every rule would send it the same arbitrary way.
  na_now <- colnames(x)[vapply(seq_len(ncol(x)),
                               function(j) anyNA(x[, j]), logical(1L))]
  na_then <- {
    known <- object[["has_na"]]
    if (is_null(known)) character()
    else names(known)[known]
  }
  unsupported <- setdiff(na_now, na_then)

  if (!is_null(unsupported)) {
    arg::err("{.arg newdata} has missing values in {.val {unsupported}}, which
              had none when the model was fit, so no splitting rule knows where
              to send them")
  }

  x <- apply_unit_maps(x, object[["unit_maps"]])

  eta <- .bartisan_predict(X = x,
                           forest_flat = object[["forest_flat"]],
                          tree_start = object[["tree_start"]],
                          bandwidth = object[["bandwidth"]],
                          num_forest = object[["num_forest"]],
                          num_trees = object[["num_trees"]],
                          num_save = nrow(object[["sigma_mu"]]),
                          soft = object[["soft"]],
                          gate = gate_code(object[["gate"]] %or% "logistic"),
                          iterations = as.integer(iterations) - 1L)

  intercept <- object[["intercept"]]
  n_new <- nrow(x)

  # An offset is not a function of the predictors, so it cannot be rebuilt from
  # `newdata`. Silently treating it as zero would return predictions that are
  # wrong by exactly the omitted offset, so insist on being given one.
  if (is_null(offset) && isTRUE(object[["has_offset"]])) {
    arg::err("this model was fit with an offset, so predicting new data
              requires {.arg offset} to be supplied for those rows")
  }

  user_offset <- {
    if (is_null(offset)) matrix(0, nrow = n_new, ncol = object[["num_forest"]])
    else as_offset_matrix(offset, n_new, object[["num_forest"]])
  }

  # The random part, which for a level not seen at fitting time is zero.
  ranef <- random_predict(object, newdata, iterations)

  for (h in seq_along(eta)) {
    eta[[h]] <- eta[[h]] + intercept[h] +
      rep(user_offset[, h], each = nrow(eta[[h]]))

    if (!is_null(ranef)) {
      eta[[h]] <- eta[[h]] + ranef[[h]]
    }
  }

  names(eta) <- names(object[["eta"]])

  eta
}

as_offset_matrix <- function(offset, n, n_forest) {
  arg::arg_numeric(offset)

  if (is.matrix(offset)) {
    if (nrow(offset) != n || ncol(offset) != n_forest) {
      arg::err("{.arg offset} must be a matrix with {n} row{?s} and {n_forest} column{?s}")
    }
    return(offset)
  }

  if (length(offset) != n) {
    arg::err("{.arg offset} must have one value per row of {.arg newdata}")
  }

  matrix(offset, nrow = n, ncol = n_forest)
}

shape_link <- function(eta, draws, object) {
  if (draws) {
    if (length(eta) == 1L) {
      return(eta[[1L]])
    }
    return(eta)
  }

  out <- vapply(eta, colMeans, numeric(ncol(eta[[1L]])))

  if (length(eta) == 1L) {
    return(drop(out))
  }

  colnames(out) <- names(eta)

  out
}

# The survival function at named times: the predictive distribution of a
# survival response, and the analogue of `type = "prob"` for a categorical one.
#
# Both survival families report a *point* summary through `type = "response"` --
# the median survival time -- which is not enough for the questions survival
# analysis is usually asked. The curve is what answers them, and computing it
# from the draws by hand means depending on how the baseline is stored, which is
# the package's business rather than the caller's.
survival_scale <- function(object, eta, aux, times, draws) {
  family <- object[["family"]][["family"]]
  link <- object[["family"]][["link"]]
  e <- eta[[1L]]

  labels <- format(times, trim = TRUE)
  out <- array(NA_real_, c(nrow(e), ncol(e), length(times)),
               dimnames = list(NULL, NULL, labels))

  if (identical(family, "dpm_aft")) {
    # log T = eta + W, so S(t | x) is the mixture's own survival at the
    # standardized log time. The mixture is reported centered, so the predictor
    # is already the conditional mean of log T.
    for (k in seq_along(times)) {
      for (s in seq_len(nrow(e))) {
        out[s, , k] <- dpm_survival(object, s, log(times[k]) - e[s, ])
      }
    }
  }
  else if (identical(family, "aft")) {
    # log T = eta + sigma * epsilon, so S(t) is the error's upper tail at the
    # standardized log time. `sigma` is one value per draw and `e` is draws by
    # rows, so the division recycles down the draws, which is what is wanted.
    sigma <- aux[, "sigma"]

    for (k in seq_along(times)) {
      z <- (log(times[k]) - e) / sigma
      out[, , k] <- switch(link,
                           weibull = exp(-exp(z)),
                           loglogistic = 1 / (1 + exp(z)),
                           lognormal = stats::pnorm(z, lower.tail = FALSE))
    }
  }
  else {
    # S(t) = exp(-Lambda_0(t) exp(eta)), with the cumulative baseline piecewise
    # linear in the drawn bin hazards: each bin contributes its hazard times how
    # much of it t reaches.
    edges <- object[["family_opts"]][["edges"]]
    lambda <- aux[, grep("^lambda[0-9]+$", colnames(aux)), drop = FALSE]
    span <- c(diff(edges), Inf)

    for (k in seq_along(times)) {
      reach <- pmin(pmax(times[k] - edges, 0), span)
      out[, , k] <- exp(-drop(lambda %*% reach) * exp(e))
    }
  }

  if (draws) {
    return(out)
  }

  averaged <- apply(out, c(2L, 3L), mean)
  dimnames(averaged) <- list(NULL, labels)
  averaged
}

# The mean of the response, except for the survival families, where the median
# survival time is the interpretable summary and depends on the scale parameter.
response_scale <- function(object, eta, aux, draws) {
  family <- object[["family"]][["family"]]
  link <- object[["family"]][["link"]]
  e <- eta[[1L]]

  # A link the engine does not carry natively is composed onto the engine's own
  # scale when fitting, so the predictor here is on the caller's scale and the
  # caller's inverse link is what maps it to the mean.
  inv <- object[["family"]][["custom_link"]][["linkinv"]]
  supplied <- is.function(inv)

  out <- switch(family,
                # The error distribution is estimated and is not centered, so the
                # conditional mean is the predictor plus whatever mean the
                # mixture came out with. `type = "link"` is the sum of trees on
                # its own.
                # The mixture is reported centered, so the predictor already
                # is the conditional mean and the response scale is the identity.
                dpm = e,
                gaussian = ,
                location_scale = if (supplied) inv(e) else e,
                # The beta mean is the inverse link of the predictor, the same
                # as a binomial's, and on the same three links.
                beta = ,
                binomial = if (supplied) inv(e) else binomial_linkinv(e, link),
                poisson = ,
                negbin = ,
                Gamma = if (supplied) inv(e) else exp(e),
                # Solve Lambda_0(t) exp(eta) = log 2 for t. The cumulative
                # baseline is piecewise linear in the drawn bin hazards, so this
                # is an exact inversion rather than a search.
                ph = {
                  edges <- object[["family_opts"]][["edges"]]
                  cols <- grep("^lambda[0-9]+$", colnames(aux))
                  span <- c(diff(edges), Inf)

                  out <- vapply(seq_len(nrow(e)), function(s) {
                    lam <- aux[s, cols]
                    # Cumulative baseline at each edge, then locate log2 / exp(eta).
                    at_edge <- c(0, cumsum(lam[-length(lam)] *
                                             span[-length(span)]))
                    target <- log(2) / exp(e[s, ])
                    b <- pmax(findInterval(target, at_edge), 1L)
                    edges[b] + (target - at_edge[b]) / lam[b]
                  }, numeric(ncol(e)))

                  t(out)
                },
                # The median survival time, as for the other survival
                # families. The mixture's own median is found by bisection on its
                # survival function, since a mixture has no closed-form quantile.
                #
                # The median being bisected for is a property of the *error*
                # distribution, which is shared by every observation: only the
                # shift `e[s, ]` differs between them. So this is one scalar per
                # draw, not one per observation. Bisecting a vector of identical
                # values instead cost a factor of n in `pnorm` calls -- at 700
                # observations, 588 million of them to produce 500 numbers, which
                # was six times the cost of the sampler itself.
                dpm_aft = {
                  median_error <- vapply(seq_len(nrow(e)), function(s) {
                    lo <- -30
                    hi <- 30
                    comp <- mixture_at(object, s)

                    for (step in seq_len(60L)) {
                      mid <- (lo + hi) / 2

                      if (dpm_survival(object, s, mid, comp) > 0.5) {
                        lo <- mid
                      }
                      else {
                        hi <- mid
                      }
                    }

                    (lo + hi) / 2
                  }, numeric(1))

                  # One median per draw, added down the rows of `e`, which is
                  # draws by observations.
                  exp(e + median_error)
                },
                aft = {
                  sigma <- aux[, "sigma"]
                  # log T = eta + sigma * epsilon, so the median is eta shifted by the
                  # median of the error, which is zero except for the Weibull case.
                  shift <- switch(link, weibull = log(log(2)), 0)
                  exp(e + shift * sigma)
                },
                # The mean of a zero-inflated count is the count mean scaled by the
                # probability that the observation is not a structural zero.
                zip = ,
                zinb = exp(e) * (1 - stats::plogis(eta[[2L]])),
                # The mean of an ordered beta response mixes the two point masses with the
                # beta mean on the interior.
                ordbeta = {
                  p_upper <- stats::plogis(e - aux[, "cut2"])
                  p_middle <- stats::plogis(e - aux[, "cut1"]) - p_upper
                  p_middle * stats::plogis(e) + p_upper
                },
                # For a response with categories there is no single mean to report, so the
                # response scale is the vector of category probabilities, as it is for
                # predict.polr() and predict.multinom().
                ordinal = ,
                mnp = ,
                multinomial = {
                  probs <- category_probs(object, eta, aux)

                  if (draws) {
                    return(probs)
                  }

                  out <- apply(probs, c(2L, 3L), mean)
                  dimnames(out) <- list(NULL, dimnames(probs)[[3L]])
                  return(out)
                },
                e)

  if (draws) {
    return(out)
  }

  colMeans(out)
}

# The covariance columns of `aux`, which for a multinomial probit fit are the
# lower triangle of the latent covariance matrix. A two-category fit has none:
# the trace constraint pins the single variance at one, so there is nothing to
# record and nothing to read back.
mnp_sigma <- function(object, aux) {
  latent <- object[["num_cat"]] - 1L
  wanted <- latent * (latent + 1L) / 2L

  if (latent == 1L) {
    return(matrix(1, nrow = nrow(object[["eta"]][[1L]]), ncol = 1L))
  }

  columns <- grep("^sigma[0-9]+$", colnames(aux), value = TRUE)

  if (length(columns) != wanted) {
    arg::err("the fit does not carry the {wanted} covariance element{?s} that a
              {object[['num_cat']]}-category probit model needs")
  }

  aux[, columns, drop = FALSE]
}

# Category probabilities, returned as draws by observations by categories.
category_probs <- function(object, eta, aux) {
  family <- object[["family"]][["family"]]
  link <- object[["family"]][["link"]]

  if (identical(family, "binomial")) {
    inv <- object[["family"]][["custom_link"]][["linkinv"]]
    p <- {
      if (is.function(inv)) inv(eta[[1L]])
      else binomial_linkinv(eta[[1L]], link)
    }
    out <- array(c(1 - p, p), dim = c(dim(p), 2L))
    levels <- object[["levels"]] %or% c("0", "1")
    dimnames(out) <- list(NULL, NULL, levels)
    return(out)
  }

  if (identical(family, "ordinal")) {
    e <- eta[[1L]]
    num_cat <- object[["num_cat"]]
    cuts <- {
      if (num_cat > 2L) aux[, sprintf("cut%d", seq_len(num_cat - 1L)),
                            drop = FALSE]
      else matrix(0, nrow = nrow(e), ncol = 1L)
    }

    cdf <- switch(link,
                  probit = stats::pnorm,
                  cloglog = function(z) -expm1(-exp(z)),
                  stats::plogis)

    out <- array(0, dim = c(dim(e), num_cat))
    previous <- array(0, dim = dim(e))

    for (k in seq_len(num_cat)) {
      current <- {
        if (k == num_cat) array(1, dim = dim(e))
        else cdf(cuts[, k] - e)
      }
      out[, , k] <- pmax(current - previous, 0)
      previous <- current
    }

    dimnames(out) <- list(NULL, NULL, object[["levels"]])
    return(out)
  }

  # Multinomial probit: the probability that the argmax of a correlated Gaussian
  # vector falls in each category, which is an orthant probability with no closed
  # form and so is simulated. See the Details of this page.
  if (identical(family, "mnp")) {
    out <- .bartisan_mnp_probs(eta, mnp_sigma(object, aux),
                               object[["family_opts"]][["replicates"]])
    dimnames(out) <- list(NULL, NULL, object[["levels"]])
    return(out)
  }

  # Multinomial. Computed by shifting out the largest predictor before
  # exponentiating, so that a large predictor cannot overflow the way a direct
  # exp() would. Under the symmetric coding every category has a predictor;
  # under reference coding the first category's is fixed at zero.
  num_cat <- object[["num_cat"]]
  symmetric <- isTRUE(object[["family_opts"]][["symmetric"]])
  dims <- dim(eta[[1L]])
  out <- array(0, dim = c(dims, num_cat))

  shift <- {
    if (symmetric) eta[[1L]]
    else array(0, dim = dims)
  }

  for (h in seq_along(eta)) {
    shift <- pmax(shift, eta[[h]])
  }

  denominator <- {
    if (symmetric) array(0, dim = dims)
    else exp(-shift)
  }

  for (h in seq_along(eta)) {
    denominator <- denominator + exp(eta[[h]] - shift)
  }

  offset <- as.integer(!symmetric)

  if (!symmetric) {
    out[, , 1L] <- exp(-shift) / denominator
  }

  for (h in seq_along(eta)) {
    out[, , h + offset] <- exp(eta[[h]] - shift) / denominator
  }

  dimnames(out) <- list(NULL, NULL, object[["levels"]])

  out
}

# Posterior mean of the response-scale fit, stored on the object by bartisan().
fitted_from_eta <- function(object, eta, average = TRUE) {
  family <- object[["family"]][["family"]]

  if (family %in% c("ordinal", "multinomial", "mnp")) {
    probs <- category_probs(object, eta, object[["aux"]])
    out <- apply(probs, c(2L, 3L), mean)
    dimnames(out) <- list(NULL, dimnames(probs)[[3L]])
    return(out)
  }

  response_scale(object, eta, object[["aux"]], draws = !average)
}

# Conditional density of the outcome given the predictors, at each posterior
# draw. The family's own C++ log density is reused rather than reimplemented in
# R, so there is one definition of each distribution in the package.
conditional_density <- function(object, newdata, eta, aux, weights, draws, log,
                                iterations) {

  # A multinomial probit likelihood is a Gaussian orthant probability with no
  # closed form, so there is no C++ log density to reuse: the probability of the
  # observed category comes from the same simulator predictions use.
  if (identical(object[["family"]][["family"]], "mnp")) {
    return(mnp_density(object, newdata, eta, aux, draws, log))
  }

  # A Dirichlet process mixture's likelihood is the mixture's own, which lives in
  # the drawn components rather than in a fixed set of nuisance parameters, so it
  # is not reachable through the engine's log density either.
  if (identical(object[["family"]][["family"]], "dpm")) {
    return(dpm_density(object, newdata, eta, iterations, draws, log))
  }

  if (identical(object[["family"]][["family"]], "dpm_aft")) {
    return(dpm_aft_density(object, newdata, eta, iterations, draws, log))
  }

  parts <- {
    if (is_null(newdata)) {
      list(y = object[["y"]], weights = weights %or% object[["prior_weights"]],
           opts = object[["family_opts"]])
    }
    else {
      density_response(object, newdata, weights)
    }
  }

  num_draws <- nrow(eta[[1L]])

  aux_matrix <- {
    if (is_null(aux)) matrix(0, nrow = num_draws, ncol = 0L)
    else as.matrix(aux)
  }

  out <- .bartisan_logdens(y = parts[["y"]],
                           weights = parts[["weights"]],
                          eta_draws = eta,
                          family_name = object[["family"]][["family"]],
                          link = object[["family"]][["link"]],
                          family_opts = parts[["opts"]],
                          aux = aux_matrix)

  warn_undefined_density(out, object)

  if (draws) {
    return(if (log) out else exp(out))
  }

  # Averaging is done on the density scale, so the summarized value is the
  # pointwise predictive density: the log of it is what gets summed for a
  # held-out log score, and `log` stays a pure transform of the returned value.
  # The averaging is done through log-sum-exp because individual densities
  # underflow readily.
  num_draws <- nrow(out)
  lppd <- apply(out, 2L, function(z) {
    m <- max(z)
    if (!is.finite(m)) {
      return(m)
    }
    m + base::log(sum(exp(z - m))) - base::log(num_draws)
  })

  if (log) lppd else exp(lppd)
}

# A composed link whose inverse does not cover the whole line can send a draw's
# parameter outside the family's support, where the density is undefined rather
# than small. `bartisan()` warns about the link at fit time and rejects such
# proposals while sampling, but a *saved* draw can still land out of support at a
# new `x`, because the forest extrapolates there and was never asked to.
#
# The value stays NaN rather than becoming zero. Zero would assert that the
# outcome is impossible, which is a stronger claim than the model makes and a
# worse one: it turns into `-Inf` in a log score, which reads as a legitimately
# terrible fit rather than an undefined one. NaN propagates visibly.
#
# What is worth reporting is the amplification. The draws are averaged before the
# log is taken, so a single undefined draw out of hundreds makes that
# observation's density NaN -- a fraction of a percent of draws routinely
# accounts for a third of the returned values.
warn_undefined_density <- function(out, object) {
  bad <- is.na(out)

  if (!any(bad)) {
    return(invisible(NULL))
  }

  n_draws <- sum(bad)
  n_obs <- sum(colSums(bad) > 0)
  total_obs <- ncol(out)
  link <- object[["family"]][["link"]]

  arg::wrn(c("the conditional density is undefined for {n_draws} of
              {length(bad)} draw-by-observation values, which makes {n_obs} of
              {total_obs} returned {cli::qty(total_obs)}value{?s} {.code NaN}.",
             i = "The {.val {link}} link's inverse does not cover the whole
                  additive predictor, so at these predictors some draws imply a
                  parameter outside the family's support. Draws are averaged
                  before the log is taken, so one undefined draw is enough to
                  make an observation {.code NaN}.",
             i = "A link whose inverse is defined on the whole line, such as
                  {.val log}, avoids this. {.code type = \"link\"} and
                  {.code draws = TRUE} show which predictors are responsible."))

  invisible(NULL)
}

# The simulated probability of each observation's own category. Simulation error
# is per draw and averages down over draws, which is why `draws = FALSE` is the
# more accurate of the two.
mnp_density <- function(object, newdata, eta, aux, draws, log) {
  probs <- category_probs(object, eta, aux)
  levels <- object[["levels"]]

  observed <- {
    if (is_null(newdata)) as.integer(object[["y"]]) + 1L
    else match_levels(density_response_vector(object, newdata), levels) + 1L
  }

  if (length(observed) != dim(probs)[2L]) {
    arg::err("{.arg newdata} and the outcome disagree about the number of
              observations")
  }

  out <- vapply(seq_along(observed),
                function(i) probs[, i, observed[i]],
                numeric(dim(probs)[1L]))

  if (!draws) {
    out <- colMeans(out)
  }

  if (log) base::log(out) else out
}

# The mixture components of one draw, as a matrix of mean, standard deviation and
# weight. The count differs from draw to draw, which is why they are stored flat
# with an offset per draw rather than in a matrix.
mixture_at <- function(object, iteration) {
  start <- object[["mixture_start"]]

  if (is_null(start)) {
    arg::err("this fit carries no mixture; only {.fn dpm} and {.fn dpm_aft}
              have one")
  }

  at <- seq_len(start[iteration + 1L] - start[iteration]) + start[iteration]

  matrix(object[["mixture_flat"]][at], ncol = 3L, byrow = TRUE,
         dimnames = list(NULL, c("mean", "sd", "weight")))
}

# The predictive density of an error under one draw: the occupied components
# weighted by their sizes, plus the chance under the Dirichlet process that a new
# observation opens a component of its own, which is the baseline's own marginal.
# This is the same expression the sampler reports its log likelihood from.
dpm_predictive <- function(object, iteration, at) {
  opts <- object[["family_opts"]]
  n <- object[["n"]]
  alpha <- object[["aux"]][iteration, "alpha"]
  total <- alpha + n

  components <- mixture_at(object, iteration)

  # The occupied components are reported centered, but the baseline a fresh
  # component would be drawn from is centered on the raw chart, so it is
  # evaluated `center` further along. `center` is the shift that was taken out.
  center <- object[["aux"]][iteration, "center"]
  scale <- sqrt(opts[["lambda"]] * (1 + 1 / opts[["k_0"]]))
  out <- alpha / total *
    stats::dt((at + center) / scale, df = opts[["nu"]]) / scale

  for (k in seq_len(nrow(components))) {
    out <- out + components[k, "weight"] * n / total *
      stats::dnorm(at, components[k, "mean"], components[k, "sd"])
  }

  out
}

# The survival function of one error under the current draw, matching
# `dpm_predictive()` term for term. This is what a censored observation
# contributes, and what turns the mixture into a survival curve.
dpm_survival <- function(object, iteration, at, components = NULL) {
  opts <- object[["family_opts"]]
  n <- object[["n"]]
  alpha <- object[["aux"]][iteration, "alpha"]
  total <- alpha + n

  # `components` is accepted already built for callers that evaluate the same
  # draw repeatedly -- the median bisection does it sixty times -- since
  # rebuilding the matrix each time costs more than the arithmetic does.
  if (is_null(components)) {
    components <- mixture_at(object, iteration)
  }

  center <- object[["aux"]][iteration, "center"]
  scale <- sqrt(opts[["lambda"]] * (1 + 1 / opts[["k_0"]]))

  out <- alpha / total *
    stats::pt((at + center) / scale, df = opts[["nu"]], lower.tail = FALSE)

  for (k in seq_len(nrow(components))) {
    out <- out + components[k, "weight"] * n / total *
      stats::pnorm(at, components[k, "mean"], components[k, "sd"],
                   lower.tail = FALSE)
  }

  out
}

# The predictive density of each observation's own outcome.
dpm_density <- function(object, newdata, eta, iterations, draws, log) {
  y <- {
    if (is_null(newdata)) object[["y"]]
    else check_numeric_response(density_response_vector(object, newdata), "dpm")
  }

  predictor <- eta[[1L]]

  if (length(y) != ncol(predictor)) {
    arg::err("{.arg newdata} and the outcome disagree about the number of
              observations")
  }

  out <- vapply(seq_len(nrow(predictor)), function(s) {
    dpm_predictive(object, iterations[s], y - predictor[s, ])
  }, numeric(length(y))) |>
    t()

  if (!draws) {
    out <- colMeans(out)
  }

  if (log) base::log(out) else out
}

# The observed-data likelihood of a `dpm_aft()` fit: the mixture's density at an
# observed event, its survival at a censoring. The complete-data version the
# sampler imputes with does not belong here, because it would credit a censored
# observation with a failure it did not have.
dpm_aft_density <- function(object, newdata, eta, iterations, draws, log) {
  surv <- {
    if (is_null(newdata)) {
      list(log_time = object[["y"]],
           event = object[["family_opts"]][["event"]])
    }
    else {
      prepare_surv(density_response_vector(object, newdata), nrow(newdata))
    }
  }

  predictor <- eta[[1L]]

  if (length(surv[["log_time"]]) != ncol(predictor)) {
    arg::err("{.arg newdata} and the outcome disagree about the number of
              observations")
  }

  seen <- surv[["event"]] > 0

  out <- vapply(seq_len(nrow(predictor)), function(s) {
    at <- surv[["log_time"]] - predictor[s, ]
    value <- dpm_survival(object, iterations[s], at)
    value[seen] <- dpm_predictive(object, iterations[s], at[seen])
    value
  }, numeric(length(surv[["log_time"]]))) |>
    t()

  if (!draws) {
    out <- colMeans(out)
  }

  if (log) base::log(out) else out
}

#' Error distribution of a Dirichlet process mixture fit
#'
#' @description
#' The estimated density of the errors of a [dpm()] fit, which is the object the
#' method exists to produce: BART commits to one normal, and this says what
#' shape the errors actually have. It is the posterior of the density of a new
#' error, evaluated on a grid.
#'
#' @param object a fitted model from [bartisan()] with `family = dpm()`.
#' @param at the grid to evaluate on. Defaults to 201 points spanning four
#'   posterior-mean error standard deviations either side of zero.
#' @param level width of the pointwise interval.
#' @param iterations optional integer vector selecting which stored draws to
#'   use. Defaults to all of them.
#'
#' @returns
#' A data frame with one row per grid point and columns `at`, `mean`, `lower` and
#' `upper` -- the posterior mean density and a pointwise interval, ready to plot
#' against the normal density a Gaussian fit would have assumed.
#'
#' @seealso [dpm()], [bartisan()]
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' d <- data.frame(x = runif(n, -1, 1))
#' d$y <- 10 * d$x^3 + rt(n, 3)
#'
#' fit <- bartisan(y ~ x, data = d, family = dpm(),
#'                control = bartisan_control(num_trees = 20, num_burn = 100,
#'                                          num_save = 100, verbose = FALSE))
#'
#' density <- error_density(fit)
#' plot(density$at, density$mean, type = "l", xlab = "error", ylab = "density")
#' lines(density$at, density$lower, lty = 2)
#' lines(density$at, density$upper, lty = 2)
#'
#' @export
error_density <- function(object, at = NULL, level = 0.95,
                          iterations = NULL) {

  if (!inherits(object, "bartisan_fit")) {
    arg::err("{.arg object} must be a fit from {.fn bartisan}")
  }

  if (!object[["family"]][["family"]] %in% c("dpm", "dpm_aft")) {
    arg::err(c("{.fn error_density} needs a fit with an estimated error
                distribution, which is what {.fn dpm} has.",
               i = "Every other family fixes the error distribution, so its
                    density is a closed form rather than something to estimate."))
  }

  arg::arg_number(level)
  arg::arg_between(level, c(0, 1))

  iterations <- resolve_iterations(iterations, nrow(object[["aux"]]))

  if (is_null(at)) {
    reach <- 4 * mean(object[["aux"]][iterations, "error_sd"])
    at <- seq(-reach, reach, length.out = 201L)
  }
  else {
    arg::arg_numeric(at)
  }

  # `matrix()` rather than the shape vapply() returns, so that a single grid
  # point is still a one-row matrix rather than a vector.
  drawn <- matrix(vapply(iterations,
                         function(s) dpm_predictive(object, s, at),
                         numeric(length(at))),
                  nrow = length(at))

  summarized <- apply(drawn, 1L, post_summary, level = level)

  summarized <- matrix(summarized, nrow = 4L,
                       dimnames = list(c("mean", "sd", "lower", "upper"), NULL))

  data.frame(at = at,
             mean = summarized["mean", ],
             lower = summarized["lower", ],
             upper = summarized["upper", ])
}

# The response column of `newdata`, for the families whose density is computed in
# R rather than through the engine.
density_response_vector <- function(object, newdata) {
  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  absent <- setdiff(all.vars(object[["terms"]][[2L]]), names(newdata))

  if (!is_null(absent)) {
    arg::err("{.arg type} {.val density} is the density of the outcome, so
              {.arg newdata} must contain the outcome variable{?s}
              {.val {absent}}")
  }

  mf <- stats::model.frame(object[["terms"]], newdata,
                           na.action = stats::na.pass,
                          xlev = object[["xlevels"]])

  stats::model.response(mf, "any")
}

# Coerce the outcome in `newdata` the same way it was coerced at fitting time,
# reusing the levels recorded then rather than whatever happens to appear in the
# new data.
density_response <- function(object, newdata, weights) {

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  response_vars <- all.vars(object[["terms"]][[2L]])
  absent <- setdiff(response_vars, names(newdata))

  if (!is_null(absent)) {
    arg::err("{.arg type} {.val density} is the density of the outcome, so
              {.arg newdata} must contain the outcome variable{?s}
              {.val {absent}}")
  }

  mf <- stats::model.frame(object[["terms"]], newdata,
                           na.action = stats::na.pass,
                           xlev = object[["xlevels"]])
  y <- stats::model.response(mf, "any")

  n <- nrow(mf)
  family <- object[["family"]][["family"]]
  levels <- object[["levels"]]
  opts <- object[["family_opts"]]
  w <- weights %or% rep.int(1, n)

  if (length(w) != n) {
    arg::err("{.arg weights} must have one value per row of {.arg newdata}")
  }

  y_out <- switch(family,
                  gaussian = ,
                  location_scale = check_numeric_response(y, family),
                  Gamma = check_numeric_response(y, family),
                  poisson = ,
                  negbin = ,
                  zip = ,
                  zinb = check_count_response(y, family),
                  beta = ,
                  ordbeta = check_numeric_response(y, family),
                  binomial = {
                    b <- prepare_binomial_levels(y, w, levels)
                    w <- b[["weights"]]
                    b[["y"]]
                  },
                  ordinal = ,
                  multinomial = match_levels(y, levels),
                  aft = {
                    a <- prepare_surv(y, n)
                    # The event indicator is part of the outcome, so it comes from the new
                    # data rather than from the stored options.
                    opts[["event"]] <- a[["event"]]
                    a[["log_time"]]
                  },
                  dpm_aft = {
                    a <- prepare_surv(y, n)
                    opts[["event"]] <- a[["event"]]
                    a[["log_time"]]
                  },
                  ph = {
                    a <- prepare_surv(y, n)
                    opts[["event"]] <- a[["event"]]
                    # The bin edges are structure fitted to the training times, so
                    # they stay as they were; only the times and events are new.
                    a[["time"]]
                  },
                  check_numeric_response(y, family))

  list(y = y_out, weights = w, opts = opts)
}

# The binomial coercion of prepare_binomial(), but mapping a factor through the
# levels seen at fitting time.
prepare_binomial_levels <- function(y, weights, levels) {
  if (is.matrix(y) && ncol(y) == 2L) {
    trials <- rowSums(y)

    if (any(trials <= 0)) {
      arg::err("every row of a two-column binomial response must have at least
                one trial")
    }

    return(list(y = y[, 1L] / trials, weights = weights * trials))
  }

  if (is.matrix(y) && ncol(y) == 1L) {
    y <- drop(y)
  }

  if (is.factor(y) || is.character(y)) {
    if (is_null(levels)) {
      arg::err("the model was not fit to a factor response, so a factor cannot
                be scored against it")
    }
    codes <- match(as.character(y), levels)

    if (anyNA(codes)) {
      arg::err("the outcome has values the model was not fit with:
                {.val {setdiff(as.character(y), levels)}}")
    }

    return(list(y = as.numeric(codes > 1L), weights = weights))
  }

  if (is.logical(y)) {
    return(list(y = as.numeric(y), weights = weights))
  }

  list(y = check_numeric_response(y, "binomial"), weights = weights)
}

match_levels <- function(y, levels) {
  if (is_null(levels)) {
    arg::err("the model records no response categories")
  }

  codes <- match(as.character(y), levels)

  if (anyNA(codes)) {
    arg::err("the outcome has categories the model was not fit with:
              {.val {setdiff(as.character(y), levels)}}")
  }

  as.numeric(codes - 1L)
}
