#' Predictions from a generalized BART model
#'
#' @description
#' Evaluates the stored posterior draws of the forests, either at the data used
#' to fit the model or at new data. Because every draw of every tree is kept,
#' predictions carry full posterior uncertainty rather than being a single point
#' estimate.
#'
#' @param object a fitted model from [genbart()].
#' @param newdata optional data frame at which to predict. Omit it to use the
#'   data the model was fit to. Missing predictor values are allowed in the
#'   columns that had them when the model was fit, since only those columns'
#'   splitting rules carry an answer for one; a missing value anywhere else is an
#'   error. See the Missing predictor values section of [genbart()].
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
#'       needs no assumption at the modelling stage: the model is still ordinal
#'       and only the reporting treats the categories as numbers. Use `values` to
#'       say what the categories are worth when the labels are not numbers, or
#'       are not the numbers you mean.}
#'     \item{`"stdlv"`}{the additive predictor divided by the standard deviation
#'       of the latent variable it indexes, for the ordinal family. An ordinal
#'       model can be written as a continuous `y* = eta + e` cut at the
#'       cutpoints, and the link fixes the distribution of `e` and so its
#'       variance; dividing by the standard deviation of `y*` puts fits with
#'       different links, or different amounts of signal, on one scale, which is
#'       what a standardized effect size on an ordinal outcome needs. See
#'       Details.}
#'     \item{`"density"`}{the conditional density of the outcome given the
#'       predictors, evaluated at the observed outcome. This requires the
#'       outcome, so `newdata` must contain it; omit `newdata` to use the data
#'       the model was fit to. The value is the likelihood contribution of the
#'       observation, so it is a density for a continuous response, a
#'       probability for a discrete one, and a survival probability for a
#'       censored survival time. Useful for held-out log scores and for
#'       posterior predictive checks.}
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
#' # The standardized latent variable
#'
#' `type = "stdlv"` reports `(eta - E[e]) / sd(y*)` for the latent
#' `y* = eta + e`, following
#' [WeightIt::predict.ordinal_weightit()]. Two parts of that need saying.
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
#' complementary log-log link, whose error is a smallest extreme value variate
#' with mean `-gamma`.
#'
#' An ordinal model is identified only up to a common shift of its cutpoints and
#' its predictor, so the location of this quantity is a convention rather than a
#' fact, and the one used here is the same one the cutpoints use -- a predictor
#' centered over the fitted sample. Against `WeightIt`, which identifies by
#' dropping the intercept column instead, the two agree on the scale and differ by
#' a constant; measured on a linear truth, the standard deviations agree to under
#' 1% and the difference is constant to three decimals. Differences on this scale,
#' which is what a standardized quantity is for, are unaffected.
#'
#' @seealso [genbart()]
#'
#' @examples
#' set.seed(1)
#'
#' n <- 150
#' d <- data.frame(x1 = runif(n), x2 = runif(n))
#' d$y <- rpois(n, exp(0.5 + d$x1))
#'
#' fit <- genbart(y ~ x1 + x2, data = d, family = poisson(),
#'                control = genbart_control(num_trees = 10, num_burn = 50,
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
predict.genbart <- function(object, newdata = NULL, type = "response",
                            draws = FALSE, iterations = NULL, offset = NULL,
                            weights = NULL, values = NULL, log = FALSE, ...) {

  type <- arg::match_arg(type, c("link", "response", "prob", "class",
                                 "mean", "stdlv", "density"))
  arg::arg_flag(draws)
  arg::arg_flag(log)

  family <- object[["family"]][["family"]]
  categorical <- family %in% c("binomial", "ordinal", "multinomial")

  if (type %in% c("prob", "class", "mean") && !categorical) {
    arg::err("{.arg type} {.val {type}} is available only for the
              {.val binomial}, {.val ordinal} and {.val multinomial} families")
  }

  if (identical(type, "stdlv") && !identical(family, "ordinal")) {
    arg::err("{.arg type} {.val stdlv} is available only for the
              {.val ordinal} family, which is the one with a latent variable to
              standardize")
  }

  if (!identical(type, "mean") && !is_null(values)) {
    cli::cli_warn("{.arg values} is ignored unless
                   {.code type = \"mean\"}")
  }

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
    return(conditional_density(object, newdata, eta, aux, weights, draws, log))
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
# An ordinal model can be written as a latent continuous response, `y* = eta +
# e`, cut at the cutpoints; the link fixes the distribution of `e` and therefore
# its variance, which is what makes the latent scale identified at all. Dividing
# by the standard deviation of `y*` puts models with different links, or with
# different amounts of signal, on one scale -- which is what a standardized
# effect size on an ordinal outcome needs. This follows
# `WeightIt::predict.ordinal_weightit(type = "stdlv")`.
#
# Two details matter. The variance of `eta` is taken over the **fitted** sample,
# so it is a property of the model rather than of whatever is being predicted;
# and `e` is shifted to have mean zero, which moves its mean into the index. That
# second point is invisible for the logit and probit links, where the error is
# already centered, and is the whole of the difference for the complementary
# log-log link, whose error has mean `-gamma`.
standardized_latent <- function(object, eta) {
  link <- object[["family"]][["link"]]

  # Variance and mean of the latent error implied by the link.
  error <- switch(link,
                  probit = c(0, 1),
                  logit = c(0, pi^2 / 3),
                  # digamma(1) is -gamma, the mean of the smallest
                  # extreme value distribution.
                  cloglog = c(digamma(1), pi^2 / 6),
                  arg::err("{.arg type} {.val stdlv} needs a link with a known
                            latent distribution, which {.val {link}} is not"))

  # The predictor's own variance, per draw, over the observations the model was
  # fitted to.
  fitted_eta <- object[["eta"]][[1L]]
  spread <- apply(fitted_eta, 1L, stats::var)
  scale <- sqrt(spread + error[2L])

  lapply(eta, function(m) (m - error[1L]) / scale[seq_len(nrow(m))])
}

# The mean of the response with the category labels read as numbers, which is
# what makes an ordinal outcome summarizable by a single number without
# pretending it is continuous at the modelling stage.
category_mean <- function(object, probs, values, draws) {
  levels <- dimnames(probs)[[3L]]

  if (is_null(values)) {
    numeric_levels <- suppressWarnings(as.numeric(levels))

    if (anyNA(numeric_levels)) {
      arg::err(c("{.arg type} {.val mean} reads the response levels as numbers,
                  and {.val {levels[is.na(numeric_levels)]}} cannot be read that
                  way",
                 i = "give {.arg values} a named numeric vector, one element per
                      level, to say what each category is worth"))
    }

    values <- stats::setNames(numeric_levels, levels)
  }
  else {
    arg::arg_numeric(values)

    if (is_null(names(values)) || !all(levels %in% names(values))) {
      arg::err("{.arg values} must be a numeric vector named for every response
                level: {.val {levels}}")
    }
  }

  weightsv <- values[levels]

  # One weighted sum per draw per observation.
  out <- apply(probs, c(1L, 2L), function(p) sum(p * weightsv))

  if (draws) {
    return(out)
  }

  colMeans(out)
}

resolve_iterations <- function(iterations, num_save) {
  if (is_null(iterations)) {
    return(seq_len(num_save))
  }

  arg::arg_whole_numeric(iterations)

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

  if (length(missing) > 0L) {
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
    if (is_null(known)) character(0L)
    else names(known)[known]
  }
  unsupported <- setdiff(na_now, na_then)

  if (length(unsupported) > 0L) {
    arg::err("{.arg newdata} has missing values in {.val {unsupported}}, which
              had none when the model was fit, so no splitting rule knows where
              to send them")
  }

  x <- apply_unit_maps(x, object[["unit_maps"]])

  eta <- .genbart_predict(X = x,
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
      arg::err("{.arg offset} must be a matrix with {n} rows and {n_forest}
                columns")
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
    gaussian = ,
    location_scale = if (supplied) inv(e) else e,
    binomial = if (supplied) inv(e) else binomial_linkinv(e, link),
    poisson = ,
    negbin = ,
    Gamma = if (supplied) inv(e) else exp(e),
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
      if (num_cat > 2L) aux[, paste0("cut", seq_len(num_cat - 1L)),
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

# Posterior mean of the response-scale fit, stored on the object by genbart().
fitted_from_eta <- function(object, eta, average = TRUE) {
  family <- object[["family"]][["family"]]

  if (family %in% c("ordinal", "multinomial")) {
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
conditional_density <- function(object, newdata, eta, aux, weights, draws,
                                log) {

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

  out <- .genbart_logdens(y = parts[["y"]],
                          weights = parts[["weights"]],
                          eta_draws = eta,
                          family_name = object[["family"]][["family"]],
                          link = object[["family"]][["link"]],
                          family_opts = parts[["opts"]],
                          aux = aux_matrix)

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

# Coerce the outcome in `newdata` the same way it was coerced at fitting time,
# reusing the levels recorded then rather than whatever happens to appear in the
# new data.
density_response <- function(object, newdata, weights) {

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  response_vars <- all.vars(object[["terms"]][[2L]])
  absent <- setdiff(response_vars, names(newdata))

  if (length(absent) > 0L) {
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
  w <- weights %or% rep(1, n)

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
