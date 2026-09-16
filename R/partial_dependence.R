#' Partial dependence on one or two predictors
#'
#' Averages the fitted surface over the sample at each value of one or two
#' predictors, so that what is left is how the prediction moves with them.
#' `partial_dependence()` returns the values and `plot()` draws them.
#'
#' @param object a `<bartisan_fit>` object; the output of a call to [bartisan()].
#' @param variables a one-sided formula naming one or two predictors, as in
#'   `~ age` or `~ age + sex`, or a character vector of their names.
#' @param newdata optional; a data frame to average over. Default is the data
#'   the model was fit to.
#' @param grid `integer`; how many values of a numeric predictor to evaluate.
#'   Default is 51. A factor is evaluated at each of its levels whatever this is.
#' @param values optional; a named list giving the values to evaluate a
#'   predictor at, which overrides `grid` for the predictors it names.
#' @param level `numeric`; the level of the credible interval. Default is `.95`.
#' @param type `string`; the prediction scale, passed to
#'   [predict.bartisan_fit()]. Default is `"response"`.
#' @param plot `logical`; whether to draw the result rather than return it.
#'   Default is `FALSE`. `plot = TRUE` calls [plot.bartisan_partial()], so the
#'   argument and the method cannot disagree.
#' @param x a `<bartisan_partial>` object; the output of a call to
#'   `partial_dependence()`. For `plot.bartisan_fit()`, a `<bartisan_fit>`.
#' @param digits `integer`; for `print()`, the number of significant digits to
#'   print the estimates and their interval to. Default is 3.
#' @param y for `plot.bartisan_fit()`, the predictors to plot, as `variables`
#'   above.
#' @param ... for `plot.bartisan_fit()`, further arguments passed to
#'   `partial_dependence()`; for `partial_dependence()`, further arguments passed to [predict.bartisan_fit()].
#'
#' @returns
#' A `<bartisan_partial>` object, a data frame with one row per grid point,
#' columns for the predictors, and `estimate`, `lower` and `upper`. `plot()`
#' returns a \pkg{ggplot2} object.
#'
#' @details
#' At each grid value every unit is assigned that value, the prediction is taken
#' for all of them, and the average over units is taken *within each posterior
#' draw*. The interval is then a quantile of those averages, so it is an interval
#' on the average prediction and not on any one unit's.
#'
#' What this is and is not worth reading as a description of the fit is the same
#' caveat that applies to any partial dependence plot. Averaging over the other
#' predictors evaluates the model at covariate combinations that may not occur,
#' so a curve over a region where the predictor's values are sparse says more
#' about the prior than about the data. It is a summary of the fitted function
#' rather than a causal claim; for a contrast that is meant causally, see
#' [estimate_effect()].
#'
#' @seealso [estimate_effect()] for treatment effects rather than fitted
#'   surfaces; [variable_importance()] for which predictors the forest uses;
#'   [bartisan-marginaleffects], since
#'   \pkgfun{marginaleffects}{plot_predictions} draws the same thing with more
#'   control over the grid
#'
#' @examples
#' data("rhc")
#' set.seed(123)
#'
#' fit <- bartisan(death ~ age + sex + meanbp + aps, data = rhc,
#'                 num_trees = 10, num_burn = 50, num_draws = 50,
#'                 verbose = FALSE)
#'
#' # How the fitted risk moves with mean blood pressure
#' pd <- partial_dependence(fit, ~ meanbp)
#' pd
#'
#' plot(pd)
#'
#' # The same thing from the fit, which is what the `plot()` method is for
#' plot(fit, ~ meanbp)
#'
#' # Two predictors, one of them a factor, which gives a curve per level
#' plot(fit, ~ meanbp + sex)
#'
#' @export
partial_dependence <- function(object, variables, newdata = NULL, grid = 51L,
                               values = NULL, level = 0.95, type = "response",
                               plot = FALSE, ...) {

  arg::arg_is(object, "bartisan_fit")
  arg::arg_whole_number(grid)
  arg::arg_gte(grid, 2)
  arg::arg_number(level)
  arg::arg_between(level, c(0, 1), inclusive = FALSE)
  arg::arg_flag(plot)

  vars <- pd_variables(variables)

  if (is_null(newdata)) {
    newdata <- object[["model"]]

    if (is_null(newdata)) {
      arg::err("this fit kept no data to average over, so {.arg newdata} is
                required")
    }
  }

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  missing <- setdiff(vars, names(newdata))

  if (!is_null(missing)) {
    arg::err("{.arg variables} names {length(missing)} column{?s} the data does
              not have: {.val {missing}}")
  }

  points <- lapply(vars, function(v) pd_grid(newdata[[v]], grid, values[[v]]))
  names(points) <- vars
  combos <- expand.grid(points, KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)

  predictor <- pd_predictor(object, newdata, vars, type, list(...))

  at_grid_point <- function(i) {
    d <- newdata

    for (v in vars) {
      d[[v]] <- pd_assign(newdata[[v]], combos[[v]][[i]])
    }

    draws <- predictor(d)

    if (!is.matrix(draws)) {
      arg::err(c("{.code type = \"{type}\"} does not give one number per
                  observation for this family, so there is nothing to average",
                 i = "Try {.code type = \"mean\"} or
                      {.code type = \"stdlv\"}."))
    }

    # Averaged within each draw, so the interval is on the average prediction.
    s <- post_summary(rowMeans(draws), level = level)
    data.frame(estimate = s[["mean"]], lower = s[["lower"]],
               upper = s[["upper"]])
  }

  # One prediction per grid point, over the whole sample each time, and none of
  # them depends on another, so they go to workers when a plan has any. The
  # streams are drawn before the branch so that both use the same ones: a family
  # whose prediction simulates rather than evaluating in closed form would
  # otherwise give different draws depending on whether a plan was set.
  seeds <- parallel_streams(nrow(combos))

  if (use_future()) {
    rows <- future.apply::future_lapply(seq_len(nrow(combos)), at_grid_point,
                                        future.seed = seeds,
                                        future.packages = "bartisan")
  }
  else {
    restore <- restore_stream()
    on.exit(restore(), add = TRUE)

    rows <- lapply(seq_len(nrow(combos)), function(i) {
      assign(".Random.seed", seeds[[i]], envir = globalenv())
      at_grid_point(i)
    })
  }

  out <- cbind(combos, do.call(rbind, rows)) |>
    unrowname()

  attr(out, "variables") <- vars
  attr(out, "level") <- level
  attr(out, "type") <- type
  attr(out, "n_units") <- nrow(newdata)
  class(out) <- c("bartisan_partial", "data.frame")

  if (!plot) {
    return(out)
  }

  plot(out)
}

pd_variables <- function(variables) {
  arg::arg_or(variables,
              arg::arg_formula(one_sided = TRUE),
              arg::arg_character)
  vars <- {
    if (rlang::is_formula(variables)) get_varnames(variables)
    else variables
  }

  if (is_null(vars) || length(vars) > 2L) {
    arg::err(c("{.arg variables} must name one or two predictors.",
               i = "Three dimensions of dependence is a table rather than a
                    plot; call this twice, or use
                    {.fn marginaleffects::predictions} with a grid of your
                    own."))
  }

  vars
}

# A factor is evaluated at its levels, a numeric predictor on an evenly spaced
# grid across its range, and anything with few distinct values at those values,
# since a grid of 25 points over 3 of them invents 22.
pd_grid <- function(z, grid, given) {
  if (!missing(given) && !is_null(given)) {
    return(given)
  }

  if (is.factor(z)) {
    return(levels(droplevels(z)))
  }

  if (is.character(z) || is.logical(z)) {
    return(sort(unique(stats::na.omit(z))))
  }

  u <- unique(stats::na.omit(z))

  if (length(u) <= grid) {
    return(sort(u))
  }

  seq(min(z, na.rm = TRUE), max(z, na.rm = TRUE), length.out = grid)
}

pd_assign <- function(z, value) {
  if (is.factor(z)) {
    return(factor(rep(as.character(value), length(z)), levels = levels(z)))
  }

  if (is.character(z)) {
    return(rep(as.character(value), length(z)))
  }

  if (is.logical(z)) {
    return(rep(as.logical(value), length(z)))
  }

  rep(as.numeric(value), length(z))
}

#' @rdname partial_dependence
#' @export
print.bartisan_partial <- function(x, digits = 3L, ...) {

  arg::arg_whole_number(digits)

  vars <- attr(x, "variables")

  cli_cat("{.underline Partial dependence}")
  cli::cat_line()
  cli_cat("{cli::qty(length(vars))}Predictor{?s}: {.val {vars}}")
  cli_cat("Averaged over {attr(x, 'n_units')} unit{?s}, on the
           {.val {attr(x, 'type')}} scale")
  cli::cat_line()

  as.data.frame(x) |>
    effect_round(digits) |>
    print(row.names = FALSE)

  cli::cat_line()
  cli_bullets_cat(c(i = "{.field lower} and {.field upper} bound the
                        {100 * attr(x, 'level')}% credible interval on the
                        {.emph average} prediction, not on any one unit's."))

  invisible(x)
}

#' @rdname partial_dependence
#' @export
plot.bartisan_partial <- function(x, ...) {
  vars <- attr(x, "variables")
  d <- as.data.frame(x)
  first <- vars[[1L]]
  discrete <- !is.numeric(d[[first]])

  # The second predictor becomes the grouping, so two numeric predictors give a
  # family of curves rather than a surface nobody can read the uncertainty off.
  group <- if (length(vars) > 1L) vars[[2L]] else NULL

  if (!is_null(group)) {
    d[[group]] <- factor(d[[group]], levels = unique(d[[group]]))
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[first]],
                                       y = .data$estimate))

  if (discrete) {
    pos <- ggplot2::position_dodge(width = 0.3)

    if (is_null(group)) {
      p <- p +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower,
                                            ymax = .data$upper),
                               width = 0, color = "grey60") +
        ggplot2::geom_point()
    }
    else {
      p <- p +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower,
                                            ymax = .data$upper,
                                            color = .data[[group]]),
                               width = 0, position = pos) +
        ggplot2::geom_point(ggplot2::aes(color = .data[[group]]),
                            position = pos)
    }
  }
  else if (is_null(group)) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                        ymax = .data$upper),
                           fill = "steelblue", alpha = 0.2) +
      ggplot2::geom_line()
  }
  else {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                                        fill = .data[[group]]), alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(color = .data[[group]]))
  }

  p +
    ggplot2::labs(x = first, y = "Average prediction") +
    ggplot2::theme_bw()
}

#' @rdname partial_dependence
#' @export
plot.bartisan_fit <- function(x, y, ...) {
  if (missing(y)) {
    arg::err(c("{.fn plot} on a fit draws partial dependence and needs to be
                told which predictors to draw it on.",
               i = "For example {.code plot(fit, ~ age)}.",
               i = "{.fn variable_importance} is where to look for which
                    predictors are worth asking about."))
  }

  partial_dependence(x, variables = y, ..., plot = TRUE)
}

# How a grid point is predicted, chosen once for the whole grid.
#
# The additive predictor is a sum over trees, and a tree that never splits on a
# column the grid varies returns the same value however that column is set. So
# those trees are evaluated once, on the data as it stands, and only the rest are
# evaluated again at each grid point and added to them. On an `rhc` fit most
# predictors are split on by well under a fifth of the trees, so most of the
# forest is walked once rather than `grid` times.
#
# Falls back to `predict()` whenever the split cannot be made safely, which is
# what `pd_tree_mask()` reports by returning `NULL`, and for the prediction types
# that need more than the predictor and the nuisance parameters to evaluate.
# The fallback is the path this replaces, so the two cannot disagree about a fit
# the optimization declines.
pd_predictor <- function(object, newdata, vars, type, dots) {
  ordinary <- function(d) {
    do.call(stats::predict,
            c(list(object, newdata = d, type = type, draws = TRUE), dots))
  }

  if (!type %in% c("link", "response", "stdlv", "mean")) {
    return(ordinary)
  }

  # The fast path reads what it needs out of `...` by name and assembles the
  # predictor itself, so an argument it does not know about would be accepted
  # and then quietly not applied. `predict()` is the thing that knows what to do
  # with one, so an unrecognized name, or a positional argument there is no name
  # to match, is a reason to hand the whole grid back to it. Slower on those
  # calls and never silently wrong.
  handled <- c("offset", "iterations", "weights", "values", "log", "times")
  supplied <- names(dots) %or% character(length(dots))

  if (!is_null(dots) && !all(supplied %in% handled)) {
    return(ordinary)
  }

  uses <- pd_tree_mask(object, vars)

  # Nothing to save when every tree moves, since the base would then be empty
  # and each grid point would do the whole forest anyway, with one wasted pass
  # on top. When no tree moves the base is the whole forest and the grid costs
  # nothing, which is the best case rather than a reason to decline.
  if (is_null(uses) || all(uses)) {
    return(ordinary)
  }

  offset <- dots[["offset"]]
  iterations <- resolve_iterations(dots[["iterations"]],
                                   nrow(object[["sigma_mu"]]))

  # The part that cannot move, including the intercept, the offset and the
  # random effects, which belong to the predictor rather than to any tree. Taken
  # at `newdata` because the columns the grid varies are the ones these trees do
  # not read.
  base <- predict_eta(object, newdata, offset, iterations,
                      tree_mask = !uses, constants = TRUE)

  aux <- {
    if (is_null(object[["aux"]])) NULL
    else object[["aux"]][iterations, , drop = FALSE]
  }

  function(d) {
    eta <- predict_eta(object, d, offset, iterations,
                       tree_mask = uses, constants = FALSE)

    for (h in seq_along(eta)) {
      eta[[h]] <- eta[[h]] + base[[h]]
    }

    # Combined after the sum rather than before it, because a varying
    # coefficient multiplies its forest by the covariate and the covariate is
    # one of the things the grid moves.
    parts <- list(eta = vc_combine(object, eta, d, iterations),
                  aux = aux, iterations = iterations)

    eta_to_type(object, parts, type, draws = TRUE, newdata = d,
                weights = dots[["weights"]], values = dots[["values"]],
                log = dots[["log"]] %or% FALSE, times = dots[["times"]])
  }
}

# Which stored trees can move across the grid, or `NULL` when that cannot be
# decided safely and the caller should predict the ordinary way.
#
# Marking a tree as moving when it does not is only slower, so everything
# uncertain is marked moving; the one case that cannot be handled that way is a
# part of the predictor that is not a tree at all.
pd_tree_mask <- function(object, vars) {
  labels <- object[["term_labels"]]
  assign <- object[["assign"]]

  if (is_null(labels) || is_null(assign) || is_null(object[["forest_flat"]])) {
    return(NULL)
  }

  # A `(1 | g)` term's intercepts are drawn parameters rather than trees, so no
  # choice of trees holds them fixed while the grid moves `g`. Refused rather
  # than special-cased, since a grid over a grouping factor is rare.
  if (any(vars %in% names(object[["random"]]))) {
    return(NULL)
  }

  # A term moves if it mentions any of the plotted variables. A label that
  # cannot be parsed is treated as moving.
  moves <- vapply(labels, function(l) {
    named <- rlang::try_fetch(get_varnames(str2lang(l)),
                              error = function(cnd) NULL)
    is_null(named) || any(vars %in% named)
  }, logical(1L))

  # A `bcf()` propensity score is a function of the covariates, so it moves with
  # the grid. Its own columns are named in the fit rather than in the formula,
  # and marking them all is the cheap side of the trade.
  score <- object[["bcf"]][["propensity"]]

  if (!is_null(score)) {
    moves[labels %in% colnames(score)] <- TRUE
  }

  groups <- sort(unique(assign))
  moving <- groups[moves[groups]]

  if (is_null(moving)) {
    return(NULL)
  }

  # Zero-based, as the engine indexes them: the design columns for the numeric
  # rules, and the level-code columns for the rules carrying a level mask.
  num_cols <- which(assign %in% moving) - 1L

  info <- object[["level_codes"]]
  cat_cols <- integer()

  if (!is_null(info) && any(info[["n_levels"]] > 0L)) {
    pos <- match(moving, groups)
    pos <- pos[info[["n_levels"]][pos] > 0L]
    cat_cols <- info[["cat_col"]][pos]
  }

  .bartisan_tree_uses(forest_flat = object[["forest_flat"]],
                      tree_start = object[["tree_start"]],
                      num_forest = object[["num_forest"]],
                      num_trees = object[["num_trees"]],
                      num_draws = nrow(object[["sigma_mu"]]),
                      num_cols = as.integer(num_cols),
                      cat_cols = as.integer(cat_cols))
}
