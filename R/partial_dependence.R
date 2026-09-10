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
#'   Default is 25. A factor is evaluated at each of its levels whatever this is.
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
#' @param y for `plot.bartisan_fit()`, the predictors to plot, as `variables`
#'   above.
#' @param ... for `plot.bartisan_fit()`, further arguments passed to
#'   `partial_dependence()`; otherwise ignored.
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
#' @examplesIf rlang::is_installed("ggplot2")
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
partial_dependence <- function(object, variables, newdata = NULL, grid = 25L,
                               values = NULL, level = 0.95, type = "response",
                               plot = FALSE) {

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

  if (length(missing) > 0L) {
    arg::err("{.arg variables} names {length(missing)} column{?s} the data does
              not have: {.val {missing}}")
  }

  points <- lapply(vars, function(v) pd_grid(newdata[[v]], grid, values[[v]]))
  names(points) <- vars
  combos <- expand.grid(points, KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)

  rows <- lapply(seq_len(nrow(combos)), function(i) {
    d <- newdata

    for (v in vars) {
      d[[v]] <- pd_assign(newdata[[v]], combos[[v]][[i]])
    }

    draws <- stats::predict(object, newdata = d, type = type, draws = TRUE)

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
  })

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
  vars <- {
    if (rlang::is_formula(variables)) {
      arg::arg_formula(variables, one_sided = TRUE)
      all.vars(variables)
    }
    else {
      arg::arg_character(variables)
      variables
    }
  }

  if (length(vars) < 1L || length(vars) > 2L) {
    arg::err(c("{.arg variables} must name one or two predictors",
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
  if (!is_null(given)) {
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

  cli_cat("{.strong Partial dependence}")
  cli::cat_line()
  cli_cat("{cli::qty(length(vars))}Predictor{?s}: {.val {vars}}")
  cli_cat("Averaged over {attr(x, 'n_units')} unit{?s}, on the
           {.val {attr(x, 'type')}} scale")
  cli::cat_line()

  print(effect_round(as.data.frame(x), digits), row.names = FALSE)
  cli::cat_line()
  cli::cli_bullets(c(i = "{.field lower} and {.field upper} bound the
                          {100 * attr(x, 'level')}% credible interval on the
                          {.emph average} prediction, not on any one unit's."))

  invisible(x)
}

#' @rdname partial_dependence
#' @export
plot.bartisan_partial <- function(x, ...) {
  require_ggplot2("partial dependence")

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
                told which predictors to draw it on",
               i = "For example {.code plot(fit, ~ age)}.",
               i = "{.fn variable_importance} is where to look for which
                    predictors are worth asking about."))
  }

  partial_dependence(x, variables = y, ..., plot = TRUE)
}
