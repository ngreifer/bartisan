#' Methods for Bayesian causal forest fits
#'
#' A fit from [bcf()] is a `<bartisan_fit>` with a treatment named, so it takes
#' every method that fit takes and adds the ones that need a treatment.
#' `summary()` reports the effect rather than the forests, and `plot()` draws the
#' conditional effects.
#'
#' @param object,x a `<bcf_fit>` object; the output of a call to [bcf()].
#' @param level `numeric`; the level of the credible interval. Default is `.95`.
#' @param estimand,comparison,interval,focal passed to [estimate_effect()],
#'   which is what computes the numbers.
#' @param digits `integer`; the number of significant digits to print.
#' @param ... for `summary()`, further arguments passed to
#'   [estimate_effect()]; otherwise ignored.
#'
#' @returns
#' `summary()` returns a `<summary.bcf_fit>` object carrying the effect, the
#' marginal mean of each potential outcome, and the quartiles of the conditional
#' effects. `plot()` returns a \pkg{ggplot2} object. Both `print()` methods
#' return their input invisibly.
#'
#' @details
#' `summary()` is [estimate_effect()] with defaults and a compact printout, so
#' anything it does not report is available by calling that function directly:
#' another estimand, a ratio rather than a difference, subgroup effects, or the
#' posterior draws themselves.
#'
#' The potential outcomes are printed beside the effect because a difference
#' means different things against different baselines, which is the point
#' `vignette("causal")` makes about earnings.
#'
#' @seealso [estimate_effect()] for the estimands and for what `summary()` does
#'   not cover; [bcf()] for the model
#'
#' @examplesIf rlang::is_installed("ggplot2")
#' data("rhc")
#' set.seed(123)
#'
#' fit <- bcf(death ~ age + sex + meanbp + aps, treatment = ~ rhc,
#'            data = rhc, num_trees = 10, num_burn = 50, num_draws = 50,
#'            verbose = FALSE)
#'
#' fit
#'
#' summary(fit)
#'
#' # The effect on the treated instead, still through `summary()`
#' summary(fit, estimand = "ATT")
#'
#' # The conditional effects, ordered, with the average behind them
#' plot(fit)
#'
#' @export
summary.bcf_fit <- function(object, level = 0.95, estimand = "ATE",
                            comparison = "difference", interval = "eti",
                            focal = NULL, ...) {

  effect <- estimate_effect(object, estimand = estimand,
                            comparison = comparison, level = level,
                            interval = interval, focal = focal, ...)

  # The conditional effects come from the same potential outcomes, so this is a
  # second pass over draws already taken rather than a second fit.
  cate <- estimate_effect(object, estimand = "CATE", comparison = comparison,
                          level = level, interval = interval, focal = focal,
                          ...)

  quartiles <- do.call(rbind, lapply(split(cate, cate[["contrast"]]),
                                     function(z) {
    q <- stats::quantile(z[["estimate"]], c(0, 0.25, 0.5, 0.75, 1),
                         names = FALSE)
    data.frame(contrast = z[["contrast"]][1L], min = q[1L], q25 = q[2L],
               median = q[3L], q75 = q[4L], max = q[5L],
               stringsAsFactors = FALSE)
  })) |>
    unrowname()

  out <- list(call = object[["call"]],
              family = object[["family"]],
              n = object[["n"]],
              treatment = object[["bcf"]][["treatment"]],
              moderators = object[["bcf"]][["moderators"]],
              effect = effect,
              potential_outcomes = attr(effect, "potential_outcomes"),
              cate = quartiles,
              level = level)

  class(out) <- "summary.bcf_fit"

  out
}

#' @rdname summary.bcf_fit
#' @export
print.summary.bcf_fit <- function(x, digits = 3L, ...) {

  arg::arg_whole_number(digits)

  print_header(x[["call"]])

  cli_cat("Family: {family_label(x[['family']])}")
  cli_cat("Observations: {x$n}")
  cli::cat_line()

  print(x[["effect"]], digits = digits)

  cli::cat_line()
  cli_cat("{.strong Potential outcomes}")
  cli::cat_line()
  print(effect_round(as.data.frame(x[["potential_outcomes"]]), digits),
        row.names = FALSE)

  cli::cat_line()
  cli_cat("{.strong Conditional effects}")
  cli::cat_line()
  print(effect_round_all(x[["cate"]], digits), row.names = FALSE)
  cli::cat_line()
  cli::cli_bullets(c(i = "Quartiles of the per-unit posterior means, which say
                          how much the effect varies rather than how well any
                          one unit is estimated."))

  invisible(x)
}

effect_round_all <- function(show, digits) {
  for (nm in setdiff(names(show), "contrast")) {
    show[[nm]] <- signif(show[[nm]], digits)
  }

  show
}

#' @rdname summary.bcf_fit
#' @export
print.bcf_fit <- function(x, digits = 3L, ...) {

  NextMethod()

  cli::cat_line()
  cli_cat("Treatment: {.val {x$bcf$treatment}}")

  if (!is_null(x[["bcf"]][["moderators"]])) {
    cli_cat("Effect moderators: {.val {x$bcf$moderators}}")
  }

  cli::cli_bullets(c(i = "{.fn summary} reports the average treatment effect;
                          {.fn estimate_effect} gives another estimand, a ratio
                          rather than a difference, or subgroup effects."))

  invisible(x)
}

#' @rdname summary.bcf_fit
#' @export
plot.bcf_fit <- function(x, level = 0.95, comparison = "difference", ...) {
  plot(estimate_effect(x, estimand = "CATE", comparison = comparison,
                       level = level, ...))
}
