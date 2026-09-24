#' Methods for Bayesian causal forest fits
#'
#' A fit from [bcf()] is a `<bartisan_fit>` with a treatment named, so it takes
#' every method that fit takes. `print()` adds the treatment and says where the
#' effect comes from, and `plot()` draws the conditional effects.
#'
#' @param x a `<bcf_fit>` object; the output of a call to [bcf()].
#' @param level `numeric`; the level of the credible interval. Default is `.95`.
#' @param comparison passed to [estimate_effect()], which is what computes the
#'   numbers.
#' @param marginal `logical`; whether to draw the marginal effect beside the
#'   conditional ones. Default is `TRUE`.
#' @param digits `integer`; the number of significant digits to print.
#' @param ... for `plot()`, further arguments passed to [estimate_effect()];
#'   otherwise ignored.
#'
#' @returns
#' `plot()` returns a \pkg{ggplot2} object. `print()` returns its input
#' invisibly.
#'
#' @details
#' The effect itself is [estimate_effect()]'s to report, and `summary()` on a
#' `<bcf_fit>` is the same summary of the forests it is on any other fit, with a
#' line at the end saying so. That way the same call prints the same thing
#' whether the model came from [bcf()] or from [bartisan()] with a [vc()] term.
#'
#' @seealso [estimate_effect()] for the effect and its estimands; [bcf()] for the
#'   model
#'
#' @examples
#' data("rhc")
#' set.seed(123)
#'
#' fit <- bcf(death ~ age + sex + meanbp + aps, treat = ~ rhc,
#'            data = rhc, num_trees = 10, num_burn = 50, num_draws = 50,
#'            verbose = FALSE)
#'
#' fit
#'
#' # The effect, with the potential outcomes it is a difference of
#' estimate_effect(fit)
#'
#' # The conditional effects, ordered, with the marginal effect beside them
#' plot(fit)
#'
#' @rdname print.bcf_fit
#' @export
print.bcf_fit <- function(x, digits = 3L, ...) {

  NextMethod()

  cli::cat_line()
  cli_cat("Treatment: {.val {x$bcf$treatment}}")

  if (!is_null(x[["bcf"]][["moderators"]])) {
    cli_cat("Effect moderators: {.val {x$bcf$moderators}}")
  }

  cli_bullets_cat(c(i = "{.fn estimate_effect} reports the treatment effect,
                        with the average potential outcomes beside it;
                        {.fn plot} draws the conditional ones."))

  invisible(x)
}

#' @rdname print.bcf_fit
#' @export
plot.bcf_fit <- function(x, level = 0.95, comparison = "difference",
                         marginal = TRUE, ...) {
  plot(estimate_effect(x, estimand = "CATE", comparison = comparison,
                       level = level, ...),
       marginal = marginal)
}
