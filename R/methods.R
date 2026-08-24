#' Summarize a generalized BART model
#'
#' @description
#' `print()` reports what was fit and how long the chain is. `summary()` adds
#' posterior summaries of the nuisance parameters and of how often each
#' predictor was used in a splitting rule, which is the model's variable
#' selection output.
#'
#' @param x,object a fitted model from [genbart()].
#' @param level width of the reported posterior intervals.
#' @param digits number of significant digits to print.
#' @param ... ignored, present for compatibility with the generics.
#'
#' @returns
#' `print()` returns its argument invisibly. `summary()` returns a list of class
#' `summary.genbart`, with a `usage` element giving the posterior summary of the
#' splitting counts for each predictor group, and an `aux` element for the
#' nuisance parameters when the family has any.
#'
#' @seealso [genbart()]
#'
#' @examples
#' set.seed(1)
#'
#' n <- 150
#' d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
#' d$y <- 2 * d$x1 + rnorm(n, sd = 0.3)
#'
#' fit <- genbart(y ~ x1 + x2 + x3, data = d,
#'                control = genbart_control(num_trees = 10, num_burn = 50,
#'                                          num_save = 50, verbose = FALSE))
#' fit
#' summary(fit)
#'
#' @export
print.genbart <- function(x, digits = 3, ...) {

  print_header(x[["call"]])

  cli::cat_line(cli::format_inline("Family: {family_label(x[['family']])}"))
  cli::cat_line(cli::format_inline("Observations: {x$n}"))
  rules <- if (x[["soft"]]) "soft" else "hard"
  cli::cat_line(cli::format_inline(
    "Forests: {x$num_forest} of {x$num_trees} trees, {rules} decision rules"))
  cli::cat_line(cli::format_inline(
    "Draws: {nrow(x$sigma_mu)} kept after {x$control$num_burn} warmup"))

  if (!is_null(x[["random"]])) {
    cli::cat_line(cli::format_inline(
      "Random intercepts: {ranef_label(x)}"))
  }

  if (!is_null(x[["aux"]])) {
    cli::cat_line()
    means <- signif(colMeans(x[["aux"]]), digits)
    cli::cat_line(cli::format_inline(
      "Posterior means: {paste0(names(means), ' = ', means, collapse = ', ')}"))
  }

  invisible(x)
}

# "g (40 levels), site (7 levels)", which is the part of the model that the
# formula says and the family label does not.
ranef_label <- function(x) {
  terms <- x[["random"]]
  paste0(vapply(terms, function(z) {
    paste0(z[["label"]], " (", z[["num_levels"]], " levels)")
  }, character(1L)), collapse = ", ")
}

# Posterior summaries of each random-effect scale, one row per grouping factor
# per additive predictor. A group effect's scale is the quantity a reader wants
# from it -- how much of the variation is between groups -- so it is reported
# next to the leaf scale rather than left in the fit.
ranef_summary <- function(object, level) {
  if (is_null(object[["tau"]])) {
    return(NULL)
  }

  rows <- list()

  for (h in seq_along(object[["tau"]])) {
    tau <- object[["tau"]][[h]]

    for (r in seq_len(ncol(tau))) {
      label <- {
        if (length(object[["tau"]]) > 1L) {
          paste0(colnames(tau)[r], " [", names(object[["tau"]])[h], "]")
        }
        else colnames(tau)[r]
      }

      rows[[label]] <- post_summary(tau[, r], level = level)
    }
  }

  do.call(rbind, rows)
}

# A family built by custom_family() reports the name the caller gave it, since
# "custom" with the "identity" link says nothing. A link the engine does not
# carry natively is flagged, because it is applied on the R side.
family_label <- function(family) {
  if (identical(family[["family"]], "custom")) {
    return(cli::format_inline("{.val {family$name}} (supplied from R)"))
  }

  supplied <- {
    if (is_null(family[["custom_link"]])) ""
    else " (supplied from R)"
  }

  cli::format_inline(
    "{.val {family$family}} with the {.val {family$link}} link{supplied}")
}

# Print methods write to stdout, so they use cli's cat_* functions rather than
# cli_text(), which emits a message on stderr and would leave the output
# invisible to capture.output() and to knitr.
print_header <- function(call) {
  cli::cat_line(cli::format_inline("{.strong Generalized BART}"))
  cli::cat_line()
  cli::cat_line(cli::format_inline(
    "Call: {.code {deparse(call, width.cutoff = 500L)[1L]}}"))
  cli::cat_line()
}

#' @rdname print.genbart
#' @export
summary.genbart <- function(object, level = 0.95, ...) {

  arg::arg_number(level)
  arg::arg_between(level, c(0, 1), inclusive = FALSE)

  usage <- lapply(object[["counts"]], function(counts) {
    out <- t(apply(counts, 2L, post_summary, level = level))
    # The proportion of draws in which a group was used at all is the more
    # readable variable-selection summary than the raw count.
    out <- cbind(out, prop_used = colMeans(counts > 0))
    out[order(out[, "prop_used"], decreasing = TRUE), , drop = FALSE]
  })

  aux <- {
    if (is_null(object[["aux"]])) NULL
    else t(apply(object[["aux"]], 2L, post_summary, level = level))
  }

  out <- list(call = object[["call"]],
              family = object[["family"]],
              n = object[["n"]],
              num_forest = object[["num_forest"]],
              num_trees = object[["num_trees"]],
              soft = object[["soft"]],
              num_save = nrow(object[["sigma_mu"]]),
              level = level,
              usage = usage,
              aux = aux,
              random = object[["random"]],
              tau = ranef_summary(object, level),
              sigma_mu = t(apply(object[["sigma_mu"]], 2L, post_summary,
                                 level = level)),
              loglik = post_summary(object[["loglik"]], level = level))

  class(out) <- "summary.genbart"

  out
}

#' @rdname print.genbart
#' @export
print.summary.genbart <- function(x, digits = 3, ...) {

  print_header(x[["call"]])

  cli::cat_line(cli::format_inline("Family: {family_label(x[['family']])}"))
  cli::cat_line(cli::format_inline("Observations: {x$n}"))
  rules <- if (x[["soft"]]) "soft" else "hard"
  cli::cat_line(cli::format_inline(
    "Forests: {x$num_forest} of {x$num_trees} trees, {rules} decision rules"))
  cli::cat_line(cli::format_inline("Draws: {x$num_save}"))

  if (!is_null(x[["random"]])) {
    cli::cat_line(cli::format_inline(
      "Random intercepts: {ranef_label(x)}"))
  }

  if (!is_null(x[["aux"]])) {
    cli::cat_line()
    cli::cat_line(cli::format_inline("{.strong Nuisance parameters}"))
    print(round(x[["aux"]], digits))
  }

  if (!is_null(x[["tau"]])) {
    cli::cat_line()
    cli::cat_line(cli::format_inline("{.strong Random-effect scales}"))
    cli::cat_line(cli::format_inline(
      "{.emph Standard deviation of the group intercepts.}"))
    print(round(x[["tau"]], digits))
  }

  cli::cat_line()
  cli::cat_line(cli::format_inline("{.strong Predictor usage}"))
  cli::cat_line(cli::format_inline(
    "{.emph Splitting rules per draw, and how often used at all.}"))

  for (h in seq_along(x[["usage"]])) {
    if (length(x[["usage"]]) > 1L) {
      cli::cat_line()
      cli::cat_line(cli::format_inline(
        "Predictor {.val {names(x$usage)[h]}}:"))
    }
    print(round(x[["usage"]][[h]], digits))
  }

  invisible(x)
}
