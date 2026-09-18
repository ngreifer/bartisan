#' Summarize a generalized BART model
#'
#' @description
#' `print()` reports what was fit and how long the chain is. `summary()` adds
#' posterior summaries of the nuisance parameters, of the random-effect scales,
#' and of how often each predictor was used in a splitting rule, which is the
#' model's variable-selection output.
#'
#' @param x,object a `<bartisan_fit>` object; the output of a call to
#'   [bartisan()]. For `print.summary.bartisan_fit()`, `x` is instead the
#'   `<summary.bartisan_fit>` object `summary()` returned.
#' @param level `numeric`; the width of the posterior intervals `summary()`
#'   reports. Default is .95 for 95% intervals.
#' @param digits `numeric`; the number of significant digits to print.
#'   Default is 3.
#' @param ... not used.
#'
#' @returns
#' `print()` returns its argument invisibly. `summary()` returns a
#' `<summary.bartisan_fit>` object, a list of the posterior summaries its own
#' `print()` method displays: the nuisance parameters when the family has any,
#' the scale of each random-effect term, and the splitting counts of each
#' predictor group.
#'
#' @seealso
#' [bartisan()]; [variable_importance()] for the splitting counts as a data
#' frame rather than printed; [diagnose()] for whether the chains the summaries
#' are computed from have converged
#'
#' @examples
#' data("rhc")
#' set.seed(123)
#'
#' # Whether a patient died, with catheterization among the predictors
#' fit <- bartisan(death ~ . - days, data = rhc, num_trees = 10,
#'                 num_burn = 50, num_draws = 50, chains = 2, verbose = FALSE)
#'
#' # What was fit, and how many draws it rests on
#' fit
#'
#' # The splitting counts, which say which predictors the forest reaches for
#' summary(fit)
#'
#' # A family with a nuisance parameter reports its posterior too, here the
#' # residual standard deviation of the log survival time
#' fit2 <- bartisan(log(days) ~ . - death, data = rhc, family = gaussian(),
#'                  num_trees = 10, num_burn = 50, num_draws = 50,
#'                  verbose = FALSE)
#'
#' summary(fit2, level = .8)
#'
#' @export
print.bartisan_fit <- function(x, digits = 3L, ...) {

  print_header(x[["call"]])

  cli_cat("Family: {family_label(x[['family']])}")
  cli_cat("Observations: {x$n}")
  rules <- if (x[["soft"]]) "soft" else "hard"
  cli_cat("Structure: {forest_label(x)}, {rules} decision rules")

  if (x$control$chains == 1) {
    cli_cat("Draws: {nrow(x$sigma_mu)} kept after {x$control$num_burn} warmup")
  }
  else {
    cli_cat("Draws: {nrow(x$sigma_mu)} kept across {x$control$chains} chains after {x$control$num_burn} warmup")
  }

  if (!is_null(x[["random"]])) {
    cli_cat("Random intercepts: {ranef_label(x)}")
  }

  if (!is_null(x[["aux"]])) {
    arg::arg_whole_number(digits)

    cli::cat_line()
    means <- signif(colMeans(x[["aux"]]), digits)
    labels <- sprintf("%s = %s", names(means), means) |>
      toString()
    cli_cat("Posterior means: {labels}")
  }

  invisible(x)
}

# "1 forest of 50 trees", or "2 forests of 50 and 10 trees" when the tree count
# differs by additive predictor.
forest_label <- function(x) {
  trees <- x[["num_trees"]]

  .lab <- {
    if (length(unique(trees)) == 1L)
      "{x$num_forest} forest{?s} of {trees[1L]} tree{?s}"
    else
      "{x$num_forest} forests of {.and {trees}} trees"
  }

  cli::format_inline(.lab)
}

# "g (40 levels), site (7 levels)", which is the part of the model that the
# formula says and the family label does not.
ranef_label <- function(x) {
  vapply(x[["random"]], function(z) {
    sprintf("%s (%s levels)", z[["label"]], z[["num_levels"]])
  }, character(1L)) |>
    toString()
}

# Posterior summaries of each random-effect scale, one row per grouping factor
# per additive predictor. A group effect's scale is the quantity a reader wants
# from it -- how much of the variation is between groups -- so it is reported
# next to the leaf scale rather than left in the fit.
ranef_summary <- function(object, level) {
  if (is_null(object[["tau"]])) {
    return(NULL)
  }

  # One row per random-effect standard deviation, which is one per column of
  # each `tau`. The labels are collected alongside rather than used as the
  # index, since `rbind()` reads them off the names at the end.
  n <- sum(vapply(object[["tau"]], ncol, integer(1L)))
  rows <- vector("list", n)
  labels <- character(n)
  at <- 0L

  for (h in seq_along(object[["tau"]])) {
    tau <- object[["tau"]][[h]]

    for (r in seq_len(ncol(tau))) {
      at <- at + 1L

      labels[at] <- {
        if (length(object[["tau"]]) > 1L) {
          sprintf("%s [%s]", colnames(tau)[r], names(object[["tau"]])[h])
        }
        else colnames(tau)[r]
      }

      rows[[at]] <- post_summary(tau[, r], level = level)
    }
  }

  rows |>
    setNames(labels) |>
    do_rbind()
}

# A family built by custom_family() reports the name the caller gave it, since
# "custom" with the "identity" link says nothing. A link the engine does not
# carry natively is flagged, because it is applied on the R side.
family_label <- function(family) {
  if (identical(family[["family"]], "custom")) {
    return(cli::format_inline("{.val {family$name}} (supplied from R)"))
  }

  # The engine calls it "mnp"; the caller asked for a multinomial with a probit
  # link, and that is what the fit should say it is.
  if (identical(family[["family"]], "mnp")) {
    return(cli::format_inline("{.val multinomial} with the {.val probit} link"))
  }

  supplied <- {
    if (is_null(family[["custom_link"]])) ""
    else " (supplied from R)"
  }

  cli::format_inline("{.val {family$family}} with the {.val {family$link}} link{supplied}")
}

# Print methods write to stdout, so they use cli's cat_* functions rather than
# cli_text(), which emits a message on stderr and would leave the output
# invisible to capture.output() and to knitr.
print_header <- function(call) {
  cli_cat("{.underline Generalized BART}")
  cli::cat_line()
  cli::cat_line("Call:\n", paste(deparse(call), collapse = "\n"))
  cli::cat_line()
}

#' @rdname print.bartisan_fit
#' @export
summary.bartisan_fit <- function(object, level = 0.95, ...) {

  arg::arg_number(level)
  arg::arg_between(level, c(0, 1), inclusive = FALSE)

  usage <- lapply(object[["counts"]], function(counts) {
    out <- t(apply(counts, 2L, post_summary, level = level))
    # The proportion of draws in which a group was used at all is the more
    # readable variable-selection summary than the raw count.
    out <- cbind(out, prop_used = colMeans(counts > 0))
    out[order(out[, "prop_used"], decreasing = TRUE), , drop = FALSE]
  })

  aux <- if (!is_null(object[["aux"]])) t(apply(object[["aux"]], 2L, post_summary, level = level))

  out <- list(call = object[["call"]],
              family = object[["family"]],
              n = object[["n"]],
              # Recorded so that `print()` can point a reader at the function
              # that reports the effect, which this summary deliberately does
              # not: a named treatment is the whole of what makes that relevant.
              treatment = object[["bcf"]][["treatment"]],
              num_forest = object[["num_forest"]],
              num_trees = object[["num_trees"]],
              soft = object[["soft"]],
              num_draws = nrow(object[["sigma_mu"]]),
              level = level,
              usage = usage,
              aux = aux,
              random = object[["random"]],
              tau = ranef_summary(object, level),
              sigma_mu = t(apply(object[["sigma_mu"]], 2L, post_summary,
                                 level = level)),
              loglik = post_summary(object[["loglik"]], level = level))

  class(out) <- "summary.bartisan_fit"

  out
}

#' @rdname print.bartisan_fit
#' @export
print.summary.bartisan_fit <- function(x, digits = 3, ...) {

  print_header(x[["call"]])

  cli_cat("Family: {family_label(x[['family']])}")
  cli_cat("Observations: {x$n}")
  rules <- if (x[["soft"]]) "soft" else "hard"
  cli_cat("Structure: {forest_label(x)}, {rules} decision rules")
  cli_cat("Draws: {x$num_draws}")

  if (!is_null(x[["random"]])) {
    cli_cat("Random intercepts: {ranef_label(x)}")
  }

  if (!is_null(x[["aux"]])) {
    cli::cat_line()
    cli_cat("{.underline Nuisance parameters}")

    # A baseline hazard can have one entry per event time, which is too many to
    # read. Printing the ends and saying how many were left out keeps the block
    # legible without hiding that they are all there in `fit$aux`.
    shown <- 12L

    if (nrow(x[["aux"]]) > shown) {
      keep <- c(seq_len(shown %/% 2L),
                seq(nrow(x[["aux"]]) - shown %/% 2L + 1L, nrow(x[["aux"]])))
      print(round(x[["aux"]][keep, , drop = FALSE], digits))
      cli_cat("{.emph {nrow(x[['aux']]) - length(keep)} more, omitted;
               all of them are in {.code fit$aux}.}")
    }
    else {
      print(round(x[["aux"]], digits))
    }
  }

  if (!is_null(x[["tau"]])) {
    cli::cat_line()
    cli_cat("{.underline Random-effect scales}")
    cli_cat("{.emph Standard deviation of the group intercepts.}")
    print(round(x[["tau"]], digits))
  }

  cli::cat_line()
  cli_cat("{.underline Predictor usage}")
  cli_cat("{.emph Splitting rules per draw, and how often used at all.}")

  for (h in seq_along(x[["usage"]])) {
    if (length(x[["usage"]]) > 1L) {
      cli::cat_line()
      cli_cat("Predictor {.val {names(x$usage)[h]}}:")
    }
    print(round(x[["usage"]][[h]], digits))
  }

  # A fit with a named treatment summarizes the same way as any other, since the
  # forests are the same object; the effect is a different question and
  # `estimate_effect()` is where it is asked.
  if (!is_null(x[["treatment"]])) {
    cli::cat_line()
    cli_bullets_cat(c(i = "This fit has a treatment, {.val {x$treatment}}.
                          {.fn estimate_effect} reports its effect, with the
                          average potential outcomes beside it."))
  }

  invisible(x)
}

#' Varying coefficients
#'
#' The coefficient functions of a model fitted with [vc()] terms, evaluated at
#' each observation. A forest has no coefficient vector, so a fit with no `vc()`
#' term has no coefficient to report and this errors rather than returning
#' anything.
#'
#' @param object a `<bartisan_fit>` object; the output of a call to [bartisan()],
#'   fitted with at least one [vc()] term.
#' @param newdata optional; a data frame at which to evaluate the coefficients.
#'   Default is `NULL` to use the data the model was fitted to.
#' @param draws `logical`; whether to return every posterior draw of each
#'   coefficient rather than its posterior mean at each observation. Default is
#'   `FALSE` to return the posterior means.
#' @param ... not used.
#'
#' @returns
#' With `draws = FALSE`, a matrix with one row per observation and one column per
#' coefficient. With `draws = TRUE`, a named list of draws-by-observations
#' matrices, one per coefficient.
#'
#' @details
#' The coefficients are the varying ones alone. The control function is the
#' surface at the value each covariate was centered on, which is a prediction
#' rather than a coefficient, and `predict(object)` reports it.
#'
#' For a factor the coefficients are recentered to sum to zero across its levels,
#' which is what makes them the deviations they are reported as. The symmetric
#' coding carries one spare function-valued dimension, so this is exact rather
#' than an approximation, and it is the reason a factor's reference level is a
#' choice made here rather than at fitting time.
#'
#' @seealso
#' [vc()] for declaring a varying coefficient; [variable_importance()] for which
#' predictors a forest uses at all
#'
#' @examples
#' data("rhc")
#' set.seed(123)
#'
#' # The effect of catheterization is allowed to vary with the other
#' # predictors, so its coefficient is a function rather than a number
#' fit <- bartisan(death ~ age + aps + surv2m + vc(rhc), data = rhc,
#'                 num_trees = 10, num_burn = 50, num_draws = 50,
#'                 verbose = FALSE)
#'
#' # One coefficient per patient, on the link scale
#' head(coef(fit))
#'
#' # How much it varies across patients
#' quantile(coef(fit)[, "rhc"])
#'
#' @exportS3Method stats::coef
coef.bartisan_fit <- function(object, newdata = NULL, draws = FALSE, ...) {
  vc <- object[["vc"]]

  if ((vc[["slopes"]] %or% 0L) == 0L) {
    arg::err(c("This model has no varying coefficients, and a forest has no
                coefficient vector.",
               i = "Use {.fn variable_importance} for which predictors the
                  forest uses, or
                  {.code marginaleffects::avg_comparisons()} for how much one
                  moves the outcome."))
  }

  arg::arg_flag(draws)

  eta <- {
    if (is_null(newdata)) object[["eta"]]
    else predict_eta(object, newdata, offset = NULL, iterations = NULL)
  }

  # The control functions are dropped: one is a prediction, not a coefficient.
  # With several additive predictors there is one per parameter, so which
  # forests to keep comes from the map rather than from a position.
  keep <- which(vc[["column"]][seq_along(eta)] > 0L)
  slopes <- eta[keep]
  names(slopes) <- names(object[["eta"]])[keep]

  slopes <- vc_recenter(slopes, vc, object)

  if (draws) {
    return(slopes)
  }

  vapply(slopes, colMeans, numeric(ncol(slopes[[1L]]))) |>
    matrix(ncol = length(slopes),
           dimnames = list(NULL, names(slopes)))
}

#' Group intercepts from a random-effect term
#'
#' Extracts the intercept each level of a grouping factor was given by a
#' `(1 | group)` term in the formula, as a posterior mean or as every draw.
#'
#' @param object a `<bartisan_fit>` object; the output of a call to [bartisan()],
#'   fitted with at least one `(1 | group)` term.
#' @param draws `logical`; whether to return every posterior draw of each
#'   intercept rather than its posterior mean. Default is `FALSE`.
#' @param ... not used.
#'
#' @returns
#' With `draws = FALSE`, a named list with one entry per grouping factor, each a
#' data frame of one row per level of that factor and one column per additive
#' predictor that has a group intercept, with the levels as row names. This is
#' the shape \pkgfun{lme4}{ranef} returns, so anything that reads that reads
#' this.
#'
#' With `draws = TRUE`, a named list with one entry per grouping factor, each
#' itself a named list of draws-by-levels matrices, one per additive predictor.
#'
#' @details
#' The column is named `(Intercept)`, as in \pkg{lme4}, when the family has a
#' single additive predictor. A family with more than one gets a group intercept
#' on each, independent of the others, and the columns are named for the
#' predictors instead; `vignette("families")` lists them per family.
#'
#' Only the intercepts are returned. The standard deviation each grouping factor
#' was drawn under is in `object$tau`, one column per factor and one matrix per
#' additive predictor, and [`prior_summary()`][bartisan-interop] reports the prior
#' it was drawn from.
#'
#' A posterior mean is the wrong summary for a level with few observations, which
#' is the case a group intercept exists for. `draws = TRUE` is what gives an
#' interval, and a level whose interval covers zero is one the data had little to
#' say about.
#'
#' The intercepts are shrunk towards zero by their prior and are deviations from
#' the additive predictor, so they come out approximately centered without being
#' constrained to sum to zero exactly. Nothing is lost by that: the level of the
#' fitted function is the additive predictor's, and a shift common to every
#' intercept is one the forest did not take.
#'
#' @seealso
#' [bartisan()] for the `(1 | group)` syntax and what it fits;
#' [`prior_summary()`][bartisan-interop] for the prior on these
#'
#' @examplesIf rlang::is_installed("nlme")
#' set.seed(123)
#' d <- data.frame(x = runif(200),
#'                 site = factor(sample(letters[1:5], 200, TRUE)))
#' d$y <- rnorm(200, d$x + as.numeric(d$site) / 3)
#'
#' fit <- bartisan(y ~ x + (1 | site), data = d, num_trees = 10,
#'                 num_burn = 50, num_draws = 50, verbose = FALSE)
#'
#' # One intercept per site, as posterior means. The generic is \pkg{nlme}'s,
#' # which \pkg{lme4} re-exports, so either qualification reaches this.
#' nlme::ranef(fit)
#'
#' # With the draws, so the intercepts come with intervals
#' apply(nlme::ranef(fit, draws = TRUE)$site[["(Intercept)"]], 2L, quantile,
#'       c(.025, .975))
#'
#' @exportS3Method nlme::ranef
ranef.bartisan_fit <- function(object, draws = FALSE, ...) {

  arg::arg_is(object, "bartisan_fit")
  arg::arg_flag(draws)

  random <- object[["random"]]

  if (is_null(random)) {
    arg::err(c("this model has no group intercepts, so there is nothing to
                extract",
               i = "A {.code (1 | group)} term in the formula is what adds
                    them; see {.fn bartisan}."))
  }

  stored <- object[["ranef"]]

  # One additive predictor means one intercept with nothing to distinguish it,
  # so the column takes the name lme4 gives a random intercept. Several means
  # one per predictor, and then the predictor's name is the informative one.
  labels <- if (length(stored) == 1L) "(Intercept)" else names(stored)

  lapply(random, function(term) {
    # The stored columns are `label:level`, so they are rebuilt from the term
    # rather than parsed out of the names: a level whose own value contains a
    # colon would otherwise split in the wrong place.
    want <- paste0(term[["label"]], ":", term[["levels"]])
    columns <- lapply(stored, function(m) m[, want, drop = FALSE])

    if (draws) {
      columns <- lapply(columns, setColnames, term[["levels"]])

      return(setNames(columns, labels))
    }

    means <- lapply(columns, colMeans)

    data.frame(means, row.names = term[["levels"]], check.names = FALSE) |>
      setNames(labels)
  }) |>
    setNames(pluck(random, "label", character(1L)))
}
