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

  rows <- list()

  for (h in seq_along(object[["tau"]])) {
    tau <- object[["tau"]][[h]]

    for (r in seq_len(ncol(tau))) {
      label <- {
        if (length(object[["tau"]]) > 1L) {
          sprintf("%s [%s]", colnames(tau)[r], names(object[["tau"]])[h])
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
  cli_cat("{.strong Generalized BART}")
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
    cli_cat("{.strong Nuisance parameters}")

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
    cli_cat("{.strong Random-effect scales}")
    cli_cat("{.emph Standard deviation of the group intercepts.}")
    print(round(x[["tau"]], digits))
  }

  cli::cat_line()
  cli_cat("{.strong Predictor usage}")
  cli_cat("{.emph Splitting rules per draw, and how often used at all.}")

  for (h in seq_along(x[["usage"]])) {
    if (length(x[["usage"]]) > 1L) {
      cli::cat_line()
      cli_cat("Predictor {.val {names(x$usage)[h]}}:")
    }
    print(round(x[["usage"]][[h]], digits))
  }

  invisible(x)
}

#' How often each predictor is used
#'
#' Reports, for every predictor, how many splitting rules the forest spends on
#' it and how often it is used at all, which is what "variable importance" means
#' for a BART model. This is the table [summary.bartisan_fit()] prints, returned
#' as a data frame rather than displayed (i.e., ready to sort, filter, or plot).
#'
#' @param object a `<bartisan_fit>` object; the output of a call to [bartisan()].
#' @param level `numeric`; the width of the interval reported for `splits`.
#'   Default is .95 for a 95% interval.
#'
#' @returns
#' A data frame with one row per predictor, sorted by `prop_used` and then
#' `splits`, both decreasing, and with the following columns:
#' * `variable`: the predictor, under the name the formula gave it
#' * `splits`: the mean number of splitting rules per draw that use it
#' * `splits_lower` and `splits_upper`: the endpoints of the `level` interval on
#'   that number
#' * `prop_used`: the proportion of draws in which it received at least one rule
#'
#' A family with more than one additive predictor has a forest for each, and the
#' data frame then gains a leading `predictor` column naming which.
#'
#' @details
#' ## Which Column Answers Which Question
#'
#' `prop_used` (the proportion of posterior draws in which the predictor
#' received at least one splitting rule) is the one to read first. It behaves
#' like a posterior probability that the predictor belongs in the model, and it
#' separates signal from noise sharply once `sparsity = TRUE` in
#' [bartisan_control()], which puts a Dirichlet prior on how the rules are shared
#' out and lets unused predictors be dropped rather than merely used rarely.
#'
#' `splits`, the mean number of rules per draw, says how much of the forest's
#' structure a predictor accounts for. It is the more familiar number and the
#' easier one to over-read.
#'
#' ## Three Things This Is Not
#'
#' **It is not an effect size.** A predictor can be split on constantly and move
#' the prediction very little, and the reverse happens too. If the question is
#' how much a predictor moves the outcome, that is a job for
#' \pkgfun{marginaleffects}{avg_comparisons} on the fitted model, not for this
#' table. See [bartisan-marginaleffects].
#'
#' **It is not stable under correlated predictors.** When two predictors carry
#' the same information the trees split on whichever is convenient, and the usage
#' distributes between them more or less arbitrarily. A predictor can matter and
#' still show a low `prop_used` because a collinear partner absorbed it, so a
#' group of correlated predictors is best read as a group.
#'
#' **It is not causal.** A ranking of predictors by usage is a description of
#' this fitted function, not of what would happen if any of them were changed.
#'
#' ## Reading It as Variable Selection
#'
#' With `sparsity = TRUE`, `prop_used` is usable as a selection rule: predictors
#' the forest genuinely needs sit near 1 and the rest fall near 0, usually with a
#' wide gap rather than a continuum. No threshold is correct in general, so the
#' gap is the thing to look at, and a conclusion worth reporting will not depend
#' on where in it the cut is made.
#'
#' @seealso [summary.bartisan_fit()], which prints the same table;
#'   [bartisan_control()] for `sparsity`; [bartisan-marginaleffects] for effects
#'   rather than usage
#'
#' @examples
#' data("rhc")
#' set.seed(123)
#'
#' # The sparsity prior concentrates the splitting rules on the predictors that
#' # earn them, which is what makes `prop_used` readable as a selection rule
#' fit <- bartisan(death ~ . - days, data = rhc, num_trees = 10,
#'                 num_burn = 50, num_draws = 50, sparsity = TRUE,
#'                 verbose = FALSE)
#'
#' imp <- variable_importance(fit)
#' imp
#'
#' # The predictors the forest reaches for in nearly every draw
#' subset(imp, prop_used > .9)
#'
#' @export
variable_importance <- function(object, level = 0.95) {

  if (!inherits(object, "bartisan_fit")) {
    arg::err("{.arg object} must be a fit from {.fn bartisan}")
  }

  arg::arg_number(level)
  arg::arg_between(level, c(0, 1), inclusive = FALSE)

  counts <- object[["counts"]]

  if (is_null(counts)) {
    arg::err("this fit carries no splitting counts")
  }

  rows <- lapply(names(counts), function(nm) {
    m <- counts[[nm]]
    summarized <- apply(m, 2L, post_summary, level = level)

    data.frame(predictor = nm,
               variable = colnames(m),
               splits = summarized["mean", ],
               splits_lower = summarized["lower", ],
               splits_upper = summarized["upper", ],
               # The proportion of draws in which the predictor was used at all,
               # which is the variable-selection reading; see the details.
               prop_used = colMeans(m > 0),
               row.names = NULL)
  })

  out <- do.call(rbind, rows)

  # One forest needs no column saying which forest.
  if (length(counts) == 1L) {
    out[["predictor"]] <- NULL
  }

  out <- out[order(out[["prop_used"]], out[["splits"]], decreasing = TRUE), , drop = FALSE]

  rownames(out) <- NULL

  out
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
#' The control function is not among them. It is the surface at the value each
#' covariate was centered on, which is a prediction rather than a coefficient;
#' `predict(object)` is what reports predictions.
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
    arg::err(c("this model has no varying coefficients, and a forest has no
                coefficient vector",
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
