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
#' @param draws `logical`; whether to return the splitting counts of every
#'   posterior draw rather than a summary of them, which is what to use for a
#'   comparison the summary does not offer (e.g., the posterior probability that
#'   one predictor takes more rules than another). Default is `FALSE` to return
#'   the summary. Cannot be combined with `plot`.
#' @param plot `logical`; whether to return a plot of the table instead of the
#'   table, as a \CRANpkg{ggplot2} object showing `splits` and its interval for
#'   each predictor. Default is `FALSE`. Equivalent to calling `plot()` on the
#'   result, which is usually the more convenient of the two since the table can
#'   be subset first (e.g., `plot(head(imp, 10))` for a wide model).
#'   \pkg{ggplot2} must be installed for either.
#'
#' @returns
#' A `<bartisan_importance>` object, which is a data frame with one row per
#' predictor and its own `print()` and `plot()` methods, sorted by `prop_used`
#' and then
#' `splits`, both decreasing, and with the following columns:
#' * `variable`: the predictor, under the name the formula gave it
#' * `prop_used`: the proportion of draws in which it received at least one rule
#' * `prop_splits`: its share of all the splitting rules in the forest, averaged
#'   over draws
#' * `splits`: the mean number of splitting rules per draw that use it
#' * `splits_lower` and `splits_upper`: the endpoints of the `level` interval on
#'   that number, which the `print()` method leaves out of the displayed table
#'
#' A family with more than one additive predictor has a forest for each, and the
#' data frame then gains a leading `predictor` column naming which.
#'
#' With `draws = TRUE`, the draws-by-predictors matrix of counts instead, or a
#' named list of them when there is more than one forest. With `plot = TRUE`, a
#' \pkg{ggplot2} object.
#'
#' @details
#' ## Which Column Answers Which Question
#'
#' `prop_used` (the proportion of posterior draws in which the predictor
#' received at least one splitting rule) is the one to read first. It behaves
#' like a posterior probability that the predictor belongs in the model, and it
#' separates signal from noise sharply once `sparsity = TRUE` in
#' [bartisan_control()], which puts a Dirichlet prior on how the rules are shared
#' out and lets unused predictors be dropped rather than merely used rarely. With
#' `sparsity = FALSE` every predictor keeps a share of the rules and `prop_used`
#' sits near 1 throughout, so the `print()` method says so rather than leaving
#' the column to be misread.
#'
#' `prop_splits` is the one to reach for when two fits are being compared.
#' `splits` counts rules, so it scales with `num_trees` and a fit of 40 trees
#' will report four times the count of a fit of 10 without being four times as
#' informative; the share is computed within each draw before averaging, so the
#' shares add to one whatever the forest size. Note that the share is a property
#' of the forest and not only of the predictor, so it can still shift when the
#' forest size changes how the rules are allocated; the ranking is the stable
#' part.
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
#' wide gap rather than a continuum. It is the posterior inclusion probability
#' that the Dirichlet prior was introduced to make readable (Linero, 2018), and
#' cutting it at .5 gives the median probability model, which Barbieri and Berger
#' (2004) show is often a better predictive submodel under squared error loss
#' than the model of highest posterior probability. That is the same quantity and
#' the same cut \pkg{SoftBart} reports from a soft BART fit, computed from the
#' same splitting counts. No threshold is correct in general, though, so the gap
#' is the thing to look at, and a conclusion worth reporting will not depend on
#' where in it the cut is made.
#'
#' The .5 cut chooses a predictive submodel and is not a test, and it does not
#' behave like one. On data where no predictor matters at all, at \eqn{n = 300}
#' with 25 predictors and 20 trees, it selected 7 of the 25 on average with
#' `sparsity = TRUE` and all 25 with `sparsity = FALSE`. A forest asked to fit
#' noise still puts its rules somewhere, and `prop_used` reports where they went
#' rather than whether they were needed.
#'
#' @seealso [summary.bartisan_fit()], which prints the same table;
#'   [bartisan_control()] for `sparsity`; [bartisan-marginaleffects] for effects
#'   rather than usage
#'
#' @references
#' Barbieri, M. M., & Berger, J. O. (2004). Optimal predictive model selection.
#' *The Annals of Statistics*, 32(3). \doi{10.1214/009053604000000238}
#'
#' Bleich, J., Kapelner, A., George, E. I., & Jensen, S. T. (2014). Variable
#' selection for BART: an application to gene regulation. *The Annals of Applied
#' Statistics*, 8(3). \doi{10.1214/14-AOAS755}
#'
#' Linero, A. R. (2018). Bayesian regression trees for high-dimensional
#' prediction and variable selection. *Journal of the American Statistical
#' Association*, 113(522), 626--636. \doi{10.1080/01621459.2016.1264957}
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
#' # `prop_splits` is the column that survives a change of forest size, since
#' # the shares add to one however many rules there are to share
#' big <- bartisan(death ~ . - days, data = rhc, num_trees = 40,
#'                 num_burn = 50, num_draws = 50, sparsity = TRUE,
#'                 verbose = FALSE)
#'
#' merge(variable_importance(fit)[c("variable", "prop_splits")],
#'       variable_importance(big)[c("variable", "prop_splits")],
#'       by = "variable", suffixes = c("_10", "_40"))
#'
#' # The counts themselves, for a comparison the summary does not make
#' counts <- variable_importance(fit, draws = TRUE)
#' mean(counts[, "aps"] > counts[, "meanbp"])
#'
#' # The ranking, drawn. Subsetting first is what keeps a wide model readable.
#' if (rlang::is_installed("ggplot2")) {
#'   plot(head(imp, 8))
#' }
#'
#' @export
variable_importance <- function(object, level = 0.95, draws = FALSE,
                                plot = FALSE) {

  arg::arg_is(object, "bartisan_fit")
  arg::arg_number(level)
  arg::arg_between(level, c(0, 1), inclusive = FALSE)
  arg::arg_flag(draws)
  arg::arg_flag(plot)

  counts <- object[["counts"]]

  if (is_null(counts)) {
    arg::err("this fit carries no splitting counts")
  }

  if (draws) {
    if (plot) {
      arg::err(c("{.arg draws} and {.arg plot} ask for different things.",
                 i = "{.code draws = TRUE} returns the counts to summarize
                      yourself; {.code plot = TRUE} draws the summary."))
    }

    # One forest needs no list around it, matching what `coef()` does.
    if (length(counts) == 1L) {
      return(counts[[1L]])
    }

    return(counts)
  }

  rows <- lapply(names(counts), function(nm) {
    m <- counts[[nm]]
    summarized <- apply(m, 2L, post_summary, level = level)

    # The share of the forest's rules, draw by draw before averaging, so that
    # the shares add to one within a draw. This is the column to compare across
    # fits: `splits` counts rules and so moves with `num_trees`, and a fit of
    # 200 trees is not more informative than one of 50 for having four times as
    # many of them.
    total <- rowSums(m)
    share <- m / ifelse(total > 0, total, NA_real_)

    data.frame(predictor = nm,
               variable = colnames(m),
               # The proportion of draws in which the predictor was used at
               # all, which is the variable-selection reading; see the details.
               prop_used = colMeans(m > 0),
               prop_splits = colMeans(share, na.rm = TRUE),
               splits = summarized["mean", ],
               splits_lower = summarized["lower", ],
               splits_upper = summarized["upper", ],
               row.names = NULL)
  })

  out <- do.call(rbind, rows)

  # One forest needs no column saying which forest.
  if (length(counts) == 1L) {
    out[["predictor"]] <- NULL
  }

  out <- out[order(out[["prop_used"]], out[["splits"]], decreasing = TRUE), ,
             drop = FALSE]

  rownames(out) <- NULL

  # `sparsity` is recorded because it decides whether `prop_used` can be read as
  # a selection rule at all, and a reader looking at the table has no other way
  # to see which it was.
  attr(out, "level") <- level
  attr(out, "sparsity") <- isTRUE(object[["control"]][["sparsity"]])
  class(out) <- c("bartisan_importance", "data.frame")

  if (!plot) {
    return(out)
  }

  # Through the method rather than beside it, so that `plot = TRUE` and `plot()`
  # cannot become two drawings of the same table.
  plot(out)
}

#' @rdname variable_importance
#' @param x a `<bartisan_importance>` object; the output of a call to
#'   `variable_importance()`.
#' @param y not used.
#' @param ... not used.
#' @export
plot.bartisan_importance <- function(x, y, ...) {
  importance_plot(x)
}

# The ranking, drawn. Kept separate from the summary so that a caller who wants
# a different picture has the data frame to draw it from, which is the same
# split `error_density()` makes.
importance_plot <- function(x) {
  require_ggplot2("variable importance")

  x[["variable"]] <- factor(x[["variable"]],
                            levels = rev(unique(x[["variable"]])))

  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data$splits, y = .data$variable)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$splits_lower,
                                        xmax = .data$splits_upper),
                           orientation = "y", width = 0, color = "grey60") +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Splitting rules per draw", y = NULL) +
    ggplot2::theme_bw()

  if (!is_null(x[["predictor"]])) {
    p <- p + ggplot2::facet_wrap(~ .data$predictor, scales = "free_y")
  }

  p
}

#' @export
print.bartisan_importance <- function(x, digits = 3L, ...) {
  level <- attr(x, "level")
  sparse <- attr(x, "sparsity")

  cli_cat("{.underline Variable importance}")
  cli::cat_line()

  show <- as.data.frame(x)

  # The interval is dropped from the printed table rather than from the object:
  # four numbers per predictor is more than a ranking needs, and the columns are
  # still there for anyone who wants them.
  show[c("splits_lower", "splits_upper")] <- NULL

  for (nm in c("prop_used", "prop_splits")) {
    show[[nm]] <- round(show[[nm]], digits)
  }

  show[["splits"]] <- round(show[["splits"]], 1L)

  print(show, row.names = FALSE)
  cli::cat_line()

  if (isFALSE(sparse)) {
    cli_bullets_cat(c(i = "Fitted with {.code sparsity = FALSE}, so every
                          predictor keeps a share of the rules and
                          {.field prop_used} is near 1 throughout. Refit with
                          {.code sparsity = TRUE} to read it as a selection
                          rule."))
  }

  if (!is_null(level)) {
    cli_bullets_cat(c(i = "{.field splits_lower} and {.field splits_upper} hold
                          the {level * 100}% interval, not shown above."))
  }

  invisible(x)
}
