# Bayesian causal forests, as a wrapper over the varying-coefficient interface.
#
# Everything here is expressible in `bartisan()` with a `vc()` term, which is
# deliberate: a wrapper that needed something the general interface could not say
# would mean the general interface was wrong.

#' Bayesian causal forests
#'
#' A varying-coefficient model set up for estimating a treatment effect: a
#' control function for the outcome under no treatment and a separate forest for
#' the effect, with the prior on the effect regularized more heavily than the
#' prior on the control function. This is Hahn, Murray and Carvalho (2020) for a
#' binary treatment and Woody, Carvalho, Hahn and Murray (2020) for a continuous
#' one.
#'
#' @param formula a model formula. The right-hand side lists the covariates; the
#'   treatment is named in `treatment` rather than here, and is removed from the
#'   covariates if it appears among them.
#' @param treatment a one-sided formula naming the treatment, as in `~ z`.
#' @param data a data frame.
#' @param family the outcome distribution, as in [bartisan()].
#' @param moderators a one-sided formula naming the covariates the treatment
#'   effect may vary with. The default, `NULL`, is all of them.
#' @param propensity what to do about the probability of treatment. `TRUE`, the
#'   default, fits a model for it and adds the fitted values to the **control
#'   function only**; `FALSE` fits nothing; a numeric vector or matrix is used as
#'   given; a one-sided formula fits it with that model rather than the outcome's
#'   covariates.
#' @param propensity_args a list of [bartisan_control()] settings for the
#'   propensity model, which is a prediction problem and so keeps the sparsity
#'   prior the outcome model turns off.
#' @param ... passed to [bartisan()], including [bartisan_control()] settings.
#'
#' @returns
#' A `bartisan_fit`, with the treatment's coefficient forest named for the
#' treatment. [coef()] gives the conditional effect for each observation and
#' `marginaleffects::avg_comparisons()` the average.
#'
#' @details
#' # What the wrapper decides
#'
#' Four things, all of which can be written out in [bartisan()] directly:
#'
#' * The treatment gets a `vc()` term, so the effect is a forest of its own with
#'   its own prior rather than whatever difference a single forest with the
#'   treatment among its predictors happens to produce.
#' * The propensity score goes in the control function and **not** in the effect
#'   forest. That is Hahn et al.'s recommendation and the flag `bcf`, `stochtree`
#'   and this package all provide; the point is to let the control function
#'   absorb the selection without letting the effect vary with it.
#' * The effect forest gets fewer trees than the control function, since patterns
#'   of effect heterogeneity are usually simpler than prognostic surfaces.
#' * `sparsity = FALSE` on the outcome model. A variable-selection prior on the
#'   variable whose contrast is the estimand puts a point mass at exactly zero in
#'   the posterior of the effect; see the measurements in [bartisan_control()].
#'   The propensity model keeps the default, because predicting who was treated
#'   is a prediction problem.
#'
#' # The treatment's type
#'
#' The treatment decides the model for the propensity score and what that score
#' even is.
#'
#' | treatment | propensity score | model |
#' | --- | --- | --- |
#' | binary | one column, the probability of treatment | [binomial()] |
#' | `K` categories | the whole vector of assignment probabilities | [multinomial()] |
#' | continuous | a conditional density, not a regression | not fitted |
#'
#' For a treatment with more than two categories the balancing score is the
#' *vector* of assignment probabilities (Imbens, 2000; Imai and van Dyk, 2004),
#' not any one of them, so all of them go into the control function. They sum to
#' one and are therefore collinear, which costs a tree ensemble nothing.
#'
#' For a continuous treatment the analogue is the conditional density of the
#' treatment given the covariates evaluated at the observed dose (Hirano and
#' Imbens, 2004), which needs a density model rather than a regression, so
#' `propensity = TRUE` is refused and a score supplied as a number is used as
#' given.
#'
#' # The assumption a continuous treatment carries
#'
#' With a continuous treatment this fits `f0(x) + z * f1(x)`: a dose response
#' that is **linear in the dose**, with a slope that varies. For a binary
#' treatment that is no assumption at all. For a continuous one it is a real one,
#' and it is the assumption Woody et al. (2020) make and diagnose. If the dose
#' response itself might be curved, either put the treatment in as an ordinary
#' predictor or give `f1` the treatment among its moderators, which makes the
#' effect vary across the dose; see [vc()].
#'
#' # Predicting for new data
#'
#' The propensity score is a predictor of the control function, and it is one the
#' caller never named, so `newdata` taken from their own frame does not carry it.
#' `predict()` rebuilds it from the model this kept, which is what makes that
#' work at all.
#'
#' One consequence is worth knowing. A predictor goes through the quantile
#' transform, which is a step function, so a score that rebuilds to within 1e-11
#' can still land on the other side of a step: predictions for the *training*
#' data come back to within a few percent of the response's spread rather than
#' exactly. Passing the stored score in `newdata` -- it is in
#' `fit$bcf$propensity` -- removes the reconstruction and reproduces the fit to
#' machine precision. Supplying `propensity` as a number rather than fitting it
#' has the same effect, and then `newdata` must carry the column.
#'
#' @seealso [bartisan()] and [vc()] for the general interface this is written in
#'   terms of, and `vignette("causal")`.
#'
#' @references
#' Hahn, P. R., Murray, J. S., & Carvalho, C. M. (2020). Bayesian regression tree
#' models for causal inference: regularization, confounding, and heterogeneous
#' effects. *Bayesian Analysis*, 15(3), 965--1056. \doi{10.1214/19-BA1195}
#'
#' Woody, S., Carvalho, C. M., Hahn, P. R., & Murray, J. S. (2020). Estimating
#' heterogeneous effects of continuous exposures using Bayesian tree ensembles.
#' \doi{10.48550/arXiv.2007.09845}
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' d$z <- rbinom(n, 1, plogis(d$x1))
#' d$y <- d$x1 + d$z * (1 + d$x2) + rnorm(n)
#'
#' fit <- bcf(y ~ x1 + x2, treatment = ~ z, data = d, family = gaussian(),
#'            num_trees = c(10, 5), num_burn = 50, num_draws = 50,
#'            verbose = FALSE)
#'
#' head(coef(fit))
#'
#' @export
bcf <- function(formula, treatment, data, family = NULL, moderators = NULL,
                propensity = TRUE, propensity_args = list(), ...) {

  arg::arg_formula(formula, one_sided = FALSE)
  arg::arg_formula(treatment, one_sided = TRUE)

  name <- all.vars(treatment)

  if (length(name) != 1L) {
    arg::err("{.arg treatment} must name exactly one variable, as in
              {.code ~ z}")
  }

  if (!name %in% names(data)) {
    arg::err("{.arg data} has no column {.val {name}}")
  }

  if (!is_null(moderators)) {
    arg::arg_formula(moderators, one_sided = TRUE)
  }

  # The treatment is not one of the covariates. Removing it rather than
  # complaining is the friendly reading of `y ~ .`, which is how most callers
  # will write the covariates.
  covariates <- setdiff(
    attr(stats::terms(formula, data = data), "term.labels"), name)

  if (length(covariates) == 0L) {
    arg::err("{.arg formula} must name at least one covariate besides the
              treatment")
  }

  scored <- bcf_propensity(propensity, name, covariates, data, propensity_args)
  score <- scored[["score"]]

  if (!is_null(score)) {
    data <- cbind(data, score)
    covariates <- c(covariates, colnames(score))
  }

  # The propensity score reaches the control function and not the effect forest,
  # which is the whole point of estimating it: it lets the control function
  # absorb the selection without letting the effect vary with it.
  moderator_terms <- {
    if (is_null(moderators)) {
      setdiff(covariates, colnames(score))
    }
    else {
      attr(stats::terms(moderators, data = data), "term.labels")
    }
  }

  outcome <- stats::reformulate(
    c(covariates,
      sprintf("vc(%s, ~ %s)", name, paste(moderator_terms, collapse = " + "))),
    response = formula[[2L]])
  environment(outcome) <- environment(formula)

  dots <- list(...)

  # Fewer trees for the effect than for the control function, since patterns of
  # effect heterogeneity are usually simpler than prognostic surfaces, and no
  # variable-selection prior, since the estimand is a contrast on the treatment.
  #
  # A categorical treatment gets one effect forest per level, so the count has to
  # be built from the treatment rather than assumed to be two.
  if (is_null(dots[["num_trees"]])) {
    z <- data[[name]]
    effects <- if (is.factor(z) || is.character(z)) {
      length(levels(as.factor(z)))
    } else {
      1L
    }

    dots[["num_trees"]] <- c(50L, rep.int(25L, effects))
  }

  if (is_null(dots[["sparsity"]])) {
    dots[["sparsity"]] <- FALSE
  }

  out <- do.call(bartisan,
                 c(list(formula = outcome, data = data, family = family), dots))

  # The propensity model is kept so that `predict()` on data the caller has can
  # rebuild the score. Without it the fit depends on a column the caller's data
  # does not have, and every prediction on new data fails on a missing variable
  # whose name they never chose.
  out[["bcf"]] <- list(treatment = name, propensity = score,
                       model = scored[["model"]], moderators = moderator_terms)
  out
}

# The propensity columns for new data, from the model the fit kept. `NULL` when
# there is nothing to add, either because no score was used or because the data
# already carries it.
bcf_newdata_score <- function(object, newdata) {
  spec <- object[["bcf"]]

  if (is_null(spec) || is_null(spec[["propensity"]])) {
    return(NULL)
  }

  wanted <- colnames(spec[["propensity"]])

  if (all(wanted %in% names(newdata))) {
    return(NULL)
  }

  if (is_null(spec[["model"]])) {
    arg::err(c("{.arg newdata} is missing {length(wanted)} column{?s} this fit
                uses as a predictor: {.val {wanted}}",
               i = "The propensity score was supplied rather than fitted, so it
                  cannot be rebuilt. Add the same column to {.arg newdata}."))
  }

  score <- as.matrix(stats::fitted(spec[["model"]], newdata = newdata))
  colnames(score) <- wanted
  score
}

# The propensity score, or nothing. The model follows the treatment's type,
# because what the score *is* follows the treatment's type.
bcf_propensity <- function(propensity, name, covariates, data, args) {
  if (isFALSE(propensity) || is_null(propensity)) {
    return(list(score = NULL, model = NULL))
  }

  z <- data[[name]]
  kind <- treatment_kind(z)

  if (is.numeric(propensity) || is.matrix(propensity)) {
    score <- as.matrix(propensity)

    if (nrow(score) != nrow(data)) {
      arg::err("{.arg propensity} must have one row per observation")
    }

    colnames(score) <- bcf_score_names(ncol(score))
    return(list(score = score, model = NULL))
  }

  if (identical(kind, "continuous")) {
    arg::err(c("a continuous treatment has no propensity score that is a
                probability",
               i = "Its analogue is the conditional density of the treatment
                  given the covariates at the observed dose (Hirano and Imbens,
                  2004), which needs a density model rather than a regression.",
               i = "Pass {.code propensity = FALSE}, or supply one as a numeric
                  vector."))
  }

  terms <- {
    if (rlang::is_formula(propensity)) {
      attr(stats::terms(propensity, data = data), "term.labels")
    }
    else {
      covariates
    }
  }

  model <- stats::reformulate(terms, response = as.symbol(name))
  fit_family <- if (identical(kind, "binary")) stats::binomial() else multinomial()

  # The propensity model keeps the sparsity prior: predicting who was treated is
  # a prediction problem, and nothing here is a contrast.
  fit <- do.call(bartisan,
                 c(list(formula = model, data = data, family = fit_family),
                   args))

  score <- as.matrix(stats::fitted(fit))
  colnames(score) <- bcf_score_names(ncol(score))
  list(score = score, model = fit)
}

treatment_kind <- function(z) {
  if (is.factor(z) || is.character(z)) {
    return(if (length(unique(stats::na.omit(z))) == 2L) "binary" else "categorical")
  }

  values <- unique(stats::na.omit(z))

  if (length(values) == 2L && all(values %in% c(0, 1))) {
    return("binary")
  }

  "continuous"
}

bcf_score_names <- function(k) {
  if (k == 1L) ".propensity" else sprintf(".propensity%d", seq_len(k))
}
