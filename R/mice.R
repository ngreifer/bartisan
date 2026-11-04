#' Multiple imputation with generalized BART
#'
#' @description
#' `mice.impute.bartisan()` imputes one incomplete variable by fitting a
#' [bartisan()] model to the rows where it is observed and drawing each missing
#' value from the posterior predictive distribution of that model. It is called
#' by \pkgfun{mice}{mice} when `"bartisan"` is named in its `method` argument and
#' is not usually called directly.
#'
#' @param y the incomplete variable, supplied by \pkg{mice}. It can be numeric,
#'   logical, or a factor, ordered or not.
#' @param ry `logical`; which values of `y` are observed and so used to fit the
#'   model, supplied by \pkg{mice}.
#' @param x a numeric matrix of predictors with no missing values, one row per
#'   element of `y`, supplied by \pkg{mice}. Factors among the predictors arrive
#'   here already coded as dummy variables.
#' @param wy `logical`; which values of `y` to impute, supplied by \pkg{mice}.
#'   Default is `NULL` to impute those where `ry` is `FALSE`.
#' @param family optional; the family of the imputation model, in any of the
#'   forms the `family` argument of [bartisan()] accepts. Default is `NULL` to
#'   choose one from `y`; see Details.
#' @param ... further arguments passed to [bartisan_control()], such as
#'   `num_trees` or `num_burn`, or a `control` argument holding the output of a
#'   call to [bartisan_control()]. Arguments that are not settings of the model
#'   are dropped, because \pkg{mice} passes every imputation method the same
#'   arguments, including some meant for other methods. See Details.
#'
#' @returns
#' A vector with one imputed value for each `TRUE` element of `wy`, of the same
#' type as `y`: a factor with the levels of `y` for a factor, a logical vector
#' for a logical one, and a numeric vector otherwise.
#'
#' @details
#' Each missing value is drawn from the posterior predictive distribution of the
#' imputation model at one draw from its posterior, chosen at random from those
#' retained, so the imputations reflect uncertainty in the regression function
#' and in the parameters of the family as well as the residual variation. This
#' is the draw a proper multiple imputation requires (Rubin, 1987), and it
#' differs from methods such as `"cart"` and `"rf"`, which draw from the
#' observed values in a leaf rather than from a fitted distribution. The model is
#' refit from scratch for every variable at every iteration of every
#' imputation, as for the other imputation methods in \pkg{mice}.
#'
#' When `family` is `NULL`, the family is chosen from `y` by the rules
#' [bartisan()] uses when its own `family` argument is left out: [ordinal()] for
#' an ordered factor, [binomial()] for a factor with two levels, a logical
#' vector, or a numeric vector of zeros and ones, [multinomial()] for any other
#' factor, and [dpm()] for any other numeric variable. The [dpm()] family models
#' the error distribution as a mixture, so the imputations of a continuous
#' variable keep the skewness, heavy tails, or multimodality its residuals have
#' rather than being made Gaussian. A count or a non-negative variable may be
#' better imputed with a family that respects its support, such as `poisson()`
#' or `Gamma()`; since the default for a count is [dpm()], its imputations will
#' not be whole numbers. A family can be given for one variable alone through
#' the `blots` argument of \pkgfun{mice}{mice} (e.g.,
#' `blots = list(visits = list(family = poisson()))`).
#'
#' Settings of the model are passed the same way, through `blots` for one
#' variable or through `...` of \pkgfun{mice}{mice} for every variable imputed
#' with this method. Because only one posterior draw is used per imputation,
#' `num_draws` is 1 unless it is supplied (directly or in `control`), which
#' makes each fit several times faster than one at the default of
#' [bartisan_control()] without affecting the distribution of the imputations;
#' `num_burn` is left at its default, since the burn-in period is the reason the
#' one retained draw can be taken as a draw from the posterior.
#'
#' \pkg{mice} finds an imputation method by its name, so \pkg{bartisan} must be
#' attached with [library()] rather than only loaded. A variable that has no
#' predictors in the `predictorMatrix` of \pkgfun{mice}{mice} cannot be imputed
#' with this method, since a forest needs something to split on.
#'
#' @seealso
#' \pkgfun{mice}{mice} for multiple imputation by chained equations;
#' [bartisan()] for the imputation model and the settings that can be passed
#' through `...`; [SL.bartisan()] and [`bartisan-parsnip`] for the other
#' prediction workflows \pkg{bartisan} fits into; `vignette("workflows")` for
#' a worked example of all three
#'
#' @references
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*. John
#' Wiley & Sons. \doi{10.1002/9780470316696}
#'
#' @examplesIf rlang::is_installed("mice")
#' library("mice")
#' data("rhc")
#' set.seed(123)
#'
#' # Make some values of a skewed continuous variable and
#' # of a binary one missing
#' d <- rhc[c("death", "age", "sex", "aps", "meanbp", "crea")]
#' is.na(d$crea[sample(nrow(d), 300)]) <- TRUE
#' is.na(d$sex[sample(nrow(d), 300)]) <- TRUE
#'
#' # Impute both with bartisan, using a smaller model to
#' # keep the example fast
#' imp <- mice(d, method = "bartisan", m = 3, maxit = 2,
#'             num_trees = 10, num_burn = 50,
#'             printFlag = FALSE)
#'
#' # Pool an analysis of the imputed data sets
#' fits <- with(imp, glm(death ~ age + sex + aps + crea,
#'                       family = binomial))
#' summary(pool(fits))
#'
#' @export mice.impute.bartisan
mice.impute.bartisan <- function(y, ry, x, wy = NULL, family = NULL, ...) {
  wy <- wy %or% !ry

  if (ncol(x) == 0L) {
    arg::err(c("{.fn mice.impute.bartisan} needs at least one predictor to
                fit a forest to.",
               i = "Give the variable at least one predictor in the
                    {.arg predictorMatrix} of {.fn mice::mice}."))
  }

  x <- as.data.frame(x)

  # The response gets a name that cannot collide with a predictor's.
  response <- make.unique(c(names(x), "y"))[ncol(x) + 1L]

  observed <- x[ry, , drop = FALSE]
  observed[[response]] <- y[ry]

  # Chosen here rather than left to `bartisan()` so that the message it prints
  # when it chooses does not repeat for every variable at every iteration; the
  # rules are the same and are stated on the help page instead.
  family <- family %or% suppressMessages(default_family(y[ry]))

  formula <- stats::reformulate(".", response = response)

  fit <- bartisan(formula,
                  data = observed, family = family,
                  control = imputation_control(list(...)))

  iteration <- sample.int(nrow(fit[["sigma_mu"]]), 1L)

  posterior_predict.bartisan_fit(fit,
                                 newdata = x[wy, , drop = FALSE],
                                 iterations = iteration) |>
    as.vector() |>
    imputed_values(y, fit[["family"]][["family"]])
}

# What `mice()` passes along includes its own settings and those of every other
# method (`eps`, `donors`, `ridge`, `type`, ...), so the model's settings are
# picked out by name and the rest dropped rather than reported as errors.
imputation_control <- function(dots) {
  control <- dots[["control"]] %or% bartisan_control()

  if (!inherits(control, "bartisan_control")) {
    arg::err("{.arg control} must be the result of {.fn bartisan_control}")
  }

  settings <- dots[names(dots) %in% names(formals(bartisan_control))]

  if (!"num_draws" %in% c(names(attr(control, "supplied")), names(settings))) {
    settings[["num_draws"]] <- 1L
  }

  merge_control(control, settings)
}

# A replicate outcome comes back on the scale the likelihood was written on:
# 0/1 for a binomial response, 1-based category codes for one with more
# categories, and the value itself otherwise. mice wants the type of `y` back.
imputed_values <- function(draw, y, family) {
  if (is.factor(y)) {
    codes <- {
      if (identical(family, "binomial")) draw + 1L
      else draw
    }

    return(factor(levels(y)[codes], levels = levels(y),
                  ordered = is.ordered(y)))
  }

  if (is.logical(y)) {
    return(draw == 1)
  }

  draw
}
