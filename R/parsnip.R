#' Generalized BART in tidymodels
#'
#' @name bartisan-parsnip
#' @usage NULL
#'
#' @description
#' This page explains the details of fitting a [bartisan()] model as a
#' \CRANpkg{parsnip} model specification, which lets it be used in
#' \pkg{tidymodels} workflows: preprocessed with a recipe, combined with one in a
#' workflow, resampled, and tuned. Setting `engine = "bartisan"` in a call to
#' \pkgfun{parsnip}{bart}, or calling `set_engine("bartisan")` on its output,
#' selects it. The engine supports the `"regression"` and `"classification"`
#' modes, and it adds the `"censored regression"` mode, which
#' \pkgfun{parsnip}{bart} does not otherwise have, for the survival families.
#'
#' @section Tuning Parameters:
#' The model arguments of \pkgfun{parsnip}{bart} correspond to arguments of
#' [bartisan_control()] as follows, and each can be tuned with
#' \pkgfun{tune}{tune}:
#' \describe{
#'   \item{`trees`}{`num_trees`, the number of trees. Default is 50.}
#'   \item{`prior_terminal_node_coef`}{`gamma`, the base of the probability that
#'     a node at depth `d` splits, `gamma * (1 + d)^(-beta)`. Default is .95.}
#'   \item{`prior_terminal_node_expo`}{`beta`, the exponent in that probability.
#'     Default is 2.}
#'   \item{`prior_outcome_range`}{`k`, which sets the prior standard deviation
#'     of the forest to `3 / k` times the natural scale of its additive
#'     predictor, so that larger values shrink the fit harder toward the
#'     intercept-only model. Default is 2. Its tuning range is .5 to 5 rather
#'     than the \pkg{dials} default, which extends to 0, where the prior is so
#'     diffuse that the fit is of no use.}
#' }
#' Any other argument of [bartisan()] or [bartisan_control()] is an engine
#' argument and is supplied to \pkgfun{parsnip}{set_engine} (e.g.,
#' `set_engine("bartisan", family = gaussian(), num_draws = 400)`); these can be
#' tuned too.
#'
#' @section Families:
#' The family is chosen from the outcome by the rules [bartisan()] uses when its
#' `family` argument is left out, unless it is supplied as an engine argument. In
#' the `"regression"` mode, a numeric outcome is fit with [dpm()] by default; in
#' the `"classification"` mode, a factor with two levels is fit with [binomial()],
#' an ordered factor with [ordinal()], and any other factor with [multinomial()];
#' and in the `"censored regression"` mode, where the outcome is a
#' \pkgfun{survival}{Surv} object with right censoring, the outcome is fit with
#' [dpm_aft()]. [bartisan()] prints a message naming the family it chose each
#' time a model is fit, which during resampling or tuning means once per fit;
#' supplying the family as an engine argument silences it. Any family whose
#' response is compatible with the mode can be requested this way, such as
#' `poisson()` for a count in the `"regression"` mode or [ph()] in the
#' `"censored regression"` mode.
#' Note that \pkg{parsnip} evaluates engine arguments before passing them on, so a
#' family supplied as an object (e.g., `family = poisson()`) is stored in the
#' call of the fit in full, which makes it long to print; supplying the family by
#' name (e.g., `family = "poisson"`) gives the same model with a short call.
#'
#' @section Predictions:
#' The following prediction types are available through
#' \pkgfun{parsnip}{predict.model_fit}:
#' \describe{
#'   \item{`"numeric"`}{the posterior mean of the response (i.e.,
#'     [predict.bartisan_fit()] with `type = "response"`), in the `"regression"`
#'     mode.}
#'   \item{`"class"`}{the category with the highest posterior mean probability,
#'     in the `"classification"` mode.}
#'   \item{`"prob"`}{the posterior mean probability of each category, in the
#'     `"classification"` mode.}
#'   \item{`"conf_int"`}{a credible interval for the quantity reported by
#'     `"numeric"` or for each of those reported by `"prob"`, in either mode.}
#'   \item{`"pred_int"`}{a posterior predictive interval for the outcome, in the
#'     `"regression"` mode. This requires \pkg{rstantools}.}
#'   \item{`"time"`}{the posterior mean of the median survival time (i.e.,
#'     [predict.bartisan_fit()] with `type = "response"`), in the
#'     `"censored regression"` mode.}
#'   \item{`"survival"`}{the posterior mean of the survival function at the
#'     times given in `eval_time`, in the `"censored regression"` mode. The
#'     result has one row per row of `new_data`, each holding a tibble with
#'     columns `.eval_time` and `.pred_survival`; setting
#'     `interval = "confidence"` adds `.pred_lower` and `.pred_upper`, the limits
#'     of a credible interval for the survival probability.}
#'   \item{`"hazard"`}{the posterior mean of the hazard at the times given in
#'     `eval_time`, in the `"censored regression"` mode, laid out as for
#'     `"survival"` with the column `.pred_hazard`.}
#'   \item{`"linear_pred"`}{the posterior mean of the additive predictor, in the
#'     `"censored regression"` mode. This is the location of the log survival
#'     time for the accelerated failure time families and the negative of the log
#'     hazard ratio for [ph()], so that it increases with the survival time for
#'     every family, as \pkg{censored} reports it; setting `increasing = FALSE`
#'     reverses its sign. For [dpm_aft_ls()] and [lognormal_aft_ls()] it is the
#'     predictor of the location only, so it does not order the units by their
#'     survival where their scales differ.}
#'   \item{`"raw"`}{the posterior draws of the mean of the response, as a matrix
#'     with one row per draw and one column per row of `new_data` (or, in the
#'     `"classification"` mode, the draws of the category probabilities, as an
#'     array with a third dimension for the categories, and in the
#'     `"censored regression"` mode, the draws of the median survival time).}
#' }
#' \pkg{parsnip} allows an evaluation time of 0, at which [predict.bartisan_fit()]
#' reports nothing. There the survival probability is reported as 1, since every
#' unit has survived to it, and the hazard as `NA`, since for most of the
#' families it is not finite.
#'
#' The intervals, including those `interval = "confidence"` adds to the
#' `"survival"` predictions, are equal-tailed: their limits are the `(1 - level) / 2` and
#' `(1 + level) / 2` quantiles of the posterior draws, where `level` is the
#' argument of that name to \pkgfun{parsnip}{predict.model_fit} (default .95).
#' They are not highest posterior density intervals, which for a skewed
#' posterior, such as that of a probability near 0 or 1, would be narrower and
#' shifted toward the mode. Setting `std_error = TRUE` adds a `.std_error`
#' column holding the posterior standard deviation.
#'
#' The underlying `<bartisan_fit>` object is available with
#' \pkgfun{parsnip}{extract_fit_engine}, and every method in the package applies
#' to it, including [estimate_effect()], [variable_importance()], and
#' [diagnose()].
#'
#' @details
#' The engine is registered with \pkg{parsnip} when \pkg{bartisan} is loaded, or
#' when \pkg{parsnip} is loaded if \pkg{bartisan} was loaded first, so
#' \pkg{bartisan} needs to be loaded (e.g., with `library("bartisan")`) before
#' the specification is fit. When resampling or tuning in parallel,
#' \pkg{tune} loads it on the workers itself.
#'
#' The model is fit through the formula interface of [bartisan()], so factor
#' predictors reach it as factors rather than dummy variables, and its splitting
#' rules for categorical predictors apply to them; a recipe that does create
#' dummy variables is respected. Case weights, set with
#' \pkgfun{workflows}{add_case_weights}, are passed as prior weights, with the
#' consequences described for the `weights` argument of [bartisan()].
#'
#' The `"censored regression"` mode does not require \CRANpkg{censored}, which
#' supplies the other engines for that mode, but its predictions are laid out as
#' those of the engines there are. In particular, \pkgfun{parsnip}{augment} with
#' `eval_time` adds the censoring weights that the dynamic survival metrics of
#' \pkg{yardstick} (e.g., \pkgfun{yardstick}{brier_survival}) need, so those
#' metrics can be used to compare a [bartisan()] model with the models
#' \pkg{censored} fits, in resampling and in tuning.
#'
#' @seealso
#' \pkgfun{parsnip}{bart} for the model specification; [bartisan()] and
#' [bartisan_control()] for the engine arguments; [SL.bartisan()] and
#' [mice.impute.bartisan()] for the other prediction workflows \pkg{bartisan}
#' fits into
#'
#' @examplesIf rlang::is_installed(c("parsnip", "survival"))
#' library("parsnip")
#' library("bartisan")
#' data("rhc")
#' set.seed(123)
#'
#' rhc$died <- factor(rhc$death, labels = c("no", "yes"))
#'
#' # A classification model, with a smaller forest to keep
#' # the example fast and the family named to silence the
#' # message reporting the choice
#' spec <- bart(trees = 10, mode = "classification") |>
#'   set_engine("bartisan", family = "binomial",
#'              num_burn = 50, num_draws = 50)
#'
#' fit <- spec |>
#'   fit(died ~ age + sex + aps + meanbp + pafi + crea,
#'       data = rhc)
#'
#' # Probabilities and 90% credible intervals for them
#' predict(fit, rhc[1:5, ], type = "prob")
#' predict(fit, rhc[1:5, ], type = "conf_int", level = .9)
#'
#' # The bartisan fit, for anything beyond predictions
#' extract_fit_engine(fit)
#'
#' # A survival model for the time to death, with the
#' # survival probability at 30 and 180 days
#' surv_fit <- bart(trees = 10, mode = "censored regression") |>
#'   set_engine("bartisan", family = "weibull_aft",
#'              num_burn = 50, num_draws = 50) |>
#'   fit(survival::Surv(days, death) ~ age + sex + aps +
#'         meanbp + pafi + crea,
#'       data = rhc)
#'
#' pred <- predict(surv_fit, rhc[1:5, ], type = "survival",
#'                 eval_time = c(30, 180))
#' pred$.pred[[1]]
NULL

register_parsnip_engine <- function() {
  if ("bartisan" %in% parsnip::get_model_env()[["bart"]][["engine"]]) {
    return(invisible())
  }

  modes <- c("regression", "classification", "censored regression")

  # parsnip defines `bart()` for the first two modes only; the third is the
  # survival families'.
  parsnip::set_model_mode("bart", "censored regression")

  for (mode in modes) {
    parsnip::set_model_engine("bart", mode = mode, eng = "bartisan")
    parsnip::set_dependency("bart", eng = "bartisan", pkg = "bartisan",
                            mode = mode)
  }

  args <- list(trees = list("num_trees", list(pkg = "dials", fun = "trees",
                                               range = c(50L, 500L))),
               prior_terminal_node_coef = list("gamma",
                                               list(pkg = "dials",
                                                    fun = "prior_terminal_node_coef")),
               prior_terminal_node_expo = list("beta",
                                               list(pkg = "dials",
                                                    fun = "prior_terminal_node_expo")),
               # The dials default reaches down to 0, where the prior on the
               # leaves is so wide that a grid point near it is a wasted fit.
               prior_outcome_range = list("k",
                                          list(pkg = "dials",
                                               fun = "prior_outcome_range",
                                               range = c(0.5, 5))))

  for (nm in names(args)) {
    parsnip::set_model_arg(model = "bart", eng = "bartisan", parsnip = nm,
                           original = args[[nm]][[1L]],
                           func = args[[nm]][[2L]],
                           has_submodel = FALSE)
  }

  for (mode in modes) {
    # The formula interface, so that factors arrive as factors and the
    # categorical splitting rules see them.
    parsnip::set_fit(model = "bart", eng = "bartisan", mode = mode,
                     value = list(interface = "formula",
                                  protect = c("formula", "data", "weights"),
                                  func = c(pkg = "bartisan", fun = "bartisan"),
                                  defaults = list()))

    parsnip::set_encoding(model = "bart", eng = "bartisan", mode = mode,
                          options = list(predictor_indicators = "none",
                                         compute_intercept = FALSE,
                                         remove_intercept = FALSE,
                                         allow_sparse_x = FALSE))
  }

  # Every prediction goes through `predict()` or `posterior_predict()`, the only
  # entry points parsnip can reach by name; the reshaping into the columns it
  # expects happens in the `post` functions, which can be the package's own.
  prediction <- function(type, post = NULL, draws = FALSE) {
    list(pre = NULL, post = post,
         func = c(pkg = "stats", fun = "predict"),
         args = list(object = quote(object$fit), newdata = quote(new_data),
                     type = type, draws = draws))
  }

  parsnip::set_pred(model = "bart", eng = "bartisan", mode = "regression",
                    type = "numeric", value = prediction("response"))
  parsnip::set_pred(model = "bart", eng = "bartisan", mode = "regression",
                    type = "conf_int",
                    value = prediction("response", parsnip_interval("conf_int"),
                                       draws = TRUE))
  parsnip::set_pred(model = "bart", eng = "bartisan", mode = "regression",
                    type = "pred_int",
                    value = list(pre = NULL,
                                 post = parsnip_interval("pred_int"),
                                 func = c(pkg = "rstantools",
                                          fun = "posterior_predict"),
                                 args = list(object = quote(object$fit),
                                             newdata = quote(new_data))))
  parsnip::set_pred(model = "bart", eng = "bartisan", mode = "regression",
                    type = "raw", value = prediction("response", draws = TRUE))

  parsnip::set_pred(model = "bart", eng = "bartisan", mode = "classification",
                    type = "class", value = prediction("class"))
  parsnip::set_pred(model = "bart", eng = "bartisan", mode = "classification",
                    type = "prob", value = prediction("prob", parsnip_prob))
  parsnip::set_pred(model = "bart", eng = "bartisan", mode = "classification",
                    type = "conf_int",
                    value = prediction("prob", parsnip_interval("conf_int"),
                                       draws = TRUE))
  parsnip::set_pred(model = "bart", eng = "bartisan", mode = "classification",
                    type = "raw", value = prediction("prob", draws = TRUE))

  censored <- "censored regression"

  parsnip::set_pred(model = "bart", eng = "bartisan", mode = censored,
                    type = "time", value = prediction("response"))
  parsnip::set_pred(model = "bart", eng = "bartisan", mode = censored,
                    type = "linear_pred",
                    value = prediction("link", parsnip_linear_pred))
  parsnip::set_pred(model = "bart", eng = "bartisan", mode = censored,
                    type = "raw", value = prediction("response", draws = TRUE))

  # The evaluation times and the interval are arguments of parsnip's own
  # prediction functions rather than of the call it builds, and the `post`
  # function sees only that call's result. Calling `list()` hands all of them
  # to `parsnip_survival()`, which makes the predictions itself.
  for (type in c("survival", "hazard")) {
    args <- list(object = quote(object$fit), new_data = quote(new_data),
                 eval_time = quote(eval_time))

    if (identical(type, "survival")) {
      args <- c(args, list(interval = quote(interval), level = quote(level)))
    }

    parsnip::set_pred(model = "bart", eng = "bartisan", mode = censored,
                      type = type,
                      value = list(pre = NULL, post = parsnip_survival(type),
                                   func = c(pkg = "base", fun = "list"),
                                   args = args))
  }

  invisible()
}

# `predict(type = "prob")` names its columns for the levels already; they are
# set again from the levels parsnip recorded, which are the ones it will check
# the columns against.
parsnip_prob <- function(x, object) {
  as.data.frame(x) |>
    setNames(object[["lvl"]])
}

# Equal-tailed posterior intervals from a matrix of draws by rows, or, for
# category probabilities, an array of draws by rows by categories. The level is
# not an argument of the prediction call; parsnip leaves it on the object for
# this function to read.
parsnip_interval <- function(type) {
  function(x, object) {
    extras <- object[["spec"]][["method"]][["pred"]][[type]][["extras"]]
    probs <- c(1 - extras[["level"]], 1 + extras[["level"]]) / 2

    draws <- {
      if (length(dim(x)) == 2L) list(x)
      else lapply(seq_len(dim(x)[3L]), function(j) x[, , j])
    }

    # Unsuffixed for a single quantity, and suffixed with the level otherwise,
    # with all the lower limits first, as parsnip's own engines lay them out.
    suffix <- {
      if (length(draws) == 1L) ""
      else paste0("_", object[["lvl"]])
    }

    limits <- lapply(draws, apply, 2L, stats::quantile, probs = probs,
                     names = FALSE)

    out <- c(lapply(limits, function(l) l[1L, ]) |>
               setNames(paste0(".pred_lower", suffix)),
             lapply(limits, function(l) l[2L, ]) |>
               setNames(paste0(".pred_upper", suffix)))

    if (isTRUE(extras[["std_error"]])) {
      out <- c(out,
               lapply(draws, apply, 2L, stats::sd) |>
                 setNames(paste0(".std_error", suffix)))
    }

    as.data.frame(out, check.names = FALSE)
  }
}

# The survival function or the hazard at parsnip's evaluation times, as one
# tibble per row of the new data, which is the layout censored's engines use.
# `x` is the list the registered call builds. The survival function is the
# posterior mean, with equal-tailed limits from the draws when an interval is
# asked for; the hazard is the posterior mean, as `predict()` reports it.
parsnip_survival <- function(type) {
  # Forced, since the registration calls this in a loop over the types.
  force(type)

  function(x, object) {
    times <- x[["eval_time"]]
    n <- nrow(x[["new_data"]])
    interval <- identical(x[["interval"]], "confidence")

    # parsnip accepts an evaluation time of 0, where `predict()` does not. Every
    # unit has survived to it, and the hazard there is not finite for most of
    # the families, so it is reported as missing.
    positive <- times > 0
    at_zero <- {
      if (identical(type, "survival")) 1
      else NA_real_
    }

    summaries <- list(matrix(at_zero, n, length(times)))
    if (interval) {
      summaries <- rep(summaries, 3L)
    }
    names(summaries) <- c(paste0(".pred_", type), ".pred_lower",
                          ".pred_upper")[seq_along(summaries)]

    if (any(positive)) {
      pred <- stats::predict(x[["object"]], newdata = x[["new_data"]],
                             type = type, times = times[positive],
                             draws = interval)

      if (interval) {
        probs <- c(1 - x[["level"]], 1 + x[["level"]]) / 2

        summaries[[1L]][, positive] <- apply(pred, c(2L, 3L), mean)
        summaries[[2L]][, positive] <- apply(pred, c(2L, 3L), stats::quantile,
                                             probs = probs[1L], names = FALSE)
        summaries[[3L]][, positive] <- apply(pred, c(2L, 3L), stats::quantile,
                                             probs = probs[2L], names = FALSE)
      }
      else {
        summaries[[1L]][, positive] <- pred
      }
    }

    rows <- lapply(seq_len(n), function(i) {
      c(list(.eval_time = times), lapply(summaries, function(s) s[i, ])) |>
        as_tbl()
    })

    as_tbl(list(.pred = rows))
  }
}

# The location of the log survival time, which increases with the time for
# every family; for `ph()` that is the negative of its log hazard. A family
# with a second forest has its scale's predictor in the second column, which is
# left out.
parsnip_linear_pred <- function(x, object) {
  if (is.matrix(x)) {
    x <- x[, 1L]
  }

  if (identical(object[["fit"]][["family"]][["family"]], "ph")) {
    x <- -x
  }

  unname(x)
}

#' @rawNamespace S3method(parsnip::predict_linear_pred, "_bartisan_fit")
predict_linear_pred._bartisan_fit <- function(object, new_data, ...,
                                              increasing = TRUE) {
  res <- NextMethod()

  if (!increasing) {
    res <- -res
  }

  res
}

# A tibble, without a dependency on tibble: these predictions go only to
# parsnip, which has it loaded, and the class and the compact row names are all
# that distinguish one from a list of columns.
as_tbl <- function(x) {
  structure(x, class = c("tbl_df", "tbl", "data.frame"),
            row.names = .set_row_names(length(x[[1L]])))
}

# Run `fun` now if `pkg` is loaded, and again whenever it is loaded later, which
# covers both orders in which the two packages can be loaded.
on_package_load <- function(pkg, fun) {
  if (isNamespaceLoaded(pkg)) {
    fun()
  }

  setHook(packageEvent(pkg, "onLoad"), function(...) fun())
}
