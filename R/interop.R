#' Interfaces to other packages
#'
#' @description
#' Methods that let a `<bartisan_fit>` object be used by the packages that assess
#' model fit. There is nothing to set up: load the other package and call its
#' function on the fit.
#'
#' @param object,model,x a `<bartisan_fit>` object; the output of a call to
#'   [bartisan()].
#' @inheritParams predict.bartisan_fit
#' @param type string; for `fitted()`, the prediction scale, passed to
#'   [predict.bartisan_fit()]; default is `"response"`. For `pp_check()`, the
#'   name of the \pkg{bayesplot} check to run without its `ppc_` prefix, so that
#'   `"dens_overlay"` (the default) calls
#'   \pkgfun{bayesplot}{ppc_dens_overlay}; \pkgfun{bayesplot}{available_ppc}
#'   lists them.
#' @param offset,weights an offset and prior weights for `newdata`, as in
#'   [predict.bartisan_fit()]. Defaults are `NULL` to use those the model was fit
#'   with. For a binomial response the weights are the numbers of trials, and so
#'   are what a replicate outcome is a fraction of; they must be given alongside
#'   `newdata` when the model was fit with more than one trial, since the number
#'   of trials is not a function of the predictors and cannot be reconstructed.
#' @param transform `logical`; for `posterior_linpred()`, whether to map the
#'   predictor through the inverse link, which is what `posterior_epred()` does.
#'   Default is `FALSE`.
#' @param nsim,ndraws `numeric`; the number of posterior draws to use, chosen at
#'   random from the retained ones. Defaults are 1 for `simulate()` and 10 for
#'   `pp_check()`. A `ppc_loo_*` check uses every retained draw whatever this is
#'   set to, and says so; see Details.
#' @param seed optional seed, set with [set.seed()] before drawing and restored
#'   afterwards, following the [stats::simulate()] convention. Default is `NULL`
#'   to leave the stream alone.
#' @param metrics `character`; for `model_performance()`, which fit statistics to
#'   report. Allowable options include `"all"` (the default), `"ELPD"`,
#'   `"LOOIC"`, `"WAIC"`, `"R2"`, `"RMSE"`, and `"SIGMA"`, and a vector of them
#'   selects several.
#' @param eta for `as_draws()`, which columns of the additive predictor to carry
#'   into the draws array alongside the scalar parameters, given as either a
#'   logical value or a numeric vector. Default is `TRUE`, which takes a
#'   representative ten spread across the range of the fitted function; `FALSE`
#'   takes none, and a numeric vector takes those observations. The default takes
#'   a handful rather than all of them because there is one column per
#'   observation, and an array with thousands of them is not something
#'   \pkgfun{posterior}{summarise_draws} or a trace plot can be pointed at. The
#'   predictor is the quantity whose convergence usually matters, and the one
#'   [diagnose()] reports on, so it is included by default.
#' @param verbose `logical`; whether to report problems that do not stop the
#'   computation, such as a family that has no mean and so no Bayesian
#'   \eqn{R^2}. Default is `TRUE`.
#' @param ... further arguments, passed to whatever the method calls (the
#'   \pkg{bayesplot} check from `pp_check()`, \pkgfun{loo}{loo} and
#'   \pkgfun{loo}{waic} from `loo()` and `waic()`, and
#'   [predict.bartisan_fit()] from the rest) and ignored where there is nowhere
#'   to pass them.
#'
#' @returns
#' `posterior_predict()`, `posterior_epred()`, `posterior_linpred()` and
#' `log_lik()` return a matrix of draws by observations. `simulate()` returns a
#' data frame of one column per replicate. `loo()` and `waic()` return the
#' `<loo>` and `<waic>` objects those functions produce, and
#' `model_performance()` a one-row data frame of class `<performance_model>`.
#' `as_draws()` returns a `<draws_array>` of iterations by chains by parameters.
#' The accessors return what their names suggest.
#'
#' @details
#' ## What Is Available
#'
#' **Posterior predictions.** \pkgfun{rstantools}{posterior_predict} draws
#' replicate outcomes from the fitted model,
#' \pkgfun{rstantools}{posterior_epred} gives their mean and
#' \pkgfun{rstantools}{posterior_linpred} the additive predictor, following the
#' \pkg{rstantools} conventions that \pkg{brms} and \pkg{rstanarm} follow, and
#' [stats::simulate()] is the same thing in the shape base R expects.
#' **Pointwise likelihood.** \pkgfun{rstantools}{log_lik} returns the
#' draws-by-observations matrix of log-likelihood contributions, which is what
#' \pkgfun{loo}{loo} and \pkgfun{loo}{waic} need; both have methods here.
#'
#' **Graphical checks.** `pp_check()` runs any of the \pkg{bayesplot}
#' posterior-predictive checks on the fit. **Summaries.**
#' \pkgfun{performance}{model_performance} collects the fit statistics in one
#' table, \pkgfun{performance}{r2} gives the Bayesian \eqn{R^2}, and
#' \pkgfun{posterior}{as_draws} hands the scalar parameters to
#' \pkgfun{posterior}{summarise_draws} or to the \pkg{bayesplot} MCMC
#' diagnostics. **Basic accessors.** [stats::fitted()], [stats::residuals()],
#' [stats::weights()] and [stats::sigma()] do what they do for a `glm`, which is
#' also most of what \pkg{insight} needs to make the fit legible to the
#' \pkg{easystats} packages.
#'
#' ## Leave-One-Out Is Approximate, and Mostly Holds Up
#'
#' \pkgfun{loo}{loo} estimates the leave-one-out predictive density by importance
#' sampling from the full-data posterior, and the estimate is trustworthy only
#' when the importance weights have a finite variance, which is what the Pareto
#' \eqn{k} diagnostic reports on. The worry for a forest is that it is a very
#' flexible function of the predictors, so one observation might carry enough
#' influence over the leaves it lands in that dropping it cannot be approximated
#' from the fit in hand.
#'
#' Measured, it usually does not. Over nine fits (`gaussian()` at \eqn{n = 100}
#' with 200 trees, `gaussian_ls()`, `poisson()`, `binomial()`, `dpm()`, and
#' [bcf()] on `lalonde` with both `gaussian()` and `tweedie()`), at most 0.2% of
#' observations exceeded \eqn{k = 0.7} and the median \eqn{k} ran between 0.03
#' and 0.32. The leaf prior is what makes the difference: it shrinks every leaf
#' towards zero and the fit is a sum over many trees, so no single observation
#' dominates the leaves it reaches. The exceptions that did turn up were about
#' the likelihood rather than the trees, and are the ones worth having; fitting
#' \eqn{t_2} errors with `gaussian()` left one observation of 400 at
#' \eqn{k = 2.4}.
#'
#' So the warning \pkg{loo} prints there is worth reading rather than expecting.
#' When it names a handful of observations, those are the influential ones, and
#' refitting without them is what shows how badly they are predicted. A log
#' score on data the model has not seen is available directly:
#' ```r
#' predict(fit, newdata = held_out, type = "density", log = TRUE)
#' ```
#'
#' The seven `ppc_loo_*` checks reweight the replicates towards the
#' leave-one-out predictive instead of comparing them with the response
#' directly, so they need those same weights. `pp_check()` computes them from
#' the fit's own pointwise log likelihood and passes them on, and `ndraws` does
#' not apply to those checks, because the weights and the replicates have to
#' line up draw for draw; supplying `lw` or `psis_object` takes over from it.
#' \pkgfun{bayesplot}{ppc_loo_calibration} wants a binary response besides,
#' which is its own requirement rather than this package's.
#'
#' The two calibration checks are the ones to reach for when the response is
#' binary, since the default check compares two distributions that can only take
#' two values and so finds nothing. `type = "loo_calibration"` is the honest
#' one, holding each observation out of the probability it is judged against;
#' `type = "calibration"` is its in-sample counterpart and reads optimistically.
#' A binned residual plot (`type = "error_binned"`) and either calibration check
#' are about the predicted probabilities rather than replicate outcomes, so they
#' are passed the mean of the predictive distribution instead of a draw from it.
#'
#' ## What a Posterior Predictive Draw Is On
#'
#' The replicate outcomes are on the scale the likelihood was written on, which
#' is the scale [bartisan()] stored the response on. A **binomial** response is a
#' proportion, so binary data come back as 0 and 1, and data given as two columns
#' or with prior weights come back as a fraction of the trials. A response with
#' **categories** comes back as an integer category index, from 1 to the number
#' of categories, because a matrix cannot hold a factor; `fit$levels` names them,
#' and [stats::simulate()] returns factors instead, since its result is a data
#' frame and can. An **accelerated failure time** response comes back as a time
#' rather than a log time, and it is an event time: the predictive distribution
#' of the outcome does not know about the censoring that may have hidden it, so
#' comparing replicates against censored observations is not like for like, and
#' `pp_check()` says so. A [custom_family()] fit has **no posterior predictive
#' distribution at all**, because a log density supplies no way to draw from it,
#' so those methods error.
#'
#' ## What Is Deliberately Absent
#'
#' There is no `logLik()` method, and that is a choice rather than a gap. The
#' generic exists so that [stats::AIC()] and [stats::BIC()] can be computed, and
#' both need a count of parameters, which a forest does not have, since the
#' number of leaves is itself drawn from the posterior. \pkgfun{loo}{loo} and
#' \pkgfun{loo}{waic} are the corresponding quantities for a model like this
#' one, and they are computed from the posterior rather than from a parameter
#' count.
#'
#' For the same reason \pkgfun{performance}{check_normality} and
#' \pkgfun{performance}{check_outliers} do not work: they ask for a
#' likelihood-ratio test and for Cook's distance, neither of which is defined
#' here. \pkgfun{performance}{check_predictions} does work, through
#' [stats::simulate()].
#'
#' ## The Bayesian R-Squared
#'
#' \pkgfun{performance}{r2} returns the quantity of Gelman et al. (2019): per
#' draw, the variance of the fitted means across observations divided by that
#' variance plus the variance of the residuals. Being a per-draw quantity it has a
#' posterior, which is why it is reported with an interval and why it can fall as
#' the model is made more flexible. It needs a mean, so it is not available for
#' `ordinal()` or `multinomial()`.
#'
#' @references
#' Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019). R-squared for
#' Bayesian regression models. *The American Statistician*, 73(3), 307--309.
#'
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. *Statistics and
#' Computing*, 27(5), 1413--1432.
#'
#' @seealso
#' [predict.bartisan_fit()] for the predictions these methods are built on;
#' [diagnose()] for the convergence and mixing diagnostics;
#' [bartisan-marginaleffects] for reading effects off a fit;
#' `vignette("diagnostics")` for the fuller treatment
#'
#' @examplesIf rlang::is_installed(c("loo", "rstantools"))
#' data("rhc")
#' set.seed(123)
#'
#' fit <- bartisan(death ~ . - days, data = rhc, num_trees = 10,
#'                 num_burn = 50, num_draws = 50, chains = 2, verbose = FALSE)
#'
#' # Replicate outcomes, one per draw per observation, whose mean is what
#' # `fitted()` reports
#' yrep <- rstantools::posterior_predict(fit)
#' range(colMeans(rstantools::posterior_epred(fit)) - fitted(fit))
#'
#' # Pointwise log likelihood, and the fit statistics built on it
#' loo::waic(rstantools::log_lik(fit))
#'
#' # Whether replicate outcomes look like the observed ones
#' if (rlang::is_installed("bayesplot")) {
#'   bayesplot::pp_check(fit, type = "bars")
#' }
#'
#' # The scalar parameters and a spread of the predictor, as a draws array
#' if (rlang::is_installed("posterior")) {
#'   posterior::summarise_draws(posterior::as_draws(fit))
#' }
#'
#' @name bartisan-interop
#' @importFrom stats fitted residuals weights sigma simulate
NULL

# ---------------------------------------------------------------------------
# Posterior predictive sampling
# ---------------------------------------------------------------------------

# One replicate outcome per draw per observation. Where a family has a mean, the
# mean comes from response_scale() rather than being derived a second time here:
# that is the function that already knows about a link the engine does not carry
# natively, and a second derivation would be a second place for it to be wrong.
posterior_sample <- function(object, eta, aux, weights = NULL,
                             iterations = NULL) {
  family <- object[["family"]][["family"]]
  link <- object[["family"]][["link"]]

  e <- eta[[1L]]
  ns <- nrow(e)
  n <- ncol(e)

  # How many numbers a replicate matrix holds. Named for what it is rather than
  # `size`, which is also the name of the negative binomial's dispersion
  # argument.
  cells <- ns * n

  if (identical(family, "custom")) {
    arg::err("a {.fn custom_family} fit supplies a log density and no way to
              draw from it, so it has no posterior predictive distribution to
              sample")
  }

  # A nuisance parameter is one number per draw, shared by every observation in
  # that draw.
  spread <- function(name) {
    matrix(aux[, name], nrow = ns, ncol = n)
  }

  uniform <- function() {
    matrix(stats::runif(cells), nrow = ns, ncol = n)
  }

  square <- function(values) {
    matrix(values, nrow = ns, ncol = n)
  }

  if (family %in% c("ordinal", "multinomial", "mnp")) {
    probs <- category_probs(object, eta, aux)
    u <- uniform()
    codes <- matrix(1L, nrow = ns, ncol = n)
    cumulative <- matrix(0, nrow = ns, ncol = n)

    # Inverse transform on the cumulative probabilities. The last category needs
    # no comparison: whatever is left over falls into it, which is also what
    # keeps a rounding error in the probabilities from producing an index out of
    # range.
    for (k in seq_len(dim(probs)[3L] - 1L)) {
      cumulative <- cumulative + probs[, , k]
      codes <- codes + (u > cumulative)
    }

    return(codes)
  }

  if (identical(family, "aft")) {
    sigma <- spread("sigma")

    # log T = eta + sigma * e, with the error distribution the link names. The
    # Weibull case is the smallest extreme value distribution, which is the log
    # of a unit exponential.
    error <- switch(link,
                    weibull = log(square(stats::rexp(cells))),
                    loglogistic = square(stats::rlogis(cells)),
                    square(stats::rnorm(cells)))

    return(exp(e + sigma * error))
  }

  if (family %in% c("zip", "zinb")) {
    mu <- exp(e)

    counts <- {
      if (identical(family, "zip")) square(stats::rpois(cells, mu))
      else square(stats::rnbinom(cells, size = spread("theta"), mu = mu))
    }

    # A structural zero replaces the count draw rather than being added to it.
    counts[uniform() < stats::plogis(eta[[2L]])] <- 0L

    return(counts)
  }

  if (identical(family, "beta")) {
    mu <- stats::plogis(e)
    phi <- spread("phi")

    return(square(stats::rbeta(cells, mu * phi, phi - mu * phi)))
  }

  if (identical(family, "ordbeta")) {
    mu <- stats::plogis(e)
    phi <- spread("phi")

    out <- square(stats::rbeta(cells, mu * phi, phi - mu * phi))

    # The two point masses are taken from the same uniform draw as each other,
    # from the ends of it, so that the three probabilities partition it.
    u <- uniform()
    out[u < 1 - stats::plogis(e - spread("cut1"))] <- 0
    out[u > 1 - stats::plogis(e - spread("cut2"))] <- 1

    return(out)
  }

  # A replicate from a Dirichlet process mixture picks a component -- one of the
  # occupied ones, or a fresh draw from the baseline -- and then a normal from
  # it. The baseline's marginal is a t, which is where the heavy tails a mixture
  # of normals reaches for come from.
  if (identical(family, "dpm")) {
    total <- object[["n"]] + aux[, "alpha"]
    out <- square(0)

    for (s in seq_len(ns)) {
      components <- mixture_at(object, iterations[s])
      probability <- c(components[, "weight"] * object[["n"]],
                       aux[s, "alpha"]) / total[s]
      picked <- sample.int(length(probability), n, replace = TRUE,
                           prob = probability)
      fresh <- picked > nrow(components)

      # A fresh component is drawn from the baseline, which lives on the raw
      # chart, so it comes back `center` below where the reported ones sit.
      offset <- -aux[s, "center"]
      center <- ifelse(fresh, offset,
                       components[pmin(picked, nrow(components)), "mean"])
      spread <- ifelse(fresh, 1, components[pmin(picked, nrow(components)),
                                            "sd"])
      noise <- stats::rnorm(n) * spread

      baseline <- sqrt(object[["family_opts"]][["lambda"]] *
                         (1 + 1 / object[["family_opts"]][["k_0"]])) *
        stats::rt(n, object[["family_opts"]][["nu"]])

      out[s, ] <- e[s, ] + center + ifelse(fresh, baseline, noise)
    }

    return(out)
  }

  mu <- response_scale(object, eta, aux, draws = TRUE)

  switch(family,
    gaussian = square(stats::rnorm(cells, mu, spread("sigma"))),
    gaussian_ls = square(stats::rnorm(cells, mu, exp(eta[[2L]]))),
    poisson = square(stats::rpois(cells, mu)),
    negbin = square(stats::rnbinom(cells, size = spread("theta"), mu = mu)),
    Gamma_ls = {
      shape <- exp(-eta[[2L]])
      square(stats::rgamma(cells, shape = shape, rate = shape / mu))
    },
    Gamma = {
      shape <- spread("shape")
      square(stats::rgamma(cells, shape = shape, rate = shape / mu))
    },
    # Drawn from the representation rather than the density: a Poisson count of
    # gamma claims, which is what the family is, and which gives the point mass
    # at zero for free whenever the count comes out zero.
    tweedie = {
      p <- spread("power")
      ph <- spread("phi")
      lambda <- mu^(2 - p) / (ph * (2 - p))
      shape <- (2 - p) / (p - 1)
      scale <- ph * (p - 1) * mu^(p - 1)
      # `spread()` gives one value per cell, so every one of these is a matrix
      # of the same shape and the subsetting has to reach all of them.
      n_claims <- stats::rpois(cells, lambda)
      out <- numeric(cells)
      hit <- n_claims > 0L
      out[hit] <- stats::rgamma(sum(hit), shape = n_claims[hit] * shape[hit],
                                scale = scale[hit])
      square(out)
    },
    binomial = {
      trials <- matrix(weights, nrow = ns, ncol = n, byrow = TRUE)

      if (any(abs(trials - round(trials)) > 1e-8)) {
        arg::err("a replicate binomial outcome needs whole numbers of trials,
                  and the prior weights are not whole numbers")
      }

      square(stats::rbinom(cells, size = round(trials), prob = mu)) / trials
    },
    arg::err("the {.val {family}} family has no posterior predictive sampler"))
}

# The response as the likelihood saw it, which is the scale a replicate outcome
# comes back on and so the thing a replicate is comparable to.
observed_response <- function(object) {
  family <- object[["family"]][["family"]]
  y <- object[["y"]]

  if (family %in% c("ordinal", "multinomial", "mnp")) {
    # Stored zero-based, because that is what the sampler indexes categories by.
    return(as.integer(y) + 1L)
  }

  if (identical(family, "aft")) {
    return(exp(y))
  }

  y
}

#' @rdname bartisan-interop
#' @exportS3Method rstantools::posterior_predict
posterior_predict.bartisan_fit <- function(object, newdata = NULL, iterations = NULL,
                                      offset = NULL, weights = NULL, ...) {

  parts <- predict_parts(object, newdata, offset, iterations)

  # For a binomial response the prior weights are the numbers of trials, and a
  # replicate outcome is a fraction of them -- so they are not a function of the
  # predictors and cannot be reconstructed for new data, the same situation an
  # offset is in. Defaulting silently to one trial would return a plausible 0/1
  # answer to a question about counts, so it is an error instead.
  if (is_null(weights) && !is_null(newdata) &&
        identical(object[["family"]][["family"]], "binomial") &&
        any(object[["prior_weights"]] != 1)) {
    arg::err("the model was fit to a binomial response with more than one trial,
              so {.arg weights} must give the number of trials for each row of
              {.arg newdata}")
  }

  trials <- {
    if (!is_null(weights)) weights
    else if (is_null(newdata)) object[["prior_weights"]]
    else rep.int(1, ncol(parts[["eta"]][[1L]]))
  }

  posterior_sample(object, parts[["eta"]], parts[["aux"]], trials,
                   parts[["iterations"]])
}

#' @rdname bartisan-interop
#' @exportS3Method rstantools::posterior_epred
posterior_epred.bartisan_fit <- function(object, newdata = NULL, ...) {
  stats::predict(object, newdata = newdata, type = "response", draws = TRUE,
                 ...)
}

#' @rdname bartisan-interop
#' @exportS3Method rstantools::posterior_linpred
posterior_linpred.bartisan_fit <- function(object, transform = FALSE, newdata = NULL,
                                      ...) {
  arg::arg_flag(transform)

  type <- if (transform) "response" else "link"

  stats::predict(object, newdata = newdata, type = type, draws = TRUE, ...)
}

#' @rdname bartisan-interop
#' @exportS3Method rstantools::log_lik
log_lik.bartisan_fit <- function(object, newdata = NULL, ...) {
  stats::predict(object, newdata = newdata, type = "density", draws = TRUE,
                 log = TRUE, ...)
}

#' @rdname bartisan-interop
#' @export
simulate.bartisan_fit <- function(object, nsim = 1, seed = NULL, ...) {
  arg::arg_count(nsim)
  arg::arg_gte(nsim, 1)

  # The stats::simulate() contract: `seed` is set for the duration and the
  # stream is left as it was found, and the seed used is recorded on the result.
  if (!is_null(seed)) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1L)
    }
    saved <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }

  num_draws <- nrow(object[["sigma_mu"]])
  iterations <- sample.int(num_draws, size = nsim,
                           replace = nsim > num_draws)

  draws <- posterior_predict.bartisan_fit(object, iterations = iterations, ...)

  family <- object[["family"]][["family"]]
  levels <- object[["levels"]]

  # A replicate comes back as a factor whenever the response was one, which is
  # what stats::simulate() does for a glm and what makes the result a drop-in for
  # the observed column. A binomial response given as proportions or with trial
  # counts is not a factor, so it stays numeric.
  binary_factor <- identical(family, "binomial") && !is_null(levels) &&
    all(object[["prior_weights"]] == 1)

  as_factor <- {
    if (family %in% c("ordinal", "multinomial", "mnp")) function(v) {
      factor(levels[v], levels = levels, ordered = identical(family, "ordinal"))
    }
    else if (binary_factor) function(v) {
      factor(levels[v + 1L], levels = levels)
    }
    else identity
  }

  out <- lapply(seq_len(nrow(draws)), function(s) as_factor(draws[s, ])) |>
    list2DF() |>
    setNames(sprintf("sim_%d", seq_len(nsim)))

  attr(out, "seed") <- seed

  out
}

# ---------------------------------------------------------------------------
# Basic accessors
# ---------------------------------------------------------------------------

#' @rdname bartisan-interop
#' @export
fitted.bartisan_fit <- function(object, type = "response", ...) {
  stats::predict(object, type = type, ...)
}

#' @rdname bartisan-interop
#' @export
residuals.bartisan_fit <- function(object, ...) {
  family <- object[["family"]][["family"]]

  if (family %in% c("ordinal", "multinomial", "mnp")) {
    arg::err(c("A response with more than two categories has no mean, so it has
                no residual on the response scale.",
               i = "Use {.code predict(object, type = \"prob\")}, or
                    {.code type = \"mean\"} if the category labels are numbers."))
  }

  if (identical(family, "custom")) {
    arg::err("a {.fn custom_family} fit supplies a log density and no mean, so
              there is no residual to take")
  }

  # The predictor rather than the median survival time: an accelerated failure
  # time model is a linear model for the log time, and the residual that means
  # anything there is on that scale.
  if (identical(family, "aft")) {
    return(object[["y"]] - colMeans(object[["eta"]][[1L]]))
  }

  observed_response(object) - stats::predict(object, type = "response")
}

#' @rdname bartisan-interop
#' @export
weights.bartisan_fit <- function(object, ...) {
  object[["prior_weights"]]
}

#' @rdname bartisan-interop
#' @export
sigma.bartisan_fit <- function(object, ...) {
  aux <- object[["aux"]]

  # A Dirichlet process mixture has no single scale parameter, but its error
  # distribution has a standard deviation, and that is the comparable number.
  if (identical(object[["family"]][["family"]], "dpm")) {
    return(mean(aux[, "error_sd"]))
  }

  # Present for the families whose scale is a single number: the Gaussian
  # residual standard deviation and the accelerated failure time scale. A
  # location-scale fit has a scale per observation, and it is the second
  # predictor rather than a nuisance parameter.
  if (is_null(aux) || !any(colnames(aux) == "sigma")) {
    return(NULL)
  }

  mean(aux[, "sigma"])
}

# ---------------------------------------------------------------------------
# loo
# ---------------------------------------------------------------------------

# Which chain each stored draw came from. The chains are stacked one after
# another, so this is the block structure of that stacking; loo needs it to
# estimate the efficiency of the draws, and treating dependent draws as
# independent would understate the standard errors it reports.
chain_ids <- function(object) {
  chains <- object[["chains"]]
  per <- nrow(object[["sigma_mu"]]) / chains

  rep(seq_len(chains), each = per)
}

#' @rdname bartisan-interop
#' @exportS3Method loo::loo
loo.bartisan_fit <- function(x, ...) {
  ll <- log_lik.bartisan_fit(x)
  r_eff <- loo::relative_eff(exp(ll), chain_id = chain_ids(x))

  loo::loo.matrix(ll, r_eff = r_eff, ...)
}

#' @rdname bartisan-interop
#' @exportS3Method loo::waic
waic.bartisan_fit <- function(x, ...) {
  log_lik.bartisan_fit(x) |>
    loo::waic.matrix(...)
}

# ---------------------------------------------------------------------------
# bayesplot
# ---------------------------------------------------------------------------

# The Pareto-smoothed importance weights that turn posterior predictive draws
# into leave-one-out predictive draws, which is what the `ppc_loo_*` checks
# reweight by. `loo()` computes these on its way to an ELPD; only the weights
# are wanted here, so `psis()` is called directly.
loo_weights <- function(object) {
  if (!rlang::is_installed("loo")) {
    arg::err(c("A leave-one-out check needs the {.pkg loo} package.",
               i = "Install it, or use a check that compares the replicates
                    against the response directly."))
  }

  ll <- log_lik.bartisan_fit(object)
  r_eff <- loo::relative_eff(exp(ll), chain_id = chain_ids(object))

  # The weights are the reciprocal of the density of each observation, so the
  # log ratios are the negated log likelihood. `loo` warns about the Pareto
  # diagnostic itself, in wording that is right for what it is and says nothing
  # about which observations or what to do; the warning below replaces it.
  psis <- withCallingHandlers(
    loo::psis(-ll, r_eff = r_eff),
    warning = function(w) {
      if (grepl("Pareto", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })

  bad <- which(loo::pareto_k_values(psis) > 0.7)

  if (!is_null(bad)) {
    arg::wrn(c("The leave-one-out weights did not converge for
                {length(bad)} observation{?s}, at {?index/indices}
                {.val {utils::head(bad, 5L)}}.",
               i = "Those observations are influential enough that dropping
                    them cannot be approximated from this fit, so the check
                    understates how badly they are predicted. Refitting
                    without them is the way to see it."))
  }

  psis
}

# Which observations are events, for the Kaplan-Meier checks. The response
# reached `bartisan()` as a survival object and the model frame still holds it,
# so the check does not have to be told something the fit already knows.
survival_status <- function(object) {
  y <- stats::model.response(object[["model"]])

  if (!inherits(y, "Surv")) {
    arg::err(c("A Kaplan-Meier check needs to know which observations are
                events, and this fit's response is not a survival object.",
               i = "Pass {.arg status_y} to say which are."))
  }

  as.numeric(y[, "status"])
}

#' @rdname bartisan-interop
#' @exportS3Method bayesplot::pp_check
pp_check.bartisan_fit <- function(object, type = "dens_overlay", ndraws = 10, ...) {
  arg::arg_string(type)
  arg::arg_count(ndraws)
  arg::arg_gte(ndraws, 1)

  fun <- sprintf("ppc_%s", type)

  if (!fun %in% getNamespaceExports("bayesplot")) {
    arg::err("{.val {type}} is not a {.pkg bayesplot} posterior predictive
              check; {.fn bayesplot::available_ppc} lists them, without the
              {.val ppc_} prefix")
  }

  if (identical(object[["family"]][["family"]], "aft")) {
    arg::wrn(c("The replicates are event times and the observed times may
                be censored, so the comparison is not like for like.",
               i = "Restrict to the uncensored observations, or check the
                    predictor instead."))
  }

  dots <- list(...)

  # A `ppc_loo_*` check does not compare the replicates with the response
  # directly: it reweights them towards the leave-one-out predictive first, and
  # so needs the importance weights as well. It also needs every draw of them,
  # because *bayesplot* requires the weights and the replicates to be the same
  # shape -- which is why `ndraws` cannot apply here.
  loo_check <- startsWith(type, "loo_")
  supplied <- any(c("lw", "psis_object") %in% names(dots))

  if (loo_check && !supplied) {
    if (!missing(ndraws)) {
      arg::wrn(c("{.arg ndraws} does not apply to a leave-one-out check, and
                  every retained draw is used.",
                 i = "The weights are computed per draw and have to line up
                      with the replicates draw for draw."))
    }

    dots[["psis_object"]] <- loo_weights(object)
  }

  # `ppc_km_overlay` overlays the replicate survival curves on the observed
  # Kaplan-Meier curve, which takes the censoring indicator as well as the
  # times. It is the one check written for a censored response, so the fit
  # supplies it rather than making the caller repeat it.
  if (startsWith(type, "km_overlay") && !"status_y" %in% names(dots)) {
    dots[["status_y"]] <- survival_status(object)
  }

  num_draws <- nrow(object[["sigma_mu"]])

  # All of them for a leave-one-out check, and for one whose weights the caller
  # computed themselves, since those came from every draw too.
  iterations <- if (loo_check || supplied) {
    NULL
  }
  else {
    sample.int(num_draws, size = min(ndraws, num_draws))
  }

  # A binned residual plot and a calibration plot are about the predicted
  # probabilities, not about replicate outcomes: both bin the second argument
  # and read the outcome within each bin, which a vector of zeros and ones
  # gives two degenerate bins of. So those checks get the mean of the
  # predictive distribution rather than a draw from it, which is what rstanarm
  # passes for the same two. `ppc_loo_calibration()` is not among them: it
  # takes replicates and forms the leave-one-out probabilities itself.
  epred_check <- type %in% c("error_binned", "calibration",
                             "calibration_grouped", "calibration_overlay",
                             "calibration_overlay_grouped")

  reps <- if (epred_check) {
    epred <- posterior_epred.bartisan_fit(object)

    if (is.null(iterations)) epred else epred[iterations, , drop = FALSE]
  }
  else {
    posterior_predict.bartisan_fit(object, iterations = iterations)
  }

  # `ppc_calibration()` names that argument `prep`, and takes `yrep` as an
  # alternative it converts; the two differ in position between the members of
  # its own family, so it is named rather than passed along.
  args <- if (startsWith(type, "calibration")) {
    c(list(observed_response(object), prep = reps), dots)
  }
  else {
    c(list(observed_response(object), reps), dots)
  }

  do.call(getExportedValue("bayesplot", fun), args)
}

# ---------------------------------------------------------------------------
# posterior
# ---------------------------------------------------------------------------

#' @rdname bartisan-interop
#' @exportS3Method posterior::as_draws
as_draws.bartisan_fit <- function(x, eta = TRUE, ...) {
  scalars <- scalar_draws(x)

  # The additive predictor is the quantity whose convergence actually matters --
  # `diagnose()` already reports it -- so a handful of its columns belong here too,
  # or the diagnostics that read this object can only see the nuisance
  # parameters. A handful rather than all of them: there is one per observation,
  # and a `draws_array` with thousands of columns is not something
  # `summarise_draws()` or a trace plot can be pointed at.
  if (!isFALSE(eta)) {
    scalars <- c(scalars, eta_draws(x, eta))
  }

  chains <- x[["chains"]]
  per <- nrow(x[["sigma_mu"]]) / chains

  array(unlist(scalars, use.names = FALSE),
        dim = c(per, chains, length(scalars)),
        dimnames = list(NULL, NULL, names(scalars))) |>
    posterior::as_draws_array()
}

# The columns of `eta` to carry into a draws array. `which` is TRUE for a
# representative spread, or observation indices to take exactly those.
#
# The spread is taken over the posterior mean of the predictor rather than at
# random, so the selection covers the range of the fitted function: the
# observations that mix worst are usually the ones at its edges, where the fewest
# observations inform the leaves.
eta_draws <- function(object, which = TRUE, size = 10L) {
  out <- list()

  for (h in seq_along(object[["eta"]])) {
    draws <- object[["eta"]][[h]]
    n <- ncol(draws)

    index <- {
      if (isTRUE(which)) {
        if (n <= size) seq_len(n)
        else {
          ordered <- order(colMeans(draws))
          ordered[unique(round(seq(1, n, length.out = size)))]
        }
      }
      else {
        arg::arg_numeric(which)

        if (any(which < 1) || any(which > n)) {
          arg::err("{.arg eta} must be observation indices between 1 and {n}")
        }

        as.integer(which)
      }
    }

    # With one forest there is nothing to disambiguate, so `eta[3]` rather than
    # the `eta.eta[3]` the multi-predictor families need.
    label <- {
      if (length(object[["eta"]]) == 1L) "eta"
      else sprintf("eta.%s", names(object[["eta"]])[h])
    }

    for (i in index) {
      out[[sprintf("%s[%d]", label, i)]] <- draws[, i]
    }
  }

  out
}

# ---------------------------------------------------------------------------
# performance
# ---------------------------------------------------------------------------

#' @rdname bartisan-interop
#' @exportS3Method performance::r2_posterior
r2_posterior.bartisan_fit <- function(model, verbose = TRUE, ...) {
  family <- model[["family"]][["family"]]

  if (family %in% c("ordinal", "multinomial", "mnp", "custom")) {
    if (verbose) {
      arg::wrn("the {.val {family}} family has no mean, so it has no
                Bayesian {.field R2}")
    }
    return(NULL)
  }

  mu <- posterior_epred.bartisan_fit(model)
  y <- observed_response(model)

  # Gelman et al. (2019): both variances are taken across observations within a
  # draw, so the ratio has a posterior of its own.
  fit_var <- apply(mu, 1L, stats::var)
  residual_var <- apply(mu, 1L, function(m) stats::var(y - m))

  list(R2_Bayes = fit_var / (fit_var + residual_var))
}

#' @rdname bartisan-interop
#' @exportS3Method performance::r2
r2.bartisan_fit <- function(model, ...) {
  performance::r2_bayes(model, ...)
}

#' @rdname bartisan-interop
#' @exportS3Method performance::model_performance
model_performance.bartisan_fit <- function(model, metrics = "all", verbose = TRUE,
                                      ...) {

  all_metrics <- c("ELPD", "LOOIC", "WAIC", "R2", "RMSE", "SIGMA")

  if (identical(metrics, "all")) {
    metrics <- all_metrics
  }
  else {
    metrics <- intersect(all_metrics, toupper(as.character(metrics)))
  }

  out <- list()

  # Every one of these three comes out of the pointwise likelihood through loo,
  # and model_performance() is reachable from performance alone, so this is the
  # one place in the file where a suggested package has to be asked for rather
  # than assumed: the other methods are only callable through their own
  # package's generic.
  if (any(c("ELPD", "LOOIC", "WAIC") %in% metrics)) {
    rlang::check_installed("loo", "to report the ELPD, LOOIC or WAIC.")
  }

  if (any(c("ELPD", "LOOIC") %in% metrics)) {
    estimates <- suppressWarnings(loo.bartisan_fit(model))[["estimates"]]

    if ("ELPD" %in% metrics) {
      out[["ELPD"]] <- estimates["elpd_loo", "Estimate"]
      out[["ELPD_SE"]] <- estimates["elpd_loo", "SE"]
    }

    if ("LOOIC" %in% metrics) {
      out[["LOOIC"]] <- estimates["looic", "Estimate"]
      out[["LOOIC_SE"]] <- estimates["looic", "SE"]
    }
  }

  if ("WAIC" %in% metrics) {
    estimates <- suppressWarnings(waic.bartisan_fit(model))[["estimates"]]
    out[["WAIC"]] <- estimates["waic", "Estimate"]
  }

  if ("R2" %in% metrics) {
    posterior_r2 <- r2_posterior.bartisan_fit(model, verbose = FALSE)

    if (!is_null(posterior_r2)) {
      out[["R2"]] <- mean(posterior_r2[["R2_Bayes"]])
    }
  }

  if (any(c("RMSE", "SIGMA") %in% metrics)) {
    residual <- {
      if (model[["family"]][["family"]] %in%
          c("ordinal", "multinomial", "mnp", "custom")) NULL
      else stats::residuals(model)
    }

    if ("RMSE" %in% metrics && !is_null(residual)) {
      out[["RMSE"]] <- sqrt(mean(residual^2))
    }

    if ("SIGMA" %in% metrics) {
      out[["Sigma"]] <- stats::sigma(model)
    }
  }

  out <- list2DF(out)

  class(out) <- c("performance_model", class(out))

  out
}
