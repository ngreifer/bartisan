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
#' @param K `numeric`; for `kfold()`, how many folds to split the sample into.
#'   Default is 10. Ignored when `folds` is given.
#' @param folds optional; for `kfold()`, an integer vector of one fold number
#'   per observation, as \pkgfun{loo}{kfold_split_random} and its relatives
#'   return. Default is `NULL` to draw them at random. Supply them to stratify,
#'   to group, or to score two models on the same split.
#' @param save_fits `logical`; for `kfold()`, whether to keep the \eqn{K} refits
#'   in the result's `fits` element. Default is `FALSE`, since each is a whole
#'   fit.
#' @param scale `string`; for `loo()`, `waic()` and `kfold()` on a survival fit,
#'   the
#'   measure to report the pointwise densities with respect to: `"time"` for the
#'   density of \eqn{T} and `"log_time"` for the density of \eqn{\log T}.
#'   Default is `NULL` to use the family's own, which is \eqn{\log T} for the
#'   accelerated failure time families and \eqn{T} for [ph()]. A fit already on
#'   the scale named is left alone, so naming one scale for every model in a
#'   comparison is enough. See Details.
#' @param metrics `character`; for `model_performance()`, which fit statistics to
#'   report. Allowable options include `"all"` (the default), `"ELPD"`,
#'   `"LOOIC"`, `"WAIC"`, `"R2"`, `"RMSE"`, and `"SIGMA"`, and a vector of them
#'   selects several.
#' @param digits `integer`; for `print()` on the output of `prior_summary()`,
#'   how many digits to round the prior's scales to. Default is 3.
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
#' `kfold()` returns a `<kfold>` object, a list whose `estimates` holds
#' `elpd_kfold`, `p_kfold` and `kfoldic` with their standard errors, whose
#' `pointwise` holds the same three per observation, and whose `folds` records
#' the split; `save_fits = TRUE` adds the \eqn{K} refits in `fits`.
#'
#' `posterior_predict()`, `posterior_epred()`, `posterior_linpred()` and
#' `log_lik()` return a matrix of draws by observations. `simulate()` returns a
#' data frame of one column per replicate. `loo()` and `waic()` return the
#' `<loo>` and `<waic>` objects those functions produce, and
#' `model_performance()` a one-row data frame of class `<performance_model>`.
#' `as_draws()` returns a `<draws_array>` of iterations by chains by parameters.
#'
#' `prior_summary()` returns a `<bartisan_prior_summary>` object, a list whose
#' `forests` is a data frame of one row per additive predictor and one column per
#' setting the prior is made of, whose `estimated` says in the same shape which
#' of them were drawn rather than held, and whose `family` holds the family's own
#' parameters with the prior each was given. `random` and `response` record the
#' group-intercept scale and what was read off the response.
#'
#' The accessors return what their names suggest.
#'
#' @details
#' ## Available Methods
#'
#' \pkgfun{rstantools}{posterior_predict} draws replicate outcomes from the
#' fitted model, \pkgfun{rstantools}{posterior_epred} gives their mean and
#' \pkgfun{rstantools}{posterior_linpred} the additive predictor, following the
#' \pkg{rstantools} conventions that \pkg{brms} and \pkg{rstanarm} follow, and
#' [stats::simulate()] is the same thing in the shape base R expects.
#' \pkgfun{rstantools}{log_lik} returns the draws-by-observations matrix of
#' log-likelihood contributions, which is what \pkgfun{loo}{loo} and
#' \pkgfun{loo}{waic} need.
#'
#' `pp_check()` runs any of the \pkg{bayesplot} posterior-predictive checks on
#' the fit, \pkgfun{performance}{model_performance} collects the fit statistics
#' in one table, \pkgfun{performance}{r2} gives the Bayesian \eqn{R^2}, and
#' \pkgfun{posterior}{as_draws} hands the scalar parameters to
#' \pkgfun{posterior}{summarise_draws} or to the \pkg{bayesplot} MCMC
#' diagnostics. \pkgfun{rstantools}{prior_summary} writes out every prior the fit
#' was given, on the scale it was given on, which is the companion to
#' `prior_only = TRUE` in [bartisan()]: one says what the prior is and the other
#' what it implies about the outcome. [stats::fitted()], [stats::residuals()],
#' [stats::weights()] and [stats::sigma()] do what they do for a `glm`, which is
#' also most of what \pkg{insight} needs to make the fit legible to the
#' \pkg{easystats} packages.
#'
#' ## Accuracy of the Leave-One-Out Approximation
#'
#' \pkgfun{loo}{loo} estimates the leave-one-out predictive density by importance
#' sampling from the full-data posterior, and the estimate is trustworthy only
#' when the importance weights have a finite variance, which is what the Pareto
#' \eqn{k} diagnostic reports on. A forest is a flexible function of the
#' predictors, so the worry is that one observation carries enough influence over
#' the leaves it lands in that dropping it cannot be approximated from the fit in
#' hand. In practice it rarely does: the leaf prior shrinks every leaf toward zero
#' and the fit is a sum over many trees, so no single observation dominates.
#'
#' The warning \pkg{loo} prints there is therefore worth reading rather than
#' expecting, and the exceptions that do turn up are usually about the likelihood
#' rather than the trees. When it names a handful of observations, those are the
#' influential ones, and refitting without them shows how badly they are
#' predicted. A log score on data the model has not seen is available directly:
#'
#' ```r
#' predict(fit, newdata = held_out, type = "density", log = TRUE)
#' ```
#'
#' ## Cross-Validation (`kfold()`)
#'
#' Where `loo()` estimates the leave-one-out density from one fit, `kfold()`
#' splits the sample, refits \eqn{K} times, and scores each part under a fit that
#' never saw it. That costs \eqn{K} fits and owes nothing to an approximation,
#' which makes it the thing to reach for when the Pareto diagnostics say the
#' weights cannot be trusted.
#'
#' It returns a `<kfold>` object that \pkgfun{loo}{loo_compare} accepts beside a
#' `<loo>` one, so two models can be compared on one split by passing the folds
#' from the first to the second:
#'
#' ```r
#' folds <- loo::kfold_split_random(K = 10, N = nobs(fit))
#'
#' loo_compare(list(full = kfold(fit, folds = folds),
#'                  small = kfold(other, folds = folds)))
#' ```
#'
#' The refits run under a `future` plan when one is set, and one `set.seed()`
#' reproduces them either way. Each is refitted from the original call, so the
#' `data` argument has to still name the data the fit was made from. Prior weights
#' and an offset are carried into both the refits and the held-out scores.
#' `p_kfold` is the gap between what the model predicts for an observation it was
#' fitted to and what it predicts for the same one held out, which is the price of
#' having used it. `vignette("comparison")` reads an example.
#'
#' ## Setting `scale` for Survival Families
#'
#' The accelerated failure time families report the density of \eqn{\log T} and
#' [ph()] the density of \eqn{T}. Both are correct for the model that produced
#' them, and they differ by the Jacobian of the change of variable, so a log score
#' taken across that boundary is off by \eqn{\sum \log t} over the events, which
#' can reverse which family looks better.
#'
#' `scale` puts them on one measure, and reads the same from either side, since a
#' fit already on the scale named is returned untouched:
#'
#' ```r
#' loo_compare(list(aft = loo(aft_fit, scale = "time"),
#'                  ph = loo(ph_fit, scale = "time")))
#' ```
#'
#' It is left to the caller rather than applied automatically so that `loo()`
#' keeps reporting the model's own predictive density and keeps agreeing with
#' `log_lik()`. Censored observations are not adjusted, a survival probability
#' being a probability on either scale. `vignette("comparison")` works through the
#' comparison and `vignette("survival")` through the families.
#'
#' ## Posterior Predictive Checks
#'
#' The seven `ppc_loo_*` checks reweight the replicates towards the
#' leave-one-out predictive, so they need those weights. `pp_check()` computes
#' them from the fit's own pointwise log likelihood and passes them on, and
#' `ndraws` does not apply to those checks, because the weights and the
#' replicates have to line up draw for draw; supplying `lw` or `psis_object`
#' takes over from it.
#'
#' The two calibration checks are the ones to reach for when the response is
#' binary, since the default check compares two distributions that can only take
#' two values. `type = "loo_calibration"` is the honest one, holding each
#' observation out of the probability it is judged against, where
#' `type = "calibration"` is its in-sample counterpart and reads optimistically.
#' Those two and a binned residual plot (`type = "error_binned"`) are about the
#' predicted probabilities rather than replicate outcomes, so they are passed the
#' mean of the predictive distribution instead of a draw from it.
#'
#' ## The Scale of a Posterior Predictive Draw
#'
#' Replicate outcomes are on the scale the likelihood was written on, which is
#' the scale [bartisan()] stored the response on. A binomial response is a
#' proportion, so binary data come back as 0 and 1, and data given as two columns
#' or with prior weights come back as a fraction of the trials. A response with
#' categories comes back as an integer category index, since a matrix cannot hold
#' a factor; `fit$levels` names them, and [stats::simulate()] returns factors
#' instead, its result being a data frame. An accelerated failure time response
#' comes back as an event time rather than a log time, and the predictive
#' distribution of the outcome knows nothing of the censoring that may have hidden
#' it, so `pp_check()` says so when replicates are compared against censored
#' observations. A [custom_family()] fit supplies a log density and no way to draw
#' from it, so these methods error on one.
#'
#' ## AIC, BIC and Normality Checks
#'
#' [stats::AIC()] and [stats::BIC()] need a count of parameters, which a forest
#' does not have, since the number of leaves is itself drawn from the posterior;
#' there is accordingly no `logLik()` method. \pkgfun{loo}{loo} and
#' \pkgfun{loo}{waic} are the corresponding quantities for a model like this one,
#' computed from the posterior rather than from a parameter count. For the same
#' reason \pkgfun{performance}{check_normality} and
#' \pkgfun{performance}{check_outliers}, which ask for a likelihood-ratio test
#' and for Cook's distance, are unavailable, where
#' \pkgfun{performance}{check_predictions} works through [stats::simulate()].
#'
#' ## The Bayesian R-Squared
#'
#' \pkgfun{performance}{r2} returns, per draw, the variance of the fitted means
#' across observations divided by that variance plus the variance of the
#' residuals. Being a per-draw quantity it has a posterior, which is why it is
#' reported with an interval and why it can fall as the model is made more
#' flexible. It needs a mean, so it is available for every family except
#' `ordinal()` and `multinomial()`.
#'
#' @seealso
#' [predict.bartisan_fit()] for the predictions these methods are built on;
#' [diagnose()] for the convergence and mixing diagnostics;
#' [`bartisan-marginaleffects`] for reading effects off a fit;
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
#' # Every prior the fit was given, on the scale it was given on
#' rstantools::prior_summary(fit)
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

  # The response scale of the other two survival families is a median survival
  # time, and the observed time is not a draw from it: a censored time is a
  # bound, not a value. `dpm_aft()` also stores the log time, so the difference
  # was being taken across two scales.
  if (family %in% c("dpm_aft", "ph")) {
    arg::err(c("A censored response has no residual on the response scale: a
                censored time is a lower bound rather than a value, and the
                fitted median survival time is not what it is a draw from.",
               i = "Compare the fit with {.code predict(object, type = \"survival\")}
                    against a Kaplan-Meier estimate, or with
                    {.code bayesplot::pp_check(object, type = \"km_overlay\")}."))
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
# prior_summary
# ---------------------------------------------------------------------------

#' @rdname bartisan-interop
#' @exportS3Method rstantools::prior_summary
prior_summary.bartisan_fit <- function(object, ...) {
  arg::arg_is(object, "bartisan_fit")

  prior <- object[["prior"]]

  if (is_null(prior)) {
    arg::err(c("this fit kept no record of the prior it was given",
               i = "It was made before {.fn prior_summary} existed; refit it to
                    summarize the prior."))
  }

  forests <- predictor_names(object)

  # One row per forest, because nearly every setting is per-forest even when the
  # caller wrote a scalar, and a forest that was given its own is exactly the
  # case worth being able to see.
  table <- data.frame(forest = forests,
                      num_trees = object[["num_trees"]][seq_along(forests)],
                      gamma = prior[["gamma"]],
                      beta = prior[["beta"]],
                      leaf_scale = prior[["sigma_mu"]],
                      k = prior[["k"]],
                      candidates = prior[["candidates"]],
                      alpha = prior[["alpha"]],
                      alpha_scale = prior[["alpha_scale"]],
                      shape_1 = prior[["alpha_shape_1"]],
                      shape_2 = prior[["alpha_shape_2"]],
                      bandwidth = prior[["bandwidth"]],
                      row.names = NULL)

  estimated <- data.frame(forest = forests,
                          leaf_scale = prior[["update_sigma_mu"]],
                          splits = prior[["update_s"]],
                          alpha = prior[["update_alpha"]],
                          bandwidth = prior[["update_bandwidth"]],
                          row.names = NULL)

  out <- list(forests = table,
              estimated = estimated,
              rules = list(soft = isTRUE(object[["soft"]]),
                           gate = object[["gate"]]),
              sparsity = prior[["sparsity"]],
              split_prior = isTRUE(prior[["split_prior"]]),
              share_sparsity = isTRUE(prior[["share_sparsity"]]),
              random = random_prior(object, prior),
              family = list(name = object[["family"]][["family"]],
                            link = object[["family"]][["link"]],
                            parameters = family_prior(object)),
              response = list(intercept = object[["intercept"]][[1L]],
                              scale = prior[["eta_scale"]]),
              prior_only = isTRUE(object[["prior_only"]]))

  class(out) <- "bartisan_prior_summary"

  out
}

# The group intercepts, when the formula has a bar. Each is Gaussian with a
# scale shared by the term's levels, and that scale is given the same
# half-Cauchy the leaf scale gets, at the first forest's leaf scale; see
# `make_random_effects()` in `src/random.cpp`.
random_prior <- function(object, prior) {
  random <- object[["random"]]

  if (is_null(random)) {
    return(NULL)
  }

  list(terms = pluck(random, "label", character(1L)),
       levels = pluck(random, "num_levels", integer(1L)),
       scale = prior[["sigma_mu"]][[1L]],
       estimated = isTRUE(prior[["update_tau"]]))
}

# The family's own parameters, as one row each of what they are and what they
# are given.
#
# Most of them follow one of two naming conventions in `family_opts`, either
# `<p>_prior_shape` and `<p>_prior_rate` or `<p>_shape` and `<p>_rate`, and both
# mean a gamma prior on `<p>` drawn when `update_<p>` is set. The rest are named
# one at a time below. The `_prior_` pair is matched first and its names are
# then withheld from the plain pattern, since `shape_prior_shape` would
# otherwise also read as a gamma prior on something called `shape_prior`.
family_prior <- function(object) {
  opts <- object[["family_opts"]]
  family <- object[["family"]][["family"]]

  rows <- list()

  add <- function(parameter, prior, estimated) {
    rows[[length(rows) + 1L]] <<- data.frame(parameter = parameter,
                                             prior = prior,
                                             estimated = estimated,
                                             row.names = NULL)
  }

  value <- function(nm) {
    signif(as.numeric(opts[[nm]])[[1L]], 4L)
  }

  drawn <- function(nm) {
    flag <- opts[[paste0("update_", nm)]]

    is_null(flag) || isTRUE(as.logical(flag)[[1L]])
  }

  # The Gaussian and accelerated failure time scales, both half-Cauchy at a
  # scale read off the response.
  if (!is_null(opts[["sigma_hat"]])) {
    add("sigma", sprintf("HalfCauchy(0, %s)", value("sigma_hat")),
        drawn("sigma"))
  }

  long <- grep("_prior_shape$", names(opts), value = TRUE)
  short <- setdiff(grep("_shape$", names(opts), value = TRUE), long)

  for (nm in c(long, short)) {
    stem <- sub("_(prior_)?shape$", "", nm)
    rate <- sub("shape$", "rate", nm)

    if (is_null(opts[[rate]])) {
      next
    }

    # The proportional hazards baseline is one rate per time bin under the same
    # prior, which is worth saying since the others are single numbers.
    each <- {
      if (stem == "lambda" && !is_null(opts[["edges"]]))
        sprintf(", one for each of %s time bins", length(opts[["edges"]]) - 1L)
      else ""
    }

    add(stem, sprintf("Gamma(shape = %s, rate = %s)%s", value(nm), value(rate),
                      each),
        drawn(stem))
  }

  # The Dirichlet process mixture's base measure, which is where a weightless
  # observation's atom comes from, and the concentration that governs how many
  # atoms there are.
  if (family %in% c("dpm", "dpm_aft")) {
    add("atom variance",
        sprintf("%s * %s / ChiSq(%s)", value("nu"), value("lambda"),
                value("nu")),
        TRUE)
    add("atom mean",
        sprintf("Normal(%s, variance / %s)", value("mu_0"), value("k_0")),
        TRUE)
    grid <- opts[["alpha_grid"]]

    add("concentration",
        if (!drawn("alpha")) sprintf("fixed at %s", value("alpha"))
        else sprintf("tapered over a grid of %s values from %s to %s",
                     length(grid), signif(min(grid), 4L),
                     signif(max(grid), 4L)),
        drawn("alpha"))
  }

  if (family == "tweedie") {
    add("power",
        if (drawn("power")) "Uniform(1, 2)"
        else sprintf("fixed at %s", value("power")),
        drawn("power"))
  }

  if (family == "mnp") {
    add("covariance",
        sprintf("InverseWishart(%s, I), rescaled to unit mean variance",
                value("nu")),
        drawn("sigma"))
  }

  # The one thing worth saying about these is that there is nothing to say: the
  # cutpoints are drawn from the likelihood and have no prior, which is why
  # `prior_only = TRUE` refuses both families.
  if (family %in% c("ordinal", "ordbeta")) {
    add("cutpoints",
        if (family == "ordinal" && !drawn("cuts")) "fixed at the values given"
        else "none; drawn from the likelihood alone",
        family != "ordinal" || drawn("cuts"))
  }

  if (is_null(rows)) {
    return(NULL)
  }

  do_rbind(rows)
}

# Whether a setting was drawn, said in a way that survives a fit whose forests
# disagree about it: naming the ones it was drawn for is the only accurate thing
# to say when a per-forest argument was used.
estimated_phrase <- function(flags, forests) {
  if (all(flags)) {
    return("estimated")
  }

  if (!any(flags)) {
    return("held fixed")
  }

  sprintf("estimated for %s and held fixed for the rest",
          toString(forests[flags]))
}

# A per-forest setting written into the prose: the number itself when every
# forest was given the same one, and the name of the column carrying it when
# they were not. The column is in the table printed just above in that case, so
# the sentence reads as the formula it is and the numbers are one line away.
shared_value <- function(column, name, digits) {
  values <- unique(column)

  if (length(values) != 1L) {
    return(name)
  }

  format(round(values, digits))
}

#' @rdname bartisan-interop
#' @export
print.bartisan_prior_summary <- function(x, digits = 3L, ...) {

  arg::arg_whole_number(digits)

  table <- x[["forests"]]
  est <- x[["estimated"]]
  forests <- table[["forest"]]

  val <- function(nm) {
    shared_value(table[[nm]], nm, digits)
  }

  varies <- function(...) {
    any(vapply(c(...), function(nm) length(unique(table[[nm]])) > 1L, logical(1L)))
  }

  cli_cat("{.underline Priors}")
  cli::cat_line()

  # The table earns its place only when the forests were given different
  # settings, which is what a per-forest argument is for. Otherwise every number
  # in it appears in the prose below, and printing both says it twice.
  if (nrow(table) > 1L && varies(setdiff(names(table), "forest"))) {
    show <- table[c(TRUE, vapply(table[-1L], function(z) length(unique(z)) > 1L,
                                 logical(1L)))]

    for (nm in setdiff(names(show), "forest")) {
      show[[nm]] <- round(show[[nm]], digits)
    }

    print(show, row.names = FALSE)
    cli::cat_line()
  }

  cli_cat("{.underline Trees}")

  tree_qty <- if (varies("num_trees")) 2L else table[["num_trees"]][[1L]]

  # Worth an illustration only when one number governs every forest; two
  # different branching probabilities do not have one root probability.
  depth_note <- {
    if (varies("gamma", "beta")) ""
    else sprintf(", so the root splits with probability %s and a node at depth 3 with %s",
                 val("gamma"),
                 round(table[["gamma"]][[1L]] * 4^-table[["beta"]][[1L]],
                       digits))
  }

  cli_bullets_cat(c(
    "*" = "{val('num_trees')}{cli::qty(tree_qty)} tree{?s} per additive predictor, summed. A node at
           depth {.emph d} branches with probability {val('gamma')} *
           (1 + {.emph d})^-{val('beta')}{depth_note}."))
  cli::cat_line()

  cli_cat("{.underline Leaves}")
  cli_bullets_cat(c(
    "*" = "Each leaf value is Normal(0, {val('leaf_scale')}^2), that scale being
           3 * {.emph s} / ({val('k')} * sqrt({val('num_trees')})) with
           {.emph s} the response's scale on the link scale. The scale is itself
           given a half-Cauchy prior centred there and is
           {estimated_phrase(est[['leaf_scale']], forests)}."))
  cli::cat_line()

  cli_cat("{.underline Splitting variables}")

  if (x[["split_prior"]]) {
    cli_bullets_cat(c(
      "*" = "Fixed by {.arg split_prior}, so how the rules are shared out among
             the {val('candidates')} predictors is not drawn."))
  }
  else if (isFALSE(x[["sparsity"]])) {
    cli_bullets_cat(c(
      "*" = "Each of the {val('candidates')} predictors is equally likely to be
             split on, and that is not drawn
             ({.code sparsity = FALSE})."))
  }
  else {
    # Assembled here rather than inside the message, since cli parses a `{}`
    # expression on one line and an `if` written across several is a syntax
    # error by the time it gets there.
    sparsity_drawn <- {
      if (all(est[["splits"]]) && all(est[["alpha"]])) "Both are estimated."
      else if (identical(est[["splits"]], est[["alpha"]]))
        sprintf("Both are %s.", estimated_phrase(est[["splits"]], forests))
      else sprintf("The shares are %s and the concentration is %s.",
                   estimated_phrase(est[["splits"]], forests),
                   estimated_phrase(est[["alpha"]], forests))
    }

    cli_bullets_cat(c(
      "*" = "The share of the rules each of the {val('candidates')} predictors
             receives is Dirichlet({val('alpha')} / {val('candidates')}), whose
             concentration enters as
             a / (a + {val('alpha_scale')}) ~ Beta({val('shape_1')},
             {val('shape_2')}). {sparsity_drawn}"))

    if (x[["share_sparsity"]]) {
      cli_bullets_cat(c(
        "*" = "The forests share one set of shares
               ({.code share_sparsity = TRUE})."))
    }
  }

  cli::cat_line()
  cli_cat("{.underline Decision rules}")

  if (x[["rules"]][["soft"]]) {
    cli_bullets_cat(c(
      "*" = "Soft, with {.val {x[['rules']][['gate']]}} gates. Each tree's
             bandwidth is drawn from an exponential with mean
             {val('bandwidth')}, on predictors mapped to [0, 1], and is
             {estimated_phrase(est[['bandwidth']], forests)}."))
  }
  else {
    cli_bullets_cat(c(
      "*" = "Hard, as in standard BART, so no bandwidth is used."))
  }

  if (!is_null(x[["random"]])) {
    random <- x[["random"]]

    tau_drawn <- if (random[["estimated"]]) "it is drawn"
    else "it is held at that prior's median"

    cli::cat_line()
    cli_cat("{.underline Group intercepts}")
    cli_bullets_cat(c(
      "*" = "{cli::qty(length(random[['terms']]))}Term{?s}
             {.val {random[['terms']]}}, with
             {paste(random[['levels']], collapse = ' and ')} levels. Each
             level's intercept is Normal(0, tau^2), and tau is given the
             half-Cauchy prior the leaf scale gets, centred at
             {round(random[['scale']], digits)}; {tau_drawn}."))
  }

  cli::cat_line()

  family <- x[["family"]]

  cli_cat("{.underline Family}: {family[['name']]}, {family[['link']]} link")

  params <- family[["parameters"]]

  if (is_null(params)) {
    cli_bullets_cat(c(
      "*" = "No parameters of its own beyond the additive predictors above."))
  }
  else {
    for (i in seq_len(nrow(params))) {
      row <- params[i, ]
      drawn <- if (row[["estimated"]]) "(estimated)" else "(fixed)"

      cli_bullets_cat(c(
        "*" = "{row[['parameter']]}: {row[['prior']]} {drawn}"))
    }
  }

  cli::cat_line()
  cli_bullets_cat(c(
    i = "The leaf scale, and any number above read off the response, are
         calibrated rather than fitted; that is how a BART prior is specified.",
    i = "{.code prior_only = TRUE} in {.fn bartisan} draws from all of this, so
         that what it implies can be read on the outcome's own scale."))

  if (x[["prior_only"]]) {
    cli_bullets_cat(c(i = "This fit is itself a draw from the prior."))
  }

  invisible(x)
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
loo.bartisan_fit <- function(x, scale = NULL, ...) {
  prior_only_refuse(x, "loo")

  ll <- survival_measure(x, log_lik.bartisan_fit(x), scale)
  r_eff <- loo::relative_eff(exp(ll), chain_id = chain_ids(x))

  loo::loo.matrix(ll, r_eff = r_eff, ...)
}

#' @rdname bartisan-interop
#' @exportS3Method loo::waic
waic.bartisan_fit <- function(x, scale = NULL, ...) {
  prior_only_refuse(x, "waic")

  survival_measure(x, log_lik.bartisan_fit(x), scale) |>
    loo::waic.matrix(...)
}

#' @rdname bartisan-interop
#' @exportS3Method loo::kfold
kfold.bartisan_fit <- function(x, K = 10, folds = NULL, scale = NULL,
                               save_fits = FALSE, ...) {

  prior_only_refuse(x, "kfold")

  rlang::check_installed("loo", "for K-fold cross-validation.")

  arg::arg_flag(save_fits)

  n <- nobs(x)
  data <- kfold_data(x)
  rows <- rownames(x[["model"]])

  folds <- kfold_folds(folds, K, n)
  K <- length(unique(folds))

  # Each fold is scored by a fit that never saw it, so the call is rebuilt
  # against the training rows. The original *expression* for the formula is what
  # is re-evaluated, not `x[["formula"]]`: that one has already had its bars
  # replaced and its `vc()` terms reduced to names, so refitting from it would
  # silently drop the random-effect and varying-coefficient structure.
  base <- kfold_call(x, data)

  # Prior weights and the offset come from the stored model frame rather than
  # being re-evaluated, since either may have lived in the caller's workspace
  # rather than in `data`. They have to be carried into the *score* as well:
  # `predict(type = "density")` falls back to the fit's own weights when none
  # are given, so a weighted fit scored without them is wrong rather than an
  # error.
  weights <- stats::model.weights(x[["model"]])
  offset <- stats::model.offset(x[["model"]])

  one_fold <- function(k) {
    train <- rows[folds != k]
    held <- rows[folds == k]

    call <- base
    call[["data"]] <- data[train, , drop = FALSE]

    if (!is_null(weights)) {
      call[["weights"]] <- weights[folds != k]
    }

    if (!is_null(offset)) {
      call[["offset"]] <- offset[folds != k]
    }

    fit <- eval(call)

    score <- stats::predict(fit, newdata = data[held, , drop = FALSE],
                            type = "density", log = TRUE,
                            weights = weights[folds == k],
                            offset = offset[folds == k])

    list(score = score, fit = if (save_fits) fit)
  }

  # The folds are independent refits, so they go to workers when a plan has
  # any, with the streams drawn here for the reason `estimate_effect()` gives.
  seeds <- parallel_streams(K)

  if (use_future()) {
    done <- future.apply::future_lapply(seq_len(K), one_fold, future.seed = seeds,
                                        future.packages = "bartisan")
  }
  else {
    restore <- restore_stream()
    on.exit(restore(), add = TRUE)

    done <- lapply(seq_len(K), function(k) {
      assign(".Random.seed", seeds[[k]], envir = globalenv())
      one_fold(k)
    })
  }

  elpd <- numeric(n)

  for (k in seq_len(K)) {
    elpd[folds == k] <- done[[k]][["score"]]
  }

  # `p_kfold` is the gap between what the model predicts for an observation it
  # was fitted to and what it predicts for the same one held out. Both sides
  # have to be on the same measure or the difference is not that gap, so
  # `scale` is applied to each rather than only to the held-out side.
  lpd <- stats::predict(x, type = "density", log = TRUE)

  if (!is_null(scale)) {
    shift <- survival_shift(x, scale)
    elpd <- elpd - shift
    lpd <- lpd - shift
  }

  kfold_object(elpd, lpd, folds, K, nrow(x[["sigma_mu"]]),
               if (save_fits) pluck(done, "fit"))
}

# The data the fit was made from, recovered the way `stats::update()` recovers
# it: the call's own `data` expression, evaluated in the environment the formula
# carries.
kfold_data <- function(x) {
  expr <- x[["call"]][["data"]]

  if (is_null(expr)) {
    arg::err(c("This fit's call names no {.arg data}, so the folds have nothing
                to be taken from.",
               i = "Refit with {.arg data} given as a data frame."))
  }

  out <- eval(expr, environment(stats::formula(x)))

  if (!is.data.frame(out)) {
    out <- as.data.frame(out)
  }

  missing <- setdiff(rownames(x[["model"]]), rownames(out))

  if (!is_null(missing)) {
    arg::err(c("The data this fit was made from is not the data that name now
                reaches: {length(missing)} of its rows are gone.",
               i = "K-fold refits from the original call, so the data has to be
                    what it was."))
  }

  out
}

# Fold assignments, either the caller's or drawn at random.
kfold_folds <- function(folds, K, n) {
  if (is_null(folds)) {
    arg::arg_count(K)
    arg::arg_gte(K, 2)
    arg::arg_lte(K, n)

    return(loo::kfold_split_random(K = K, N = n))
  }

  arg::arg_numeric(folds)

  if (length(folds) != n) {
    arg::err("{.arg folds} must give one fold per observation, and gives
              {length(folds)} for {n}")
  }

  folds <- as.integer(folds)

  if (anyNA(folds) || min(folds) < 1L) {
    arg::err("{.arg folds} must be whole numbers from 1 up, with no missing
              values")
  }

  if (!setequal(folds, seq_len(max(folds)))) {
    arg::err("{.arg folds} must use every fold from 1 to {max(folds)}, and
              leaves at least one empty")
  }

  folds
}

# The call a fold refits from: the original, with everything but the data
# resolved to a value so that a worker needs nothing from the caller's
# environment. The function in position one is resolved too, without which a
# worker cannot find `bartisan()` at all.
kfold_call <- function(x, data) {
  env <- environment(stats::formula(x))
  call <- x[["call"]]

  call[[1L]] <- eval(call[[1L]], env)

  for (nm in setdiff(names(call), c("", "data", "subset", "weights",
                                    "offset"))) {
    call[[nm]] <- eval(call[[nm]], env)
  }

  # The rows are chosen by name below, so a `subset` would choose them twice.
  # Assigning NULL to a name a call does not have is an error, hence the guard.
  if ("subset" %in% names(call)) {
    call[["subset"]] <- NULL
  }

  # K refits would otherwise report progress K times over.
  if (!"verbose" %in% names(call)) {
    call[["verbose"]] <- FALSE
  }

  call
}

# The three columns \pkg{loo} expects, and the shape its own
# `table_of_estimates()` produces: the estimate is the sum over observations and
# the standard error is the spread of the pointwise values scaled by the count.
#
# Three rows rather than one because `loo_compare()` flattens each object's
# estimates and binds them: a one-row matrix beside a `loo` object's three makes
# `sapply()` return a list rather than a matrix, and the comparison fails.
kfold_object <- function(elpd, lpd, folds, K, draws, fits = NULL) {
  pointwise <- cbind(elpd_kfold = elpd,
                     p_kfold = lpd - elpd,
                     kfoldic = -2 * elpd)

  estimates <- cbind(Estimate = colSums(pointwise),
                     SE = sqrt(nrow(pointwise) * apply(pointwise, 2L,
                                                       stats::var)))

  out <- list(estimates = estimates, pointwise = pointwise, folds = folds)

  if (!is_null(fits)) {
    out[["fits"]] <- fits
  }

  attr(out, "K") <- as.integer(K)
  attr(out, "dims") <- c(draws, nrow(pointwise))

  class(out) <- c("kfold", "loo")
  out
}

# A prior-only fit has no likelihood to score, so everything built on one has to
# say so rather than return a number. `model_performance()` is included because
# it reaches `loo()` on its own and the refusal would otherwise surface from a
# call that never mentioned it.
prior_only_refuse <- function(object, what) {
  if (!isTRUE(object[["prior_only"]])) {
    return(invisible(TRUE))
  }

  arg::err(c("{.fn {what}} scores a fit against the data, and this fit was made
              with {.code prior_only = TRUE}, so it was never shown any.",
             i = "Refit without {.arg prior_only} to score it."))
}

# The survival families do not all write their likelihood with respect to the
# same measure. An accelerated failure time family reports the density of
# \eqn{\log T} and `ph()` the density of \eqn{T}, so the two differ by the
# Jacobian of the change of variable and a log score taken across that boundary
# is off by \eqn{\sum \log t}, which is large enough to reverse an ordering.
#
# The correction is not applied on its own initiative, because `loo()` would
# then stop reporting the model's own predictive density: it would no longer
# agree with `log_lik()`, and a comparison against a proportional hazards fit
# from another package would silently acquire the error this is meant to
# remove. Naming the scale is what asks for it, and it reads the same from
# either side, since whichever family is already there is left alone.
survival_measure <- function(object, ll, scale) {
  shift <- survival_shift(object, scale)

  if (is_null(shift)) {
    return(ll)
  }

  sweep(ll, 2L, shift, "-")
}

# The per-observation constant the change of variable adds, or `NULL` when there
# is none to add because the fit is already on the scale asked for. Split out of
# the sweep above so that `kfold()` can apply the same shift to two vectors: the
# held-out score and the in-sample one it is differenced against have to be on
# one measure or their difference is not the quantity `p_kfold` names.
survival_shift <- function(object, scale) {
  if (is_null(scale)) {
    return(NULL)
  }

  scale <- arg::match_arg(scale, c("time", "log_time"))

  family <- object[["family"]][["family"]]

  if (!family %in% c("aft", "dpm_aft", "ph")) {
    arg::err(c("{.arg scale} names the measure a survival model's density is
                taken with respect to, and this fit's family is
                {.val {family}}.",
               i = "Leave it empty; only {.fn ph} and the accelerated failure
                    time families report on two different scales."))
  }

  on_log_time <- !identical(family, "ph")

  if (identical(scale, "log_time") == on_log_time) {
    return(NULL)
  }

  # Events only. A censored observation contributes a survival probability,
  # which is a probability on either scale and has no measure to change.
  y <- stats::model.response(object[["model"]])
  shift <- as.numeric(y[, "status"]) * log(as.numeric(y[, "time"]))

  # Going the other way is the same constant with the other sign.
  if (identical(scale, "time")) shift else -shift
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

  # A residual variance needs every observed value to be a draw from the fitted
  # mean, and a censored time is a bound on one. The number this produced was
  # therefore not an R2, and for `dpm_aft()` it also subtracted a median
  # survival time from a log time.
  if (family %in% c("aft", "dpm_aft", "ph")) {
    if (verbose) {
      arg::wrn("the {.val {family}} family's response is a censored time, so a
                residual variance, and with it a Bayesian {.field R2}, is not
                defined for it")
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

  # A prior-only fit has no likelihood, so the three scores built on one are
  # dropped rather than erroring: the rest of the table is computable and the
  # caller asked for the table, not for the ELPD.
  if (isTRUE(model[["prior_only"]])) {
    dropped <- intersect(c("ELPD", "LOOIC", "WAIC"), metrics)

    if (verbose && !is_null(dropped)) {
      arg::msg(c(i = "Leaving out {.val {dropped}}: this fit was made with
                      {.code prior_only = TRUE}, so there is no likelihood to
                      score."))
    }

    metrics <- setdiff(metrics, c("ELPD", "LOOIC", "WAIC"))
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
          c("ordinal", "multinomial", "mnp", "custom", "dpm_aft", "ph")) NULL
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
