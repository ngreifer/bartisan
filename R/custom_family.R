#' Fit a model to a likelihood written in R
#'
#' @description
#' `custom_family()` takes the log density itself, as an R function, and returns
#' a family object [bartisan()] fits the model for. It is the route to a
#' response distribution none of the compiled families in
#' [`bartisan-families`] cover.
#'
#' @param logdens the log density, given as a function of the response and the
#'   additive predictors, `function(y, eta)`, where `y` is a numeric vector of
#'   length `n` and `eta` an `n` by `num_predictors` matrix, returning a numeric
#'   vector of length `n`. With nuisance parameters it takes a third argument,
#'   `function(y, eta, aux)`, where `aux` is a numeric vector of their current
#'   values. It is the log density of one unit of prior weight, so that
#'   `weights` behave as they do elsewhere, and terms free of `eta` may be
#'   dropped.
#' @param num_predictors `numeric`; how many additive predictors the density has
#'   (i.e., how many forests to fit). Default is 1.
#' @param start `numeric`; the value each additive predictor starts at, in place
#'   of the intercept-only fit the compiled families use. One value or one per
#'   predictor. Default is 0.
#' @param derivatives optional; a `function(y, eta, h)` returning a list with
#'   elements `score` and `info`, the first derivative of `logdens` with respect
#'   to the `h`th predictor and minus its second derivative, each a vector of
#'   length `n`. Default is `NULL` to take central differences of `logdens`. It
#'   covers the additive predictors only: a nuisance parameter is always
#'   differenced, which costs three calls per sweep rather than three per leaf.
#' @param aux_names optional `character`; the names of the nuisance parameters
#'   to draw, if any. Naming them is what declares them, because the names label
#'   the columns of `fit$aux` and are what `summary()` and [diagnose()] report
#'   them under. They must be distinct and non-empty. Default is `NULL` for
#'   none, unless `aux_start` is given, in which case the parameters are named
#'   by `names(aux_start)` when it carries names and positionally (`"aux1"`,
#'   `"aux2"`, and so on) when it does not.
#' @param aux_start optional `numeric`; the value each nuisance parameter starts
#'   at, given as one value or one per parameter. Default is `NULL`, which is 0
#'   for each. The sampler will walk to the posterior from a poor start, so this
#'   need only be the right order of magnitude. Supplying it is a second way to
#'   declare the parameters, so `aux_start = c(shape = 1)` both names one and
#'   starts it at 1.
#' @param name string; a label used when printing the fit. Default is
#'   `"custom"`.
#'
#' @returns
#' A `<bartisan_family>` object, which is a list containing at least the
#' elements `family` and `link` and which inherits from `family`, so that
#' [bartisan()] recognizes it wherever it recognizes an ordinary
#' [`stats::family`] object.
#'
#' @details
#' Nothing else about the sampler changes: the leaf-level Laplace proposal needs
#' the first two derivatives of the log density with respect to each additive
#' predictor and nothing more, and central differences of the supplied function
#' produce both.
#'
#' ```r
#' # A Poisson model written out by hand. Terms free of eta may be dropped;
#' # they cancel from every acceptance ratio.
#' pois <- custom_family(function(y, eta) y * eta[, 1] - exp(eta[, 1]),
#'                       start = log(mean(d$y)))
#'
#' # Two predictors: a mean and a log standard deviation.
#' ls <- custom_family(function(y, eta) dnorm(y, eta[, 1], exp(eta[, 2]),
#'                                            log = TRUE),
#'                     num_predictors = 2, start = c(0, 0))
#' ```
#'
#' The function is called once per leaf per Fisher-scoring step with the
#' observations reaching that leaf, so it must be vectorized over `y` and the
#' rows of `eta`; it must not be vectorized *within* an observation, and it must
#' return exactly one value per row. Supplying `derivatives` cuts three calls to
#' one and removes the differencing error.
#'
#' Because it sees a subset rather than the whole sample, anything else the
#' density needs has to be a scalar or reach it through `eta`. A
#' per-observation vector captured from the enclosing environment will not line
#' up with the rows it is handed, and nothing can detect that for you: the
#' fit runs and is wrong. Where such a quantity is genuinely needed, an offset
#' belongs in the formula and a varying trial count or exposure belongs in a
#' family written for it.
#'
#' ## Nuisance Parameters
#'
#' These are drawn alongside the trees when `aux_names` names them, and
#' `logdens` then takes a third argument holding their current values:
#'
#' ```r
#' # A Gaussian written out by hand, with its scale drawn rather than fixed.
#' by_hand <- custom_family(
#'   logdens = function(y, eta, aux) dnorm(y, eta[, 1], exp(aux[1]), log = TRUE),
#'   aux_names = "log_sigma", aux_start = 0)
#' ```
#'
#' They are reported in `fit$aux` under those names, and covered by `summary()`
#' and [diagnose()] like any other family's. There is no prior argument and no
#' bounds argument, because a nuisance parameter here is carried as an additive
#' predictor whose forest is pinned at depth zero (one tree that can never
#' split, so the forest is a single scalar), and it is drawn by the same
#' Laplace-plus-Metropolis step as any leaf, under that step's Gaussian leaf
#' prior. So a parameter with a restricted range is handled the way it would be
#' for a real predictor, by writing the transform into `logdens`: the `exp()`
#' above is what keeps the scale positive.
#'
#' ## Limits of a Log Density
#'
#' A density says how likely an observed value is, not how to draw a new one, so
#' a `custom_family()` fit has no posterior predictive distribution, which is
#' what [`simulate()`][bartisan-interop], [`pp_check()`][bartisan-interop] and
#' [`r2()`][bartisan-interop] rely on. Each of them errors on such a fit
#' rather than returning something it cannot support. [loo()][bartisan-interop]
#' and [`waic()`][bartisan-interop] are unaffected, since both read the
#' pointwise log likelihood the family already computes.
#'
#' There are smaller limits as well. The response must be numeric, so a factor
#' has to be coded
#' first. And since the package cannot know what the mean of the density is,
#' `predict(type = "response")` returns the additive predictors rather than a
#' fitted mean.
#'
#' @seealso
#' [`bartisan-families`] for the compiled families, one of which is usually the
#' better answer; [bartisan()] for fitting a model with the result;
#' `vignette("families", package = "bartisan")` for the long form
#'
#' @examples
#' set.seed(123)
#'
#' d <- data.frame(x1 = runif(300), x2 = runif(300))
#' d$y <- rpois(300, exp(1 + sin(pi * d$x1)))
#'
#' # A Poisson likelihood written out by hand, with its derivatives. Terms
#' # free of eta may be dropped, since they cancel from every acceptance
#' # ratio, and the same terms are absent from the score.
#' pois <- custom_family(
#'   function(y, eta) y * eta[, 1] - exp(eta[, 1]),
#'   derivatives = function(y, eta, h) {
#'     list(score = y - exp(eta[, 1]), info = exp(eta[, 1]))
#'   },
#'   start = log(mean(d$y)))
#'
#' fit <- bartisan(y ~ x1 + x2, data = d, family = pois,
#'                 num_trees = 20, num_burn = 100, num_draws = 100,
#'                 verbose = FALSE)
#'
#' fit
#'
#' # A beta-binomial, which no built-in family covers: counts out of a known
#' # number of trials, overdispersed relative to a binomial. `phi` is drawn
#' # alongside the trees, and `size` is a scalar, so closing over it is safe.
#' size <- 30
#' d$hits <- rbinom(300, size, rbeta(300, plogis(d$x1) * 6,
#'                                   (1 - plogis(d$x1)) * 6))
#'
#' bb <- custom_family(
#'   logdens = function(y, eta, aux) {
#'     p <- plogis(eta[, 1])
#'     phi <- exp(aux[1])
#'     lbeta(y + p * phi, size - y + (1 - p) * phi) -
#'       lbeta(p * phi, (1 - p) * phi)
#'   },
#'   aux_names = "log_phi", aux_start = log(5), name = "beta-binomial")
#'
#' fit_bb <- bartisan(hits ~ x1 + x2, data = d, family = bb,
#'                    num_trees = 20, num_burn = 100, num_draws = 100,
#'                    verbose = FALSE)
#'
#' # The drawn precision, on the scale it was written on.
#' exp(mean(fit_bb$aux[, "log_phi"]))
#'
#' @export
custom_family <- function(logdens, num_predictors = 1L, start = 0,
                          derivatives = NULL, aux_names = NULL, aux_start = NULL,
                          name = "custom") {
  if (!is.function(logdens)) {
    arg::err("{.arg logdens} must be a function of the response and the
              additive predictors")
  }

  arg::arg_count(num_predictors)
  arg::arg_gte(num_predictors, 1)
  arg::arg_numeric(start)
  arg::arg_string(name)

  arg::when_not_null(
    derivatives,
    arg::arg_function
  )

  num_predictors <- as.integer(num_predictors)

  if (length(start) != 1L && length(start) != num_predictors) {
    arg::err("{.arg start} must have one value, or one per additive predictor
              ({num_predictors})")
  }

  # The nuisance parameters are declared by naming them, because their names are
  # what labels the columns of `fit$aux` and what `summary()` and `diagnose()`
  # report them under. Giving only starting values names them, from
  # `names(aux_start)` where it carries them and positionally where it does not.
  # Either route has to produce usable column names, so the same check applies
  # to both.
  arg::when_not_null(
    aux_names,
    arg::arg_character
  )

  arg::when_not_null(
    aux_start,
    arg::arg_numeric
  )

  if (!is_null(aux_names)) {
    if (anyDuplicated(aux_names) > 0L || !all(nzchar(aux_names))) {
      arg::err("{.arg aux_names} must be distinct and non-empty")
    }
  }
  else if (!is_null(aux_start)) {
    aux_names <- names(aux_start)

    if (is_null(aux_names)) {
      aux_names <- paste0("aux", seq_along(aux_start))
    }
    else if (anyDuplicated(aux_names) > 0L || !all(nzchar(aux_names))) {
      arg::err("the names of {.arg aux_start} must be distinct and non-empty")
    }
  }

  if (!is_null(aux_names)) {
    if (is_null(aux_start)) {
      aux_start <- rep.int(0, length(aux_names))
    }
    else if (length(aux_start) == 1L) {
      aux_start <- rep.int(aux_start, length(aux_names))
    }
    else if (length(aux_start) != length(aux_names)) {
      arg::err("{.arg aux_start} must have one value, or one per nuisance
              parameter ({length(aux_names)})")
    }

    aux_start <- unname(aux_start)
  }

  num_aux <- length(aux_names)

  # Empty but typed, rather than `NULL`, so that a family with no nuisance
  # parameters and one with some are the same shape of object.
  aux_names <- aux_names %or% character()
  aux_start <- aux_start %or% numeric()

  if (num_aux > 0L && length(formals(logdens)) < 3L) {
    arg::err(c("{.arg logdens} must take a third argument for the nuisance
                parameters when there are any.",
               i = "It is called as {.code logdens(y, eta, aux)}, with
                    {.arg aux} a numeric vector of length {num_aux}."))
  }

  new_bartisan_family("custom", "identity", logdens = logdens,
                      num_predictors = num_predictors,
                      start = rep(start, length.out = num_predictors),
                      derivatives = derivatives,
                      num_aux = num_aux,
                      aux_names = aux_names,
                      aux_start = aux_start,
                      name = name)
}
