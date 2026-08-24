#' Response families for generalized BART
#'
#' @description
#' `genbart()` accepts the [stats::family] objects used by [stats::glm()], so
#' `gaussian()`, `binomial("probit")`, `poisson()` and `Gamma("log")` all work
#' unchanged. The functions documented here supply the additional families that
#' have no `glm()` counterpart, in the same style, so that they can be passed to
#' the `family` argument the same way.
#'
#' @param link the link function. Each family compiles the links for which the
#'   additive predictor is the natural unconstrained scale; any other link is
#'   applied from R, for the families where that is well defined. See Details.
#' @param theta for `negbin()` and `zi_negbin()`, a fixed value for the
#'   dispersion parameter. The default, `NULL`, draws it along with everything
#'   else.
#' @param shape for `Gamma_shape()`, a fixed value for the gamma shape. The
#'   default, `NULL`, draws it.
#' @param phi for `ordbeta()`, a fixed value for the beta precision. The
#'   default, `NULL`, draws it.
#' @param reference for `multinomial()`, the response category to hold as the
#'   reference. The default, `NULL`, fits one forest per category instead and
#'   leaves the model unidentified, which is what makes the prior symmetric in
#'   the categories; see Details.
#' @param logdens for `custom_family()`, the log density. A function of the
#'   response and the additive predictors, `function(y, eta)`, where `y` is a
#'   numeric vector of length `n` and `eta` an `n` by `num_predictors` matrix,
#'   returning a numeric vector of length `n`. It is the log density of *one
#'   unit of prior weight*, so that `weights` behave as they do elsewhere, and
#'   terms free of `eta` may be dropped.
#' @param num_predictors for `custom_family()`, how many additive predictors the
#'   density has, that is, how many forests to fit.
#' @param start for `custom_family()`, the value each additive predictor starts
#'   at, in place of the intercept-only fit the compiled families use. One value
#'   or one per predictor.
#' @param derivatives for `custom_family()`, an optional
#'   `function(y, eta, h)` returning a list with elements `score` and `info`,
#'   the first derivative of `logdens` with respect to the `h`th predictor and
#'   minus its second derivative, each a vector of length `n`. The default,
#'   `NULL`, takes central differences of `logdens`.
#' @param name for `custom_family()`, a label used when printing the fit.
#'
#' @details
#' Every family reduces to a scalar additive predictor, or to several of them,
#' together with the first two derivatives of the log density with respect to
#' each. That is the whole interface the sampler needs, which is why the set of
#' available families is not restricted to the conditionally conjugate ones.
#'
#' The supported families and links are:
#'
#' | Family | Links | Additive predictors | Drawn nuisance parameters |
#' |---|---|---|---|
#' | `gaussian()` | `identity` | 1 | residual standard deviation |
#' | `binomial()` | `logit`, `probit`, `cloglog` | 1 | none |
#' | `poisson()` | `log` | 1 | none |
#' | `negbin()` | `log` | 1 | dispersion |
#' | `Gamma()`, `Gamma_shape()` | `log` | 1 | shape |
#' | `ordinal()` | `logit`, `probit`, `cloglog` | 1 | cutpoints |
#' | `multinomial()` | `logit` | one per category, or per non-reference level | none |
#' | `weibull_aft()`, `loglogistic_aft()`, `lognormal_aft()` | none | 1 | scale |
#' | `location_scale()` | `identity` | 2 | none |
#' | `zi_poisson()` | `log` | 2 | none |
#' | `zi_negbin()` | `log` | 2 | dispersion |
#' | `ordbeta()` | `logit` | 1 | 2 cutpoints, precision |
#'
#' # Links the engine does not compile
#'
#' The links listed above are the ones the sampler evaluates in compiled code.
#' Any other link is accepted for `gaussian()`, `binomial()`, `poisson()`,
#' `negbin()` and `Gamma()`, and applied from R: the additive predictor is
#' mapped to the scale the compiled family works on -- the mean, the log mean or
#' the log odds -- by composing the caller's inverse link with the family's own,
#' and the chain rule carries the derivatives back. So
#' `binomial("cauchit")` works, as does any link object of the kind
#' [stats::make.link()] returns, including one written by hand:
#'
#' ```r
#' my_link <- structure(
#'   list(linkfun = function(mu) qnorm(mu) / 2,
#'        linkinv = function(eta) pnorm(2 * eta),
#'        mu.eta = function(eta) 2 * dnorm(2 * eta),
#'        valideta = function(eta) TRUE, name = "half-probit"),
#'   class = "link-glm")
#'
#' genbart(y ~ ., data = d, family = binomial(my_link))
#' ```
#'
#' Three things are worth knowing about this route. It costs a call into R for
#' every leaf the sampler visits, so a fit is slower than one with a compiled
#' link, though not by the factor a call per observation would cost. The leaf
#' prior scale is calibrated for the compiled link, so an unusual link may want
#' `sigma_mu` set by hand in [genbart_control()]. And the additive predictor is
#' unconstrained, so a link whose inverse has a restricted range --
#' `poisson("identity")`, `Gamma("inverse")` -- will produce non-finite
#' densities for some predictors; those proposals are rejected rather than
#' breaking the chain, but the fit will be poor. Prefer links whose inverse is
#' defined on the whole line.
#'
#' The families with more than one additive predictor, or whose link enters
#' somewhere other than a single mean -- `ordinal()`, `multinomial()`, the
#' accelerated failure time families, `location_scale()`, the zero-inflated
#' families and `ordbeta()` -- take only their listed links. `custom_family()`
#' is the way to reach anything else.
#'
#' # Supplying a likelihood
#'
#' `custom_family()` takes the log density itself, as an R function, and fits
#' the model that goes with it. Nothing else about the sampler changes: the
#' leaf-level Laplace proposal needs the first two derivatives of the log
#' density with respect to each additive predictor and nothing more, and central
#' differences of the supplied function produce both.
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
#' one and removes the differencing error, and is worth doing when the
#' derivatives are easy to write down.
#'
#' What `custom_family()` does not do: it cannot draw a nuisance parameter, so a
#' dispersion has to be fixed inside the closure; the response must be numeric,
#' so a factor has to be coded first; and since the package cannot know what the
#' mean of the density is, `predict(type = "response")` returns the additive
#' predictors rather than a fitted mean.
#'
#' # Families in detail
#'
#' `ordinal()` expects an ordered factor, though a numeric or integer response
#' is accepted and its sorted unique values are taken as the categories. It uses
#' the cumulative-link parameterization of [MASS::polr()], in which
#' \eqn{P(Y \le k) = F(c_k - \eta)}, so that larger values of the additive
#' predictor shift mass towards higher categories.
#'
#' Only the differences \eqn{c_k - \eta_i} are identified, so one location has
#' to be pinned. With three or more categories the draws are reported in the
#' chart where **the additive predictor has mean zero over the fitted sample and
#' every cutpoint is free**, which is the convention [MASS::polr()] reports in
#' and makes the cutpoints directly readable as category boundaries. With exactly
#' two categories the single boundary is folded into the intercept instead, so a
#' two-category response is exactly binary regression with the matching link and
#' on the same scale.
#'
#' The choice of chart is not a choice of model: the sampler works in whichever
#' is better conditioned and the shift is applied to the recorded predictor, the
#' recorded cutpoints and the recorded leaf values together, so every category
#' probability is untouched and the stored forest still replays to the reported
#' predictor. It does mean that `cut1` is a free parameter rather than a constant
#' zero, which is a change from earlier versions.
#'
#' **Comparing with [MASS::polr()].** This is the chart `polr()` reports in when
#' its predictors are centered. `polr()` identifies the location by leaving the
#' intercept out of the design matrix rather than by centering, so its `zeta` is
#' shifted by the mean of its own linear predictor, which is not zero unless the
#' predictors happen to be. Either of these puts the two side by side:
#'
#' ```r
#' fit <- genbart(y ~ x1 + x2, data = d, family = ordinal("probit"))
#' colMeans(fit$aux)                       # cutpoints, predictor centered
#'
#' p <- MASS::polr(y ~ x1 + x2, data = d, method = "probit")
#' p$zeta - mean(p$lp)                     # the same chart
#'
#' d2 <- transform(d, x1 = x1 - mean(x1), x2 = x2 - mean(x2))
#' MASS::polr(y ~ x1 + x2, data = d2, method = "probit")$zeta   # likewise
#' ```
#'
#' On a linear truth with n = 3000 the three agree to about 0.02, which is Monte
#' Carlo error; there is a test to that effect. The *gaps* between cutpoints are
#' identified outright and match whatever chart either is in.
#'
#' All three links are fitted through a latent variable, which
#' [genbart_control()]'s `augment` uses by default: conditional on it the target
#' over a leaf is quadratic and the sampler takes the closed form instead of
#' Fisher scoring plus a Metropolis ratio. For `ordinal("probit")` the latent
#' variable is normal (Albert and Chib 1993) and the measured gain on a thousand
#' observations and fifty trees is 14 times with soft rules and 30 times with
#' hard ones. For `ordinal("logit")` it is logistic, which is a normal whose
#' precision is Polya-Gamma, and the gain is 7 times with soft rules and 15 with
#' hard ones -- smaller because the precision is an extra draw per observation.
#' For `ordinal("cloglog")` it is an exponential waiting time, because the
#' complementary log-log model *is* the discrete proportional hazards model; the
#' gain is 5 times with hard rules and about 1.6 times in effective samples per
#' second with soft ones, since that target is the exponential form rather than
#' the quadratic one and the exponential form needs hard rules. Accuracy is
#' unchanged in every case. So the link can be chosen on modelling grounds; see
#' [genbart_control()] for the details and for how to turn the rewriting off.
#'
#' `Gamma_shape()` is the gamma family again, and differs from `stats::Gamma()`
#' only in taking a `shape` argument: the shape, which acts as the inverse
#' dispersion, is drawn by default and can instead be fixed. It does *not*
#' regress the shape on the predictors. `negbin()` and `ordbeta()` take
#' `theta` and `phi` the same way.
#'
#' `multinomial()` expects an unordered factor. By default it fits one forest
#' per category and leaves the model unidentified, since adding any function of
#' the predictors to every category's forest leaves the probabilities alone.
#' This is the parameterization of Murray (2021): the leaf prior is proper, so
#' the posterior is proper too, and every identified quantity -- a probability,
#' an odds ratio -- is recovered from the draws. Its point is that the prior is
#' then symmetric in the categories, whereas reference coding makes the fit
#' depend on which category was singled out. The leaf prior scale is divided by
#' \eqn{\sqrt 2} to compensate for each log-odds contrast now being a
#' difference of two forests, which leaves the prior on the identified log odds
#' where reference coding puts it.
#'
#' Passing `reference` instead pins that category's predictor at zero and fits
#' `num_cat - 1` forests, so the fitted functions are log odds against it and
#' are interpretable the way the coefficients from [nnet::multinom()] are. That
#' is one fewer forest, and the right choice when a particular contrast is the
#' quantity of interest.
#'
#' The accelerated failure time families expect a right-censored response,
#' supplied either as a [survival::Surv()] object or as a two-column matrix of
#' times and event indicators. They model \eqn{\log T = \eta + \sigma\epsilon}
#' with \eqn{\epsilon} standard Gumbel, logistic or normal respectively, giving
#' Weibull, log-logistic and log-normal survival times.
#'
#' `location_scale()` regresses the mean and the log standard deviation of a
#' normal response on separate forests, so the variance is an unrestricted
#' function of the predictors.
#'
#' `zi_poisson()` and `zi_negbin()` are zero-inflated counts, mixing a point
#' mass at zero with a Poisson or negative binomial count component. Both parts
#' get their own forest: the first predictor is the log mean of the count
#' component and the
#' second the log odds that an observation is a structural zero, so the
#' excess-zero mechanism is free to depend on the predictors rather than being a
#' single constant. The two are reported as the `count` and `zero` predictors.
#'
#' `ordbeta()` is the ordered beta regression of Kubinec (2023), for a response
#' on the closed unit interval with point masses at zero and one, such as a
#' proportion or a slider scale. One predictor drives both the probability of
#' landing on an endpoint, through a pair of cutpoints as in an ordinal model,
#' and the mean of the beta density in between:
#' \deqn{P(Y = 0) = 1 - \mathrm{logit}^{-1}(\eta - c_1),}
#' \deqn{P(Y = 1) = \mathrm{logit}^{-1}(\eta - c_2),}
#' with the remaining mass following a
#' \eqn{\mathrm{Beta}(\mu\phi, (1 - \mu)\phi)} density for
#' \eqn{\mu = \mathrm{logit}^{-1}(\eta)}. Because the predictor also enters
#' the beta mean it is identified, so unlike `ordinal()` both cutpoints are
#' drawn.
#'
#' Proportional hazards regression is deliberately absent. Its partial
#' likelihood couples observations through risk sets, so it does not decompose
#' into a sum over the observations reaching a leaf, which is what the leafwise
#' Laplace approximation requires. The accelerated failure time families cover
#' the same ground within this framework.
#'
#' @returns
#' A list of class `genbart_family`, containing at least the elements `family`
#' and `link`. Objects of this class are recognized by `genbart()` alongside
#' ordinary [stats::family] objects.
#'
#' @references
#' Kubinec, R. (2023). Ordered beta regression: a parsimonious, well-fitting
#' model for continuous data with lower and upper bounds. *Political Analysis*,
#' 31(4), 519--536. \doi{10.1017/pan.2022.20}
#'
#' Murray, J. S. (2021). Log-linear Bayesian additive regression trees for
#' multinomial logistic and count regression models. *Journal of the American
#' Statistical Association*, 116(534), 756--769.
#' \doi{10.1080/01621459.2020.1813587}
#'
#' @seealso [genbart()]
#'
#' @examplesIf FALSE
#' genbart(y ~ ., data = d, family = negbin())
#' genbart(y ~ ., data = d, family = ordinal("probit"))
#' genbart(survival::Surv(time, status) ~ ., data = d, family = weibull_aft())
#' genbart(count ~ ., data = d, family = zi_negbin())
#' genbart(proportion ~ ., data = d, family = ordbeta())
#'
#' # A link the engine does not compile, applied from R.
#' genbart(y ~ ., data = d, family = binomial("cauchit"))
#'
#' # A likelihood supplied from R.
#' genbart(y ~ ., data = d,
#'         family = custom_family(function(y, eta) y * eta[, 1] - exp(eta[, 1])))
#'
#' @name genbart-families
NULL

new_genbart_family <- function(family, link, ...) {
  structure(c(list(family = family, link = link), list(...)),
            class = c("genbart_family", "family"))
}

#' @rdname genbart-families
#' @export
negbin <- function(link = "log", theta = NULL) {
  link <- arg::match_arg(link, "log")

  if (!is_null(theta)) {
    arg::arg_number(theta)
    arg::arg_gt(theta, 0)
  }

  new_genbart_family("negbin", link, theta = theta)
}

#' @rdname genbart-families
#' @export
Gamma_shape <- function(link = "log", shape = NULL) {
  link <- arg::match_arg(link, "log")

  if (!is_null(shape)) {
    arg::arg_number(shape)
    arg::arg_gt(shape, 0)
  }

  new_genbart_family("Gamma", link, shape = shape)
}

#' @rdname genbart-families
#' @export
ordinal <- function(link = "logit") {
  link <- arg::match_arg(link, c("logit", "probit", "cloglog"))

  new_genbart_family("ordinal", link)
}

#' @rdname genbart-families
#' @export
multinomial <- function(link = "logit", reference = NULL) {
  link <- arg::match_arg(link, "logit")

  if (!is_null(reference)) {
    if (length(reference) != 1L || is.na(reference)) {
      arg::err("{.arg reference} must be a single response category")
    }

    reference <- as.character(reference)
  }

  new_genbart_family("multinomial", link, reference = reference)
}

#' @rdname genbart-families
#' @export
weibull_aft <- function() {
  new_genbart_family("aft", "weibull")
}

#' @rdname genbart-families
#' @export
loglogistic_aft <- function() {
  new_genbart_family("aft", "loglogistic")
}

#' @rdname genbart-families
#' @export
lognormal_aft <- function() {
  new_genbart_family("aft", "lognormal")
}

#' @rdname genbart-families
#' @export
location_scale <- function(link = "identity") {
  link <- arg::match_arg(link, "identity")

  new_genbart_family("location_scale", link)
}

#' @rdname genbart-families
#' @export
zi_poisson <- function(link = "log") {
  link <- arg::match_arg(link, "log")

  new_genbart_family("zip", link)
}

#' @rdname genbart-families
#' @export
zi_negbin <- function(link = "log", theta = NULL) {
  link <- arg::match_arg(link, "log")

  if (!is_null(theta)) {
    arg::arg_number(theta)
    arg::arg_gt(theta, 0)
  }

  new_genbart_family("zinb", link, theta = theta)
}

#' @rdname genbart-families
#' @export
ordbeta <- function(link = "logit", phi = NULL) {
  link <- arg::match_arg(link, "logit")

  if (!is_null(phi)) {
    arg::arg_number(phi)
    arg::arg_gt(phi, 0)
  }

  new_genbart_family("ordbeta", link, phi = phi)
}

#' @rdname genbart-families
#' @export
custom_family <- function(logdens, num_predictors = 1L, start = 0,
                          derivatives = NULL, name = "custom") {
  if (!is.function(logdens)) {
    arg::err("{.arg logdens} must be a function of the response and the
              additive predictors")
  }

  arg::arg_whole_number(num_predictors)
  arg::arg_gte(num_predictors, 1)
  arg::arg_numeric(start)
  arg::arg_string(name)

  if (!is_null(derivatives) && !is.function(derivatives)) {
    arg::err("{.arg derivatives} must be a function or {.code NULL}")
  }

  num_predictors <- as.integer(num_predictors)

  if (length(start) != 1L && length(start) != num_predictors) {
    arg::err("{.arg start} must have one value, or one per additive predictor
              ({num_predictors})")
  }

  new_genbart_family("custom", "identity", logdens = logdens,
                     num_predictors = num_predictors,
                     start = rep(start, length.out = num_predictors),
                     derivatives = derivatives, name = name)
}

genbart_family_names <- c("gaussian", "binomial", "poisson", "negbin", "Gamma",
                          "Gamma_shape", "ordinal", "multinomial",
                          "weibull_aft", "loglogistic_aft", "lognormal_aft",
                          "location_scale", "zi_poisson", "zi_negbin",
                          "ordbeta")

# The scale each engine family's additive predictor natively lives on. A link
# the engine does not carry is handled by composing the caller's inverse link
# with the link named here, which is why only families with a single mean and a
# conventional link can take an arbitrary one.
native_links <- c(gaussian = "identity", binomial = "logit", poisson = "log",
                  negbin = "log", Gamma = "log")

valid_links <- list(custom = "identity",
                    gaussian = "identity",
                    binomial = c("logit", "probit", "cloglog"),
                    poisson = "log",
                    negbin = "log",
                    Gamma = "log",
                    ordinal = c("logit", "probit", "cloglog"),
                    multinomial = "logit",
                    aft = c("weibull", "loglogistic", "lognormal"),
                    location_scale = "identity",
                    zip = "log",
                    zinb = "log",
                    ordbeta = "logit")

# Normalize whatever the user passed to `family` into a genbart family object.
# Accepts a string, a family-generating function, a stats::family object, or one
# of the genbart families above, mirroring how glm() resolves the argument.
as_genbart_family <- function(family) {
  if (is.character(family)) {
    arg::arg_string(family)
    arg::arg_element(family, genbart_family_names)
    family <- get(family, mode = "function", envir = asNamespace("genbart"))
  }

  if (is.function(family)) {
    family <- family()
  }

  if (!inherits(family, "family")) {
    arg::err("{.arg family} must be a family name, a family function, or a
              family object, such as {.code binomial(\"logit\")}")
  }

  # stats::Gamma() and the genbart Gamma_shape() both report "Gamma", so they
  # reach the same engine family and differ only in the options they carry.
  name <- family[["family"]]
  link <- family[["link"]]

  if (!name %in% names(valid_links)) {
    arg::err("family {.val {name}} is not supported by {.fn genbart}. Supported
              families are {.val {names(valid_links)}}")
  }

  allowed <- valid_links[[name]]
  custom_link <- NULL

  if (!link %in% allowed) {
    if (!name %in% names(native_links)) {
      arg::err("the {.val {link}} link is not supported for the {.val {name}}
                family. Supported links are {.val {allowed}}")
    }

    # A link the engine does not carry is honored by composing it onto the
    # scale the engine's family works on. That needs the inverse link, which a
    # stats::family object carries; a bare name is resolved through make.link().
    custom_link <- link_functions(family, link)
  }

  # A genbart family carries its own options -- a fixed dispersion, a reference
  # category, a supplied log density -- which have to survive; a stats::family
  # object carries link machinery that has already been read off above.
  extra <- {
    if (inherits(family, "genbart_family")) {
      family[!names(family) %in% c("family", "link")]
    }
    else list()
  }

  family <- do.call(new_genbart_family,
                    c(list(family = name, link = link), extra))

  if (!is_null(custom_link)) {
    family[["custom_link"]] <- custom_link
    family[["native_link"]] <- native_links[[name]]
  }

  family
}

# The inverse link and its derivative, taken from the family object when it has
# them -- which is the case for anything stats::binomial() and friends return,
# including a link built by the caller and passed as a "link-glm" object -- and
# resolved from the name otherwise.
link_functions <- function(family, link) {
  out <- family[c("linkfun", "linkinv", "mu.eta")]

  if (!is.function(out[["linkinv"]]) || !is.function(out[["linkfun"]])) {
    resolved <- try(stats::make.link(link), silent = TRUE)

    if (inherits(resolved, "try-error")) {
      arg::err("the {.val {link}} link is not one {.fn stats::make.link} knows.
                Supply it as a link object, as in
                {.code binomial(link = my_link)}, so that the inverse link
                comes with it")
    }

    out <- resolved[c("linkfun", "linkinv", "mu.eta")]
  }

  if (!is.function(out[["mu.eta"]])) {
    out[["mu.eta"]] <- NULL
  }

  out
}

# Build the map from the caller's additive predictor to the scale the engine's
# family works on, together with its derivative when the link brought one.
compose_link <- function(custom_link, native) {
  linkinv <- custom_link[["linkinv"]]
  mu_eta <- custom_link[["mu.eta"]]

  theta <- switch(native,
    identity = function(eta) linkinv(eta),
    log = function(eta) log(linkinv(eta)),
    logit = function(eta) stats::qlogis(linkinv(eta)))

  dtheta <- {
    if (is_null(mu_eta)) NULL
    else switch(native,
      identity = function(eta) mu_eta(eta),
      log = function(eta) mu_eta(eta) / linkinv(eta),
      logit = function(eta) {
        mu <- linkinv(eta)
        mu_eta(eta) / (mu * (1 - mu))
      })
  }

  list(link_theta = theta, link_dtheta = dtheta)
}
