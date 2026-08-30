# Turn the response as it comes out of model.frame() into the numeric form the
# engine expects, and work out the family-specific starting values.
#
# Two quantities need care and are computed here rather than in the sampler.
#
# The offset holds an intercept-only fit on the link scale -- exact for most
# families, a moment approximation where censoring or a mixture makes the closed
# form awkward. The
# forests carry no intercept of their own, so anchoring the additive predictor
# at the null fit is what makes the leaf prior a statement about departures from
# that fit rather than about the absolute level of the response.
#
# eta_scale is the natural scale of each additive predictor, used to set the
# leaf prior. It is one for the families whose link already puts the predictor
# on a standardized scale, and a spread of the response for those, such as the
# Gaussian and the accelerated failure time families, whose predictor lives on
# the scale of the data.
prepare_response <- function(family, y, weights, offset, x, n) {

  name <- family[["family"]]
  link <- family[["link"]]

  out <- list(family = name, link = link, n_forest = 1L, n_aux = 0L,
              num_cat = NULL, levels = NULL, opts = list())

  if (is_null(weights)) {
    weights <- rep.int(1, n)
  }
  else {
    arg::arg_numeric(weights)

    if (length(weights) != n) {
      arg::err("{.arg weights} must have one value per observation")
    }

    if (any(weights < 0)) {
      arg::err("{.arg weights} must be non-negative")
    }
  }

  switch(name,
    gaussian = {
      y <- check_numeric_response(y, name)
      out$y <- y
      out$weights <- weights
      out$eta_scale <- stats::sd(y)
      out$opts <- list(sigma_hat = residual_scale(y, x, weights))
      intercept <- stats::weighted.mean(y, weights)
    },

    # A Dirichlet process mixture for the error distribution. The predictor's
    # own prior is the Gaussian family's, since conditional on the mixture the
    # target is the same Gaussian one; what is different is the baseline
    # distribution of the mixture, which is set up here because it is calibrated
    # from a linear fit.
    dpm_aft = {
      if (!isTRUE(all.equal(unname(weights), rep.int(1, n)))) {
        arg::err("{.fn dpm_aft} does not take prior weights, because a weight
                  would have to be a multiplicity in the Dirichlet process,
                  which is not what a fractional weight means",
                 i = "the other survival families take them")
      }

      a <- prepare_surv(y, n)
      out$y <- a$log_time
      out$weights <- weights
      out$eta_scale <- stats::sd(a$log_time)

      # The mixture's baseline is calibrated the way `dpm()` calibrates it, from
      # a linear fit -- but on the log times of the *observed events*, since a
      # censoring time understates its own failure time and would pull the
      # residual scale down.
      seen <- a$event > 0
      out$opts <- c(dpm_baseline(family, a$log_time[seen], x[seen, , drop = FALSE],
                                 sum(seen)),
                    list(event = a$event))
      intercept <- mean(a$log_time[seen])
    },

    dpm = {
      if (!isTRUE(all.equal(unname(weights), rep.int(1, n)))) {
        arg::err("{.fn dpm} does not take prior weights, because a weight
                  would have to be a multiplicity in the Dirichlet process,
                  which is not what a fractional weight means",
                 i = "{.fn gaussian}, {.fn ordinal} and {.fn location_scale} all
                      take them")
      }

      y <- check_numeric_response(y, name)
      out$y <- y
      out$weights <- weights
      out$eta_scale <- stats::sd(y)
      out$opts <- dpm_baseline(family, y, x, n)
      intercept <- stats::weighted.mean(y, weights)
    },

    # A family supplied from R. Nothing is known about the log density beyond
    # the function itself, so the starting value comes from the caller and the
    # leaf prior scale is the default one for a predictor on a standardized
    # scale.
    custom = {
      y <- check_numeric_response(y, name)
      out$y <- y
      out$weights <- weights

      # A nuisance parameter is carried as a trailing additive predictor whose
      # forest the engine pins at depth zero, so it is one scalar drawn by the
      # ordinary leaf machinery. `n_forest` counts them because the engine has to
      # build them; `n_aux` is how many of those trailing forests are nuisances,
      # and they are reported in `aux` rather than in `eta`.
      n_aux <- family[["num_aux"]] %or% 0L
      out$n_aux <- n_aux
      out$n_forest <- family[["num_predictors"]] + n_aux
      out$eta_scale <- rep.int(1, out$n_forest)
      out$opts <- list(num_predictors = family[["num_predictors"]],
                       num_aux = n_aux,
                       aux_names = as.character(family[["aux_names"]]),
                       logdens = family[["logdens"]],
                       derivatives = family[["derivatives"]],
                       name = family[["name"]])
      intercept <- {
        if (n_aux > 0L) c(family[["start"]], family[["aux_start"]])
        else family[["start"]]
      }
    },

    location_scale = {
      y <- check_numeric_response(y, name)
      out$y <- y
      out$weights <- weights
      out$n_forest <- 2L
      # The second predictor is a log standard deviation. A tighter prior scale
      # is used for it, because a forest prior as wide as the one used on the
      # mean would allow the standard deviation to vary by more than an order of
      # magnitude before seeing any data.
      out$eta_scale <- c(stats::sd(y), 0.5)
      intercept <- c(stats::weighted.mean(y, weights), log(stats::sd(y)))
    },

    binomial = {
      b <- prepare_binomial(y, weights, n)
      out$y <- b$y
      out$weights <- b$weights
      out$levels <- b$levels
      out$eta_scale <- 1
      p <- stats::weighted.mean(b$y, b$weights)
      p <- min(max(p, 1 / (2 * sum(b$weights))), 1 - 1 / (2 * sum(b$weights)))
      intercept <- binomial_link(p, link)
    },

    poisson = ,
    negbin = {
      y <- check_count_response(y, name)
      out$y <- y
      out$weights <- weights
      out$eta_scale <- 1
      mu <- max(stats::weighted.mean(y, weights), 1 / (2 * n))
      intercept <- log(mu)

      if (identical(name, "negbin")) {
        v <- stats::var(y)
        theta_start <- {
          if (v > mu) mu^2 / (v - mu)
          else 10
        }
        out$opts <- list(theta = family[["theta"]] %or% theta_start,
                         theta_prior_shape = 0.01,
                         theta_prior_rate = 0.01,
                         update_theta = is_null(family[["theta"]]))
      }
    },

    Gamma = {
      y <- check_numeric_response(y, name)

      if (any(y <= 0)) {
        arg::err("the {.val Gamma} family requires a strictly positive response")
      }

      out$y <- y
      out$weights <- weights
      out$eta_scale <- 1
      mu <- stats::weighted.mean(y, weights)
      intercept <- log(mu)
      shape_start <- max(mu^2 / stats::var(y), 0.1)
      # The engine can hold the shape fixed, and nothing exposes that: a caller
      # who knows the gamma shape is rare enough that the argument was not worth
      # a second family function. `update_shape` stays here because the engine
      # reads it.
      out$opts <- list(shape = shape_start,
                       shape_prior_shape = 0.01,
                       shape_prior_rate = 0.01,
                       update_shape = TRUE)
    },

    ordinal = {
      o <- prepare_ordered(y, name)
      out$y <- o$codes
      out$weights <- weights
      out$levels <- o$levels
      out$num_cat <- o$num_cat
      out$eta_scale <- 1

      # Place the cutpoints at the observed cumulative proportions on the link
      # scale, then shift so the first sits at zero and let the offset carry the
      # location. This makes a two-category response identical to binary
      # regression with the same link.
      props <- cumsum(tabulate(o$codes + 1L, nbins = o$num_cat)) / length(o$codes)
      props <- pmin(pmax(props[-o$num_cat], 1 / (2 * n)), 1 - 1 / (2 * n))
      raw <- binomial_link(props, link)
      out$opts <- list(num_cat = o$num_cat,
                       cuts = raw - raw[1L],
                       update_cuts = o$num_cat > 2L)
      intercept <- -raw[1L]
    },

    multinomial = {
      m <- prepare_unordered(y, name, family[["reference"]])
      symmetric <- is_null(family[["reference"]])
      out$y <- m$codes
      out$weights <- weights
      out$levels <- m$levels
      out$num_cat <- m$num_cat
      out$n_forest <- m$num_cat - as.integer(!symmetric)
      out$opts <- list(num_cat = m$num_cat, symmetric = symmetric)

      counts <- pmax(tabulate(m$codes + 1L, nbins = m$num_cat), 0.5)

      # Under the symmetric coding the two log-odds contrasts that share a
      # category are each a difference of two forests, so their prior variance
      # is twice what a single forest carries. Shrinking the leaf scale by
      # sqrt(2) leaves the prior on the identified log odds where reference
      # coding puts it (Murray 2021, sec. 4.3).
      if (symmetric) {
        out$eta_scale <- rep.int(1 / sqrt(2), m$num_cat)
        # Only contrasts are identified, so the levels are centered rather than
        # taken against a reference.
        intercept <- log(counts) - mean(log(counts))
      }
      else {
        out$eta_scale <- rep.int(1, m$num_cat - 1L)
        intercept <- log(counts[-1L] / counts[1L])
      }
    },

    # Multinomial probit. Like the reference-coded multinomial in its coding --
    # one forest per non-reference category -- and unlike it in having a
    # covariance matrix, which is what the probit link is for.
    mnp = {
      m <- prepare_unordered(y, name, family[["reference"]] %or%
                               first_level(y))
      out$y <- m$codes
      out$weights <- weights
      out$levels <- m$levels
      out$num_cat <- m$num_cat
      out$n_forest <- m$num_cat - 1L
      out$eta_scale <- rep.int(1, m$num_cat - 1L)

      # Imai and van Dyk's (2005) choice, nu = C + 1 with Psi the identity,
      # which puts a uniform prior on the correlations of the unnormalized
      # covariance matrix. Not exposed: see the family's documentation for what
      # raising it actually does.
      out$opts <- list(num_cat = m$num_cat,
                       nu = m$num_cat,
                       update_sigma = m$num_cat > 2L,
                       replicates = family[["replicates"]])

      counts <- pmax(tabulate(m$codes + 1L, nbins = m$num_cat), 0.5)
      shares <- counts / sum(counts)

      # The pairwise probit contrast of each category against the reference,
      # which is the exact intercept when there are two categories and a
      # reasonable anchor when there are more.
      intercept <- stats::qnorm(shares[-1L] / (shares[-1L] + shares[1L]))
    },

    zip = ,
    zinb = {
      y <- check_count_response(y, name)
      out$y <- y
      out$weights <- weights
      out$n_forest <- 2L
      out$eta_scale <- c(1, 1)

      # Split the observed zeros into the share a plain count model would
      # produce and the excess, which is what the inflation component absorbs.
      p_zero <- mean(y == 0)
      mu_start <- {
        if (any(y > 0)) mean(y[y > 0])
        else max(mean(y), 0.1)
      }
      excess <- min(max(p_zero - exp(-mu_start), 0.01), 0.9)
      intercept <- c(log(mu_start), stats::qlogis(excess))

      if (identical(name, "zinb")) {
        v <- stats::var(y)
        m <- max(mean(y), 1 / (2 * n))
        theta_start <- {
          if (v > m) m^2 / (v - m)
          else 10
        }
        out$opts <- list(theta = family[["theta"]] %or% theta_start,
                         theta_prior_shape = 0.01,
                         theta_prior_rate = 0.01,
                         update_theta = is_null(family[["theta"]]))
      }
    },

    beta = {
      y <- check_numeric_response(y, name)

      if (any(y <= 0) || any(y >= 1)) {
        arg::err("The {.val beta} family requires a response strictly between 0
                  and 1.",
                 i = "With observations at 0 or 1, use {.fn ordbeta}, which
                      models those as point masses.")
      }

      out$y <- y
      out$weights <- weights
      out$eta_scale <- 1

      mid <- mean(y)
      intercept <- stats::qlogis(mid)

      # Method of moments for the precision, which is what a beta regression
      # would start from: var = mu (1 - mu) / (1 + phi).
      v <- stats::var(y)
      phi_start <- {
        if (length(y) > 1L && v > 0 && v < mid * (1 - mid)) {
          mid * (1 - mid) / v - 1
        }
        else 5
      }

      out$opts <- list(phi = family[["phi"]] %or% phi_start,
                       phi_prior_shape = 0.01,
                       phi_prior_rate = 0.01,
                       update_phi = is_null(family[["phi"]]))
    },

    ordbeta = {
      y <- check_numeric_response(y, name)

      if (any(y < 0) || any(y > 1)) {
        arg::err("the {.val ordbeta} family requires a response between 0 and 1")
      }

      interior <- y > 0 & y < 1

      if (!any(interior)) {
        arg::err("the {.val ordbeta} family needs some responses strictly
                  between 0 and 1; with only 0 and 1 use {.fn binomial}")
      }

      out$y <- y
      out$weights <- weights
      out$eta_scale <- 1

      # Cutpoints are placed at the observed proportions of zeros and ones, on
      # the logit scale and relative to the intercept, so the starting values
      # reproduce the observed endpoint masses.
      floor_p <- 1 / (2 * n)
      p_zero <- min(max(mean(y == 0), floor_p), 1 - floor_p)
      p_one <- min(max(mean(y == 1), floor_p), 1 - floor_p)
      mid <- min(max(mean(y[interior]), floor_p), 1 - floor_p)
      intercept <- stats::qlogis(mid)

      cut1 <- intercept - stats::qlogis(1 - p_zero)
      cut2 <- intercept - stats::qlogis(p_one)

      if (cut1 >= cut2) {
        spread <- max(abs(intercept), 1)
        cut1 <- intercept - spread
        cut2 <- intercept + spread
      }

      # Method of moments for the beta precision on the interior responses.
      v <- stats::var(y[interior])
      phi_start <- {
        if (sum(interior) > 1L && v > 0 && v < mid * (1 - mid)) {
          mid * (1 - mid) / v - 1
        }
        else 5
      }

      out$opts <- list(cut1 = cut1,
                       cut2 = cut2,
                       phi = family[["phi"]] %or% max(phi_start, 0.5),
                       phi_prior_shape = 0.01,
                       phi_prior_rate = 0.01,
                       update_phi = is_null(family[["phi"]]))
    },

    ph = {
      a <- prepare_surv(y, n)
      out$y <- a$time
      out$weights <- weights

      if (!any(a$event > 0)) {
        arg::err("{.fn ph} needs at least one observed event")
      }

      # Bin edges at evenly spaced quantiles of the observed times, and about
      # n^(1/3) of them, which is the order the Freedman-Diaconis rule gives for
      # a histogram and what Basak et al. (2024) recommend. The first edge is
      # zero; the last bin runs to infinity.
      num_bins <- family[["num_bins"]] %or%
        max(2L, min(30L, as.integer(ceiling(n^(1 / 3)))))
      probs <- seq(0, 1, length.out = num_bins + 1L)
      edges <- unique(c(0, unname(stats::quantile(a$time,
                                                 probs[-c(1L, num_bins + 1L)]))))

      # The predictor is a log hazard ratio, identified only against the baseline,
      # so it starts at zero and the baseline carries the level.
      out$eta_scale <- 1
      hazard <- sum(a$event) / max(sum(a$time), .Machine$double.eps)
      out$opts <- list(event = a$event,
                       edges = edges,
                       lambda_shape = family[["lambda_shape"]] %or% 1,
                       lambda_rate = 1 / hazard,
                       update_lambda = family[["update_lambda"]] %or% TRUE)
      intercept <- 0
    },

    aft = {
      a <- prepare_surv(y, n)
      out$y <- a$log_time
      out$weights <- weights
      out$eta_scale <- stats::sd(a$log_time)
      out$opts <- list(event = a$event,
                       sigma_hat = stats::sd(a$log_time),
                       update_sigma = TRUE)
      intercept <- stats::weighted.mean(a$log_time, weights)
    },

    arg::err("family {.val {name}} is not supported")
  )

  if (!all(is.finite(out$eta_scale)) || any(out$eta_scale <= 0)) {
    arg::err("the response has no variation, so there is nothing to model")
  }

  # The branches above put the intercept on the scale the engine's family works
  # on. When the caller supplied a link the engine does not carry natively, the
  # additive predictor lives on the caller's scale instead, so the same
  # intercept-only fit has to be re-expressed there, and the engine is told to
  # use its native link with the composition in front of it.
  if (!is_null(family[["custom_link"]])) {
    native <- family[["native_link"]]
    mu <- switch(native,
                 identity = intercept,
                 log = exp(intercept),
                 logit = stats::plogis(intercept))
    intercept <- family[["custom_link"]][["linkfun"]](mu)

    if (!all(is.finite(intercept))) {
      arg::err("the {.val {link}} link maps the intercept-only fit to a
                non-finite value; supply {.arg offset} to start the additive
                predictor somewhere valid")
    }

    out$link <- native
    out$opts <- c(out$opts, compose_link(family[["custom_link"]], native))
  }

  out$offset <- build_offset(intercept, offset, out$n_forest, n)

  out
}

# The offset is stored with one row per additive predictor, matching the layout
# the engine uses.
build_offset <- function(intercept, offset, n_forest, n) {
  out <- matrix(rep(intercept, length.out = n_forest), nrow = n_forest,
                ncol = n)

  if (is_null(offset)) {
    return(out)
  }

  arg::arg_numeric(offset)

  if (is.matrix(offset)) {
    if (nrow(offset) != n || ncol(offset) != n_forest) {
      arg::err("{.arg offset} must be a matrix with {n} rows and {n_forest}
                columns for this family")
    }
    return(out + t(offset))
  }

  if (length(offset) != n) {
    arg::err("{.arg offset} must have one value per observation")
  }

  out + matrix(offset, nrow = n_forest, ncol = n, byrow = TRUE)
}

check_numeric_response <- function(y, name) {
  if (is.matrix(y) && ncol(y) == 1L) {
    y <- drop(y)
  }

  if (!is.numeric(y)) {
    arg::err("the {.val {name}} family requires a numeric response")
  }

  if (anyNA(y)) {
    arg::err("the response contains missing values")
  }

  as.numeric(y)
}

check_count_response <- function(y, name) {
  y <- check_numeric_response(y, name)

  if (any(y < 0) || any(abs(y - round(y)) > 1e-8)) {
    arg::err("the {.val {name}} family requires a response of non-negative
              counts")
  }

  round(y)
}

binomial_link <- function(p, link) {
  switch(link,
    logit = stats::qlogis(p),
    probit = stats::qnorm(p),
    cloglog = log(-log1p(-p)),
    stats::qlogis(p))
}

binomial_linkinv <- function(eta, link) {
  switch(link,
    logit = stats::plogis(eta),
    probit = stats::pnorm(eta),
    cloglog = -expm1(-exp(eta)),
    stats::plogis(eta))
}

# Accepts the same shapes of binomial response that glm() does, and returns the
# response as a proportion with the number of trials in the weights.
prepare_binomial <- function(y, weights, n) {
  levels <- NULL

  if (is.matrix(y) && ncol(y) == 2L) {
    trials <- rowSums(y)

    if (any(trials <= 0)) {
      arg::err("every row of a two-column binomial response must have at least
                one trial")
    }

    return(list(y = y[, 1L] / trials, weights = weights * trials,
                levels = colnames(y)))
  }

  if (is.matrix(y) && ncol(y) == 1L) {
    y <- drop(y)
  }

  if (is.factor(y)) {
    levels <- levels(y)

    if (length(levels) != 2L) {
      arg::err("a factor response for the {.val binomial} family must have
                exactly two levels; use {.fn multinomial} or {.fn ordinal} for
                more")
    }

    y <- as.numeric(y != levels[1L])
  }
  else if (is.logical(y)) {
    y <- as.numeric(y)
  }
  else {
    y <- check_numeric_response(y, "binomial")

    if (any(y < 0) || any(y > 1)) {
      arg::err("a numeric response for the {.val binomial} family must lie
                between 0 and 1")
    }
  }

  list(y = y, weights = weights, levels = levels)
}

prepare_ordered <- function(y, name) {
  if (is.ordered(y)) {
    levels <- levels(y)
    codes <- as.integer(y) - 1L
  }
  else if (is.factor(y)) {
    arg::wrn(c("The response is an unordered factor.",
               i = "The {.val ordinal} family will use the existing level
                    order: {.val {levels(y)}}",
               i = "Supply an {.cls ordered} factor to make the order
                    explicit."))
    levels <- levels(y)
    codes <- as.integer(y) - 1L
  }
  else {
    values <- sort(unique(y))
    levels <- as.character(values)
    codes <- match(y, values) - 1L
  }

  num_cat <- length(levels)

  if (num_cat < 2L) {
    arg::err("the response has only one category")
  }

  list(codes = as.numeric(codes), levels = levels, num_cat = num_cat)
}

# The first level of a response, used as the reference category when the caller
# does not name one. A multinomial probit is written as contrasts against a
# reference, so unlike the symmetric multinomial it always has one.
first_level <- function(y) {
  if (is.factor(y)) {
    return(levels(y)[1L])
  }

  as.character(sort(unique(y))[1L])
}

prepare_unordered <- function(y, name, reference = NULL) {
  if (is.factor(y)) {
    levels <- levels(y)
    codes <- as.integer(y) - 1L
  }
  else {
    values <- sort(unique(y))
    levels <- as.character(values)
    codes <- match(y, values) - 1L
  }

  num_cat <- length(levels)

  if (num_cat < 2L) {
    arg::err("the response has only one category")
  }

  # Drop unused levels rather than fitting a forest that can never fire.
  present <- tabulate(codes + 1L, nbins = num_cat) > 0

  if (!all(present)) {
    arg::wrn("dropping {sum(!present)} unused response level{?s}:
              {.val {levels[!present]}}")
    levels <- levels[present]
    codes <- match(codes, which(present) - 1L) - 1L
    num_cat <- length(levels)
  }

  # A reference category is put first, since the engine always drops the first
  # predictor when it is reference-coded.
  if (!is_null(reference)) {
    at <- match(reference, levels)

    if (is.na(at)) {
      arg::err("{.arg reference} must name one of the response categories
                {.val {levels}}")
    }

    order <- c(at, setdiff(seq_along(levels), at))
    levels <- levels[order]
    codes <- match(codes + 1L, order) - 1L
  }

  list(codes = as.numeric(codes), levels = levels, num_cat = num_cat)
}

prepare_surv <- function(y, n) {
  if (inherits(y, "Surv")) {
    type <- attr(y, "type")

    if (!identical(type, "right")) {
      arg::err("the accelerated failure time families support right-censored
                data only, but the response has type {.val {type}}")
    }

    y <- unclass(y)
  }

  if (!is.matrix(y) || ncol(y) != 2L) {
    arg::err("an accelerated failure time family needs a two-column response of
              times and event indicators, most easily supplied as
              {.fn survival::Surv}")
  }

  time <- as.numeric(y[, 1L])
  event <- as.numeric(y[, 2L])

  if (anyNA(time) || anyNA(event)) {
    arg::err("the response contains missing values")
  }

  if (any(time <= 0)) {
    arg::err("survival times must be strictly positive")
  }

  if (!all(event %in% c(0, 1))) {
    arg::err("the event indicator must be 0 for censored and 1 for observed
              events")
  }

  if (!any(event == 1)) {
    arg::err("every observation is censored, so the model is not identified")
  }

  list(time = time, log_time = log(time), event = event)
}

# The baseline distribution of the Dirichlet process mixture, and the prior on
# its concentration. Everything here is calibrated from a linear fit, following
# George et al. (2019, secs. 3.1 and 3.2), which is the same device BART uses for
# its own scale prior and is why the defaults need no tuning.
dpm_baseline <- function(family, y, x, n) {
  residuals <- linear_residuals(y, x)
  sigma_hat <- residual_scale(y, x)

  nu <- family[["nu"]]
  q <- family[["q"]]

  # lambda placed so that P(sigma < sigma_hat) = q under
  # sigma^2 ~ nu lambda / chisq_nu, which is BART's construction with the
  # paper's larger nu and higher quantile: the mixture covers small errors with
  # extra components, so a single component's prior can afford to be tighter.
  lambda <- sigma_hat^2 * stats::qchisq(1 - q, nu) / nu

  # k_0 scales the baseline's mean so that the marginal of mu, which is
  # sqrt(lambda / k_0) times a t on nu degrees of freedom, reaches the edge of
  # the residuals at `k_s` of its own scale units.
  reach <- max(abs(residuals))
  reach <- if (reach > 0) reach else 1
  k_0 <- lambda * family[["k_s"]]^2 / reach^2

  # The concentration's prior is Rossi's: pick the smallest and largest cluster
  # counts thought plausible, turn each into a concentration, and taper between
  # them. `alpha log(1 + n / alpha)` is the expected number of occupied clusters,
  # which is the standard approximation and is what makes this solvable in one
  # line each.
  most <- family[["max_clusters"]] %or% max(as.integer(0.1 * n), 2L)
  most <- min(most, n)
  alpha_min <- concentration_for(1, n)
  alpha_max <- concentration_for(most, n)

  if (!(alpha_max > alpha_min)) {
    alpha_max <- alpha_min * 10
  }

  grid <- seq(alpha_min, alpha_max, length.out = 100L)
  taper <- 1 - (grid - alpha_min) / (alpha_max - alpha_min)
  logprior <- family[["psi"]] * log(pmax(taper, 1e-12))

  list(nu = nu,
       lambda = lambda,
       mu_0 = 0,
       k_0 = k_0,
       alpha = family[["alpha"]] %or% stats::median(grid),
       update_alpha = is_null(family[["alpha"]]),
       alpha_grid = grid,
       alpha_logprior = logprior)
}

# The concentration that gives `clusters` occupied clusters on average, from
# `clusters = alpha * log(1 + n / alpha)`, solved numerically because it has no
# closed form.
concentration_for <- function(clusters, n) {
  expected <- function(alpha) alpha * log1p(n / alpha) - clusters

  if (expected(1e-8) > 0) {
    return(1e-8)
  }

  stats::uniroot(expected, c(1e-8, 10 * n), tol = 1e-8)$root
}

# Residuals of the linear fit the baseline is calibrated from, which is the same
# fit residual_scale() reads its own anchor off.
linear_residuals <- function(y, x) {
  if (is_null(x) || nrow(x) <= ncol(x) + 1L) {
    return(y - mean(y))
  }

  fit <- try(stats::lm.fit(cbind(1, x), y), silent = TRUE)

  if (inherits(fit, "try-error")) {
    return(y - mean(y))
  }

  fit[["residuals"]]
}

# Prior scale for the residual standard deviation. A linear fit gives a much
# better anchor than the marginal spread when the predictors explain anything,
# and this is the same device the original BART implementations use.
#
# The weights enter because they say how much each observation is worth, and a
# prior scale read off the wrong observations is the wrong scale. Only their
# relative sizes can matter: a residual variance is per observation, so weights
# that do not average one over the rows they keep would rescale it. Normalizing
# over the *kept* rows leaves unit weights exactly where they were and makes a
# zero weight mean "drop this row" rather than "shrink every scale".
residual_scale <- function(y, x, weights = NULL) {
  w <- weights %or% rep.int(1, length(y))
  contributes <- w > 0
  kept <- sum(contributes)

  if (kept < 2L) {
    return(stats::sd(y))
  }

  average <- mean(w[contributes])

  if (!isTRUE(average > 0)) {
    return(stats::sd(y))
  }

  w <- w / average
  s <- sqrt(sum(w * (y - stats::weighted.mean(y, w))^2) / (kept - 1L))

  if (!is.finite(s) || s <= 0) {
    return(stats::sd(y))
  }

  if (is_null(x) || kept <= ncol(x) + 1L) {
    return(s)
  }

  fit <- try(stats::lm.wfit(cbind(1, x), y, w), silent = TRUE)

  if (inherits(fit, "try-error")) {
    return(s)
  }

  df <- fit[["df.residual"]]

  if (is_null(df) || df < 1L) {
    return(s)
  }

  # `lm.wfit()` returns residuals on the response's own scale, so the weights go
  # back in here; this is what `summary.lm()` does for a weighted fit.
  out <- sqrt(sum(w * fit[["residuals"]]^2) / df)

  if (!is.finite(out) || out <= 0) s else min(out, s)
}
