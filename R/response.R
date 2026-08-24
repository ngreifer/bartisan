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

  out <- list(family = name, link = link, n_forest = 1L, num_cat = NULL,
              levels = NULL, opts = list())

  if (is_null(weights)) {
    weights <- rep(1, n)
  }

  arg::arg_numeric(weights)

  if (length(weights) != n) {
    arg::err("{.arg weights} must have one value per observation")
  }

  if (any(weights < 0)) {
    arg::err("{.arg weights} must be non-negative")
  }

  switch(name,
    gaussian = {
      y <- check_numeric_response(y, name)
      out$y <- y
      out$weights <- weights
      out$eta_scale <- stats::sd(y)
      out$opts <- list(sigma_hat = residual_scale(y, x))
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
      out$n_forest <- family[["num_predictors"]]
      out$eta_scale <- rep(1, out$n_forest)
      out$opts <- list(num_predictors = out$n_forest,
                       logdens = family[["logdens"]],
                       derivatives = family[["derivatives"]],
                       name = family[["name"]])
      intercept <- family[["start"]]
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
      out$opts <- list(shape = family[["shape"]] %or% shape_start,
                       shape_prior_shape = 0.01,
                       shape_prior_rate = 0.01,
                       update_shape = is_null(family[["shape"]]))
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
        out$eta_scale <- rep(1 / sqrt(2), m$num_cat)
        # Only contrasts are identified, so the levels are centered rather than
        # taken against a reference.
        intercept <- log(counts) - mean(log(counts))
      }
      else {
        out$eta_scale <- rep(1, m$num_cat - 1L)
        intercept <- log(counts[-1L] / counts[1L])
      }
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

      if (!(cut1 < cut2)) {
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
    cli::cli_warn(c("the response is an unordered factor",
                    i = "the {.val ordinal} family will use the existing level
                         order: {.val {levels(y)}}",
                    i = "supply an {.cls ordered} factor to make the order
                         explicit"))
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
    cli::cli_warn("dropping {sum(!present)} unused response level{?s}:
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

  list(log_time = log(time), event = event)
}

# Prior scale for the residual standard deviation. A linear fit gives a much
# better anchor than the marginal spread when the predictors explain anything,
# and this is the same device the original BART implementations use.
residual_scale <- function(y, x) {
  s <- stats::sd(y)

  if (is_null(x) || nrow(x) <= ncol(x) + 1L) {
    return(s)
  }

  fit <- try(stats::lm.fit(cbind(1, x), y), silent = TRUE)

  if (inherits(fit, "try-error")) {
    return(s)
  }

  df <- fit[["df.residual"]]

  if (is_null(df) || df < 1L) {
    return(s)
  }

  out <- sqrt(sum(fit[["residuals"]]^2) / df)

  if (!is.finite(out) || out <= 0) s else min(out, s)
}
