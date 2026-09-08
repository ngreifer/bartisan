# The Tweedie compound Poisson-gamma. What makes it worth its own file is that
# its density has no closed form at a positive response, so the thing to test is
# not the plumbing but the series: the sampler is handed a log density it can
# only get right by summing one, and the reference below sums a different one.

# The compound Poisson-gamma written straight from its representation, sharing no
# code with the engine: a Poisson count of gamma draws, summed over the count.
# This is the independent answer the series has to reproduce.
tweedie_logdens_ref <- function(y, mu, phi, power, nmax = 60000L) {
  lambda <- mu^(2 - power) / (phi * (2 - power))
  shape <- (2 - power) / (power - 1)
  scale <- phi * (power - 1) * mu^(power - 1)

  if (y == 0) {
    return(-lambda)
  }

  n <- seq_len(nmax)
  terms <- -lambda + n * log(lambda) - lgamma(n + 1) +
    (n * shape - 1) * log(y) - y / scale - n * shape * log(scale) -
    lgamma(n * shape)
  top <- max(terms)
  top + log(sum(exp(terms - top)))
}

tweedie_opts <- function(phi, power) {
  list(phi = phi, power = power, phi_prior_shape = 0.01,
       phi_prior_rate = 0.01, update_phi = FALSE, update_power = FALSE)
}

tweedie_engine_logdens <- function(y, eta, phi, power) {
  as.vector(.bartisan_logdens(y, rep(1, length(y)), list(matrix(eta, 1L)),
                              "tweedie", "log", tweedie_opts(phi, power),
                              matrix(c(phi, power), 1L, 2L)))
}

# A draw from the family, for the recovery test and for comparing against
# `simulate()`.
rtweedie_ref <- function(n, mu, phi, power) {
  lambda <- mu^(2 - power) / (phi * (2 - power))
  shape <- (2 - power) / (power - 1)
  scale <- phi * (power - 1) * mu^(power - 1)
  count <- stats::rpois(n, lambda)
  out <- numeric(n)
  hit <- count > 0L
  out[hit] <- stats::rgamma(sum(hit), shape = count[hit] * shape,
                            scale = scale[hit])
  out
}

test_that("the series reproduces the compound Poisson-gamma sum", {
  # Spread over the power, the dispersion and the response, because the number
  # of terms the series needs depends on all three and on none of them the same
  # way. The last cell is earnings-scale, which is what the family is for.
  cases <- list(
    list(y = c(0, 0.3, 1.4, 6, 25), eta = c(0.5, 0.5, 1.0, 1.5, 2.0),
         phi = 1.3, power = 1.5),
    list(y = c(0, 0.01, 2, 900), eta = c(-1, 0, 3, 6.5),
         phi = 0.4, power = 1.2),
    list(y = c(0, 5, 50, 5000), eta = c(1, 2, 4, 8), phi = 20, power = 1.8),
    list(y = c(0, 1e3, 6e4), eta = c(8, 8.5, 9), phi = 180, power = 1.433),
    list(y = c(0, 7216), eta = c(8.8, 8.8), phi = 10, power = 1.433))

  for (case in cases) {
    got <- tweedie_engine_logdens(case[["y"]], case[["eta"]], case[["phi"]],
                                  case[["power"]])
    want <- vapply(seq_along(case[["y"]]), function(i) {
      tweedie_logdens_ref(case[["y"]][i], exp(case[["eta"]][i]), case[["phi"]],
                          case[["power"]])
    }, numeric(1L))

    expect_equal(got, want, tolerance = 1e-8,
                 info = sprintf("phi = %g, power = %g", case[["phi"]],
                                case[["power"]]))
  }
})

test_that("the zero is the closed form rather than a summed one", {
  # P(y = 0) = exp(-mu^(2-p) / (phi (2-p))) exactly, and it is the whole of the
  # log density there, so the series must contribute nothing.
  for (power in c(1.2, 1.5, 1.8)) {
    for (phi in c(0.5, 4, 180)) {
      mu <- c(1, 10, 6800)
      got <- tweedie_engine_logdens(rep(0, 3L), log(mu), phi, power)
      expect_equal(got, -mu^(2 - power) / (phi * (2 - power)),
                   tolerance = 1e-12)
    }
  }
})

test_that("the score and information are the closed forms, and not negative", {
  for (g in list(c(1.3, 1.5), c(0.4, 1.2), c(20, 1.8), c(180, 1.433))) {
    phi <- g[1L]
    power <- g[2L]
    y <- c(0, 0.4, 3.2, 40, 900)
    eta <- matrix(c(0.2, 1.1, 2.4, 3.9, 6.2), 1L)
    mu <- exp(as.vector(eta))

    got <- .bartisan_derivs(y, rep(1, length(y)), list(eta), "tweedie", "log",
                            tweedie_opts(phi, power),
                            matrix(c(phi, power), 1L, 2L), 0L, FALSE, FALSE)

    expect_equal(as.vector(got[["d1"]]), mu^(1 - power) * (y - mu) / phi,
                 tolerance = 1e-10)

    # The observed second derivative, not the expected information: both of its
    # terms are positive for a power in (1, 2) and a non-negative response, so
    # the true curvature is usable as it stands.
    expect_equal(as.vector(got[["info"]]),
                 ((2 - power) * mu^(2 - power) +
                    (power - 1) * y * mu^(1 - power)) / phi,
                 tolerance = 1e-10)

    # And against the engine's own differencing, which is the check that the
    # analytic forms belong to the log density above rather than to something
    # else.
    diffed <- .bartisan_derivs(y, rep(1, length(y)), list(eta), "tweedie",
                               "log", tweedie_opts(phi, power),
                               matrix(c(phi, power), 1L, 2L), 0L, TRUE, FALSE)
    expect_equal(got[["d1"]], diffed[["d1"]], tolerance = 1e-6)
    expect_equal(got[["info"]], diffed[["info"]], tolerance = 1e-5)
  }

  set.seed(77L)
  y <- c(0, stats::runif(200L, 0, 500))
  eta <- matrix(stats::runif(201L, -3, 7), 1L)
  got <- .bartisan_derivs(y, rep(1, length(y)), list(eta), "tweedie", "log",
                          tweedie_opts(4, 1.5), matrix(c(4, 1.5), 1L, 2L), 0L,
                          FALSE, FALSE)
  expect_true(all(got[["info"]] >= 0))
})

test_that("the power is fixed by default and drawn when asked for", {
  expect_identical(tweedie()[["power"]], 1.5)
  expect_null(tweedie()[["phi"]])
  expect_null(tweedie(power = NULL)[["power"]])
  expect_identical(tweedie(power = 1.8, phi = 2)[["phi"]], 2)

  # The open interval only, since the family degenerates at either end. Matched
  # on the argument name rather than on the bounds, which are `arg`'s to word
  # and have been reworded once already; what this needs to pin down is the
  # contract, so the interior is asserted too.
  expect_error(tweedie(power = 1), "power")
  expect_error(tweedie(power = 2), "power")
  expect_error(tweedie(power = 0.5), "power")
  expect_error(tweedie(power = 3), "power")
  expect_identical(tweedie(power = 1.001)[["power"]], 1.001)
  expect_identical(tweedie(power = 1.999)[["power"]], 1.999)
  expect_error(tweedie(power = c(1.4, 1.6)), "single number")
  expect_error(tweedie(phi = 0), "must be positive")
  expect_error(tweedie(link = "identity"), "log")
})

test_that("naming the distribution any of three ways fits the same model", {
  # An ordinary `family` object carries none of our settings, so an empty option
  # list has to read as "take the defaults" rather than as "draw everything".
  # Otherwise `family = "tweedie"` and a `family` object naming a Tweedie would
  # be the same distribution and a different model.
  d <- sim_x(n = 200, seed = 71)
  set.seed(171)
  d$y <- rtweedie_ref(nrow(d), exp(1.5 + d$x1), 4, 1.5)

  plain <- structure(list(family = "tweedie", link = "log"), class = "family")

  fits <- lapply(list("tweedie", tweedie(), plain), function(f) {
    set.seed(9L)
    bartisan(y ~ ., d, family = f, control = quick_control())
  })

  for (fit in fits) {
    expect_true(all(fit[["aux"]][, "power"] == 1.5))
  }

  expect_equal(fits[[1L]][["eta"]], fits[[2L]][["eta"]])
  expect_equal(fits[[1L]][["eta"]], fits[[3L]][["eta"]])

  # And an explicit request for a drawn power is still honored.
  set.seed(9L)
  drawn <- bartisan(y ~ ., d, family = tweedie(power = NULL),
                    control = quick_control(num_burn = 100L, num_draws = 100L))
  expect_gt(stats::var(drawn[["aux"]][, "power"]), 0)
})

test_that("the series holds up across the whole interval the bound allows", {
  # `power` may be anywhere strictly inside (1, 2), and the series' shape
  # parameter (2 - p) / (p - 1) runs from 999 to 0.001 over that range, so the
  # edges are where a term selection or an overflow would show up.
  y <- c(0, 0.5, 5, 200)
  eta <- matrix(c(0.5, 1, 2, 4), 1L)

  for (power in c(1.001, 1.01, 1.1, 1.9, 1.99, 1.999)) {
    got <- tweedie_engine_logdens(y, as.vector(eta), 4, power)

    expect_true(all(is.finite(got)),
                label = sprintf("power = %g", power))

    # And the zero is still the closed form, which is the part of the density
    # that has no series in it at all.
    expect_equal(got[1L], -exp(0.5)^(2 - power) / (4 * (2 - power)),
                 tolerance = 1e-10)
  }
})

test_that("the response has to be one the family can describe", {
  d <- sim_x(n = 80, seed = 91)

  d$y <- c(-1, stats::runif(nrow(d) - 1L))
  expect_error(bartisan(y ~ ., d, family = tweedie(), control = quick_control()),
               "non-negative")

  d$y <- rep(0, nrow(d))
  expect_error(bartisan(y ~ ., d, family = tweedie(), control = quick_control()),
               "positive")

  # No zeros is not an error, since the likelihood is still proper; it is a
  # warning because the point mass is then doing nothing.
  d$y <- stats::runif(nrow(d), 1, 5)
  expect_warning(bartisan(y ~ ., d, family = tweedie(),
                          control = quick_control()),
                 "no zeros")
})

test_that("a fit recovers the mean function and the dispersion", {
  skip_on_cran()

  phi <- 4
  power <- 1.5
  d <- sim_x(n = 1200, p = 4, seed = 411)
  truth <- 2 + sin(pi * d$x1) + 0.8 * d$x2 - 0.6 * d$x3
  set.seed(1411)
  d$y <- rtweedie_ref(nrow(d), exp(truth), phi, power)

  # The family is only interesting where there are zeros to model.
  expect_gt(mean(d$y == 0), 0.05)

  fit <- bartisan(y ~ ., d, family = tweedie(), chains = 2,
                  control = quick_control(num_trees = 50L, num_burn = 400L,
                                          num_draws = 600L))

  expect_identical(colnames(fit[["aux"]]), c("phi", "power"))

  # Fixed by default means fixed, not started there.
  expect_true(all(fit[["aux"]][, "power"] == power))

  eta <- colMeans(fit[["eta"]][[1L]])
  expect_gt(stats::cor(eta, truth), 0.9)
  expect_lt(abs(mean(eta - truth)), 0.15)
  expect_equal(mean(fit[["aux"]][, "phi"]), phi, tolerance = 0.25)

  # The mean is exp(eta) with nothing else in it, which is what makes a
  # counterfactual mean cheap. Compared draw by draw rather than against
  # `exp(colMeans(eta))`, which is a different quantity: the posterior mean of
  # exp(eta) is not exp of the posterior mean.
  expect_equal(predict(fit, draws = TRUE), exp(fit[["eta"]][[1L]]),
               tolerance = 1e-8)

  expect_predictor_invariant(fit, d)
})

test_that("the posterior predictive puts back the point mass", {
  skip_on_cran()

  d <- sim_x(n = 800, p = 3, seed = 511)
  set.seed(1511)
  d$y <- rtweedie_ref(nrow(d), exp(2 + sin(pi * d$x1)), 4, 1.5)

  fit <- bartisan(y ~ ., d, family = tweedie(), chains = 2,
                  control = quick_control(num_trees = 50L, num_burn = 300L,
                                          num_draws = 400L))

  reps <- as.matrix(stats::simulate(fit, nsim = 100L))

  # A replicate is a draw from the same family, so it has to carry a point mass
  # of its own: this is the check that `simulate()` uses the representation and
  # not a continuous approximation to it.
  expect_true(all(reps >= 0))
  expect_gt(mean(reps == 0), 0)
  expect_equal(mean(reps == 0), mean(d$y == 0), tolerance = 0.05)
  expect_equal(mean(colMeans(reps)), mean(d$y), tolerance = 0.1 * mean(d$y))
})

test_that("the family is reachable by name and drawn power moves", {
  skip_on_cran()

  d <- sim_x(n = 600, p = 3, seed = 611)
  set.seed(1611)
  d$y <- rtweedie_ref(nrow(d), exp(2 + sin(pi * d$x1)), 4, 1.6)

  by_name <- bartisan(y ~ ., d, family = "tweedie",
                      control = quick_control(num_burn = 50L, num_draws = 50L))
  expect_identical(by_name[["family"]][["family"]], "tweedie")
  expect_true(all(by_name[["aux"]][, "power"] == 1.5))

  drawn <- bartisan(y ~ ., d, family = tweedie(power = NULL),
                    control = quick_control(num_trees = 50L, num_burn = 300L,
                                            num_draws = 300L))
  p <- drawn[["aux"]][, "power"]
  expect_gt(stats::var(p), 0)
  expect_true(all(p > 1 & p < 2))
})

# The cell of the audit matrix that a new family most needs filled in. A `vc()`
# model failed to refresh its inner family's augmentation for a year without
# anything looking wrong, and the tell was a *positive* log likelihood, so every
# structure gets checked for one here.
test_that("the family survives every structure that wraps a forest", {
  skip_on_cran()

  set.seed(31L)
  n <- 700L
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  z = stats::rbinom(n, 1L, 0.5),
                  g = factor(sample(letters[1:6], n, TRUE)),
                  off = stats::rnorm(n, 0, 0.1),
                  w = stats::runif(n, 0.5, 1.5))
  truth <- 1.8 + sin(pi * d$x1) + 0.7 * d$z
  set.seed(131L)
  d$y <- rtweedie_ref(n, exp(truth), 4, 1.5)

  ctrl <- function(...) {
    quick_control(num_trees = 20L, num_burn = 150L, num_draws = 200L, ...)
  }

  fits <- list(
    plain = bartisan(y ~ x1 + x2 + z, d, family = tweedie(), control = ctrl()),
    weights = bartisan(y ~ x1 + x2 + z, d, family = tweedie(),
                       weights = d$w, control = ctrl()),
    offset = bartisan(y ~ x1 + x2 + z + offset(off), d, family = tweedie(),
                      control = ctrl()),
    ranef = bartisan(y ~ x1 + x2 + z + (1 | g), d, family = tweedie(),
                     control = ctrl()),
    vc = bartisan(y ~ vc(z, ~ x1) + x1 + x2, d, family = tweedie(),
                  control = ctrl()),
    hard = bartisan(y ~ x1 + x2 + z, d, family = tweedie(),
                    control = ctrl(gate = "hard")),
    drawn = bartisan(y ~ x1 + x2 + z, d, family = tweedie(power = NULL),
                     control = ctrl()))

  for (nm in names(fits)) {
    fit <- fits[[nm]]

    # A density with a point mass at zero can exceed 1 there, so a positive log
    # likelihood is not impossible in principle; at this dispersion the point
    # mass is nowhere near 1, and a positive total is the signature of a
    # likelihood evaluated with the wrong nuisance parameters.
    expect_lt(mean(fit[["loglik"]]), 0, label = nm)
    expect_true(all(is.finite(fit[["loglik"]])), label = nm)
    expect_true(all(fit[["aux"]][, "phi"] > 0), label = nm)

    # Two are left out for the reasons `expect_predictor_invariant()` gives:
    # an offset fit needs its offset passed to `predict()`, and a `vc()` fit
    # reports one combined predictor where `eta` keeps the forests apart.
    if (!nm %in% c("offset", "vc")) {
      expect_predictor_invariant(fit, d)
    }
  }

  # The structures that leave the mean alone should still recover it.
  for (nm in c("plain", "weights", "ranef", "hard", "drawn")) {
    expect_gt(stats::cor(colMeans(fits[[nm]][["eta"]][[1L]]), truth), 0.9,
              label = nm)
  }

  # And the varying coefficient is the quantity the wrapper exists to report.
  expect_equal(mean(coef(fits[["vc"]])[, 1L]), 0.7, tolerance = 0.25)
})
