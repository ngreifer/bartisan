# The hurdle families. A zero comes from the zero part and a positive value from
# the positive part, so each check below is one the factorization makes exact:
# the density is the zero probability at a zero and the complement times the
# positive part's density elsewhere, truncated at zero for a count; the mean is
# the complement times the positive part's mean.

# Zeros from a logistic model in x2, and positive values from a truncated
# Poisson and from a gamma with a log mean in x1. The truncated count is drawn by
# inversion on the upper tail, as `posterior_predict()` draws it.
sim_hurdle <- function(n = 120, seed = 1) {
  d <- sim_x(n = n, seed = seed)
  zero <- stats::runif(n) < stats::plogis(-0.5 + 1.5 * (d$x2 - 0.5))
  mu <- exp(0.3 + 1.2 * d$x1)
  tail <- stats::ppois(0, mu, lower.tail = FALSE)
  d$count <- ifelse(zero, 0, stats::qpois(stats::runif(n, 0, tail), mu,
                                          lower.tail = FALSE))
  d$cost <- ifelse(zero, 0, stats::rgamma(n, shape = 2, rate = 2 / mu))
  attr(d, "zero_prob") <- stats::plogis(-0.5 + 1.5 * (d$x2 - 0.5))
  d
}

# The density written with R's own d* functions, draws by observations.
hurdle_reference <- function(family, y, eta, aux) {
  e <- eta[[1L]]
  mu <- exp(e)
  y <- matrix(y, nrow = nrow(e), ncol = ncol(e), byrow = TRUE)

  positive <- switch(family,
                     hurdle_poisson = stats::dpois(y, mu) /
                       stats::ppois(0, mu, lower.tail = FALSE),
                     hurdle_negbin = {
                       theta <- aux[, "theta"]
                       stats::dnbinom(y, size = theta, mu = mu) /
                         stats::pnbinom(0, size = theta, mu = mu,
                                        lower.tail = FALSE)
                     },
                     hurdle_Gamma = {
                       shape <- aux[, "shape"]
                       stats::dgamma(y, shape = shape, rate = shape / mu)
                     },
                     hurdle_lognormal = stats::dlnorm(y, e, aux[, "sigma"]))

  p_zero <- stats::plogis(eta[[2L]])
  ifelse(y == 0, p_zero, (1 - p_zero) * positive)
}

test_that("the hurdle constructors resolve and check their arguments", {
  for (name in c("hurdle_poisson", "hurdle_negbin", "hurdle_Gamma",
                 "hurdle_lognormal")) {
    expect_identical(as_bartisan_family(name)[["family"]], name)
  }

  expect_identical(as_bartisan_family(hurdle_negbin(theta = 2))[["theta"]], 2)
  expect_identical(hurdle_lognormal()[["link"]], "identity")

  expect_error(hurdle_poisson("identity"), "should be")
  expect_error(hurdle_Gamma("inverse"), "should be")
  expect_error(hurdle_negbin(theta = -1), "must be")
})

test_that("each hurdle family's density, mean and forests are the ones it claims", {
  skip_if_not_installed("patrick")

  d <- sim_hurdle()

  patrick::with_parameters_test_that(
    "{name}:",
    {
      fit <- bartisan(stats::reformulate(c("x1", "x2", "x3"), response),
                      data = d, family = family, control = quick_control())

      expect_identical(fit[["num_forest"]], 2L)
      expect_identical(names(fit[["eta"]]), forests)
      expect_identical(colnames(fit[["aux"]]), aux_names)
      expect_predictor_invariant(fit, d)

      eta <- predict(fit, type = "link", draws = TRUE)
      dens <- predict(fit, type = "density", draws = TRUE)
      expect_equal(dens, hurdle_reference(name, d[[response]], eta, fit[["aux"]]),
                   ignore_attr = TRUE)

      # The reported log likelihood is the density's, normalizing constants
      # and all.
      expect_equal(rowSums(log(dens)), fit[["loglik"]], tolerance = 1e-8)

      p_zero <- stats::plogis(eta[[2L]])
      mean_positive <- switch(name,
                              hurdle_poisson = exp(eta[[1L]]) /
                                -expm1(-exp(eta[[1L]])),
                              hurdle_negbin = {
                                theta <- fit[["aux"]][, "theta"]
                                mu <- exp(eta[[1L]])
                                mu / stats::pnbinom(0, size = theta, mu = mu,
                                                    lower.tail = FALSE)
                              },
                              hurdle_Gamma = exp(eta[[1L]]),
                              hurdle_lognormal = exp(eta[[1L]] +
                                                       fit[["aux"]][, "sigma"]^2 / 2))
      expect_equal(predict(fit, type = "response", draws = TRUE),
                   (1 - p_zero) * mean_positive, ignore_attr = TRUE)
    },
    patrick::cases(
      list(name = "hurdle_poisson", family = hurdle_poisson(),
           response = "count", forests = c("count", "zero"), aux_names = NULL),
      list(name = "hurdle_negbin", family = hurdle_negbin(),
           response = "count", forests = c("count", "zero"),
           aux_names = "theta"),
      list(name = "hurdle_Gamma", family = hurdle_Gamma(), response = "cost",
           forests = c("mean", "zero"), aux_names = "shape"),
      list(name = "hurdle_lognormal", family = hurdle_lognormal(),
           response = "cost", forests = c("mean", "zero"),
           aux_names = "sigma")
    )
  )
})

test_that("a hurdle response is checked for what the family needs", {
  d <- sim_hurdle(n = 40, seed = 2)
  ctrl <- quick_control()

  d$y <- d$count
  d$y[1L] <- -1
  expect_error(bartisan(y ~ x1, data = d, family = hurdle_poisson(),
                        control = ctrl), "non-negative")

  d$y <- d$count + 0.5
  expect_error(bartisan(y ~ x1, data = d, family = hurdle_negbin(),
                        control = ctrl), "counts")

  d$y <- 0
  expect_error(bartisan(y ~ x1, data = d, family = hurdle_Gamma(),
                        control = ctrl), "positive responses")

  d$y <- d$cost + 1
  expect_warning(bartisan(y ~ x1, data = d, family = hurdle_lognormal(),
                          control = ctrl), "no zeros")
})

test_that("the untruncated Poisson start inverts the truncated mean", {
  for (mu in c(0.01, 0.5, 3, 40)) {
    truncated <- mu / -expm1(-mu)
    expect_equal(untruncated_poisson_mean(truncated), mu, tolerance = 1e-6)
  }

  expect_gt(untruncated_poisson_mean(1), 0)
})

test_that("hurdle replicates are zero as often as the zero part says", {
  skip_on_cran()
  skip_if_not_installed("patrick")
  skip_if_not_installed("rstantools")

  d <- sim_hurdle(n = 200, seed = 3)

  patrick::with_parameters_test_that(
    "{name}:",
    {
      fit <- bartisan(stats::reformulate(c("x1", "x2", "x3"), response),
                      data = d, family = family,
                      control = quick_control(num_draws = 100L))

      set.seed(30)
      replicates <- rstantools::posterior_predict(fit)
      p_zero <- stats::plogis(predict(fit, type = "link", draws = TRUE)[[2L]])

      # 20,000 cells, so the share of zeros has a standard error near .003.
      expect_equal(mean(replicates == 0), mean(p_zero), tolerance = 0.02)
      expect_true(all(replicates >= 0))

      if (counts) {
        expect_true(all(replicates == round(replicates)))
      }
    },
    patrick::cases(
      list(name = "hurdle_poisson", family = hurdle_poisson(),
           response = "count", counts = TRUE),
      list(name = "hurdle_negbin", family = hurdle_negbin(),
           response = "count", counts = TRUE),
      list(name = "hurdle_Gamma", family = hurdle_Gamma(), response = "cost",
           counts = FALSE),
      list(name = "hurdle_lognormal", family = hurdle_lognormal(),
           response = "cost", counts = FALSE)
    )
  )
})

# The rewriting has to leave the posterior alone. Point estimates alone would
# pass a rewriting that understated the uncertainty, so the spread is compared
# too, as it is for the zero-inflated families in test-augment.R.
test_that("the hurdle rewriting agrees with the direct family", {
  skip_on_cran()
  skip_if_not_installed("patrick")

  d <- sim_hurdle(n = 300, seed = 4)
  zero_prob <- attr(d, "zero_prob")
  mu <- exp(0.3 + 1.2 * d$x1)

  patrick::with_parameters_test_that(
    "{name}:",
    {
      ctrl <- function(a) {
        bartisan_control(num_trees = 50, num_burn = 400, num_draws = 400,
                         augment = a, verbose = FALSE)
      }

      set.seed(5)
      direct <- bartisan(stats::reformulate(c("x1", "x2", "x3"), response),
                         data = d, family = family, control = ctrl(FALSE))
      set.seed(5)
      rewritten <- bartisan(stats::reformulate(c("x1", "x2", "x3"), response),
                            data = d, family = family, control = ctrl(TRUE))

      for (h in 1:2) {
        a <- colMeans(direct[["eta"]][[h]])
        b <- colMeans(rewritten[["eta"]][[h]])
        expect_gt(stats::cor(a, b), 0.95, label = paste("forest", h))
        expect_equal(mean(apply(rewritten[["eta"]][[h]], 2L, stats::sd)),
                     mean(apply(direct[["eta"]][[h]], 2L, stats::sd)),
                     tolerance = 0.3, label = paste("forest", h, "spread"))
      }

      # Both recover the zero probability, which is the same logistic model
      # whatever the positive part is.
      for (fit in list(direct, rewritten)) {
        p_zero <- colMeans(stats::plogis(fit[["eta"]][[2L]]))
        expect_gt(stats::cor(p_zero, zero_prob), 0.8)
        expect_gt(stats::cor(colMeans(fit[["eta"]][[1L]]), log(mu)), 0.9)
      }
    },
    patrick::cases(
      list(name = "hurdle_poisson", family = hurdle_poisson(),
           response = "count"),
      list(name = "hurdle_negbin", family = hurdle_negbin(),
           response = "count"),
      list(name = "hurdle_Gamma", family = hurdle_Gamma(), response = "cost"),
      list(name = "hurdle_lognormal", family = hurdle_lognormal(),
           response = "cost")
    )
  )
})
