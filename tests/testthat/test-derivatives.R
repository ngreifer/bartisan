# The sampler proposes leaf values from a Gaussian fitted by Fisher scoring, so
# a wrong analytic score sends the proposal to the wrong place. These tests
# check every family's analytic score against a central difference of its own
# log density, which is an independent route to the same quantity.
#
# The information is checked only for the families that report the observed
# second derivative. Several deliberately report the *expected* information
# instead, because the observed version can be negative -- probit and cloglog
# binomial, negative binomial, the zero-inflated mixtures, and the beta part of
# the ordered beta family. For those, disagreement is the intended behavior, so
# only the score is compared. The gamma family is not among them: its response
# is strictly positive, so its observed curvature cannot go negative, and it
# reports the true one.

derivs <- function(family, link, y, eta, opts = list(),
                   aux = matrix(0, 1L, 0L), component = 0L,
                   by_difference = FALSE) {
  weights <- rep(1, length(y))
  .genbart_derivs(y, weights, eta, family, link, opts, aux,
                  as.integer(component), by_difference, FALSE)
}

expect_score_matches_difference <- function(family, link, y, eta, opts = list(),
                                           aux = matrix(0, 1L, 0L),
                                           component = 0L,
                                           check_info = FALSE) {
  analytic <- derivs(family, link, y, eta, opts, aux, component, FALSE)
  numeric <- derivs(family, link, y, eta, opts, aux, component, TRUE)

  scale <- max(1, max(abs(numeric[["d1"]])))
  testthat::expect_lt(max(abs(analytic[["d1"]] - numeric[["d1"]])) / scale, 1e-6)

  if (check_info) {
    scale <- max(1, max(abs(numeric[["info"]])))
    testthat::expect_lt(max(abs(analytic[["info"]] - numeric[["info"]])) / scale,
                        1e-4)
  }
}

# A spread of predictor values, including the tails where the stable forms of
# the log densities take their alternative branches.
grid <- matrix(seq(-3, 3, length.out = 41L), nrow = 1L)
n_grid <- ncol(grid)

test_that("single-predictor families have correct analytic scores", {
  set.seed(11)

  expect_score_matches_difference("gaussian", "identity", stats::rnorm(n_grid),
                                  list(grid), list(sigma_hat = 1),
                                  matrix(1, 1L, 1L), check_info = TRUE)

  binary <- stats::rbinom(n_grid, 1, 0.5)
  expect_score_matches_difference("binomial", "logit", binary, list(grid),
                                  check_info = TRUE)
  expect_score_matches_difference("binomial", "probit", binary, list(grid))
  expect_score_matches_difference("binomial", "cloglog", binary, list(grid))

  counts <- stats::rpois(n_grid, 2)
  expect_score_matches_difference("poisson", "log", counts, list(grid),
                                  check_info = TRUE)
  expect_score_matches_difference("negbin", "log", counts, list(grid),
                                  list(theta = 2, theta_prior_shape = 0.01,
                                       theta_prior_rate = 0.01,
                                       update_theta = TRUE),
                                  matrix(2, 1L, 1L))
  # The gamma family reports the observed second derivative, not the expected
  # information, so its information is checked too.
  expect_score_matches_difference("Gamma", "log", stats::rgamma(n_grid, 2, 1),
                                  list(grid),
                                  list(shape = 2, shape_prior_shape = 0.01,
                                       shape_prior_rate = 0.01,
                                       update_shape = TRUE),
                                  matrix(2, 1L, 1L), check_info = TRUE)
})

test_that("the ordinal scores are correct for both links", {
  set.seed(12)

  y <- sample(0:3, n_grid, replace = TRUE)
  opts <- list(num_cat = 4L, cuts = c(0, 1, 2.5), update_cuts = TRUE)
  aux <- matrix(c(0, 1, 2.5), 1L, 3L)

  for (link in c("logit", "probit")) {
    expect_score_matches_difference("ordinal", link, y, list(grid), opts, aux,
                                    check_info = TRUE)
  }
})

test_that("the ordered beta score is correct at both endpoints and inside", {
  set.seed(13)

  y <- stats::runif(n_grid, 0.05, 0.95)
  y[1:6] <- 0
  y[7:12] <- 1

  expect_score_matches_difference(
    "ordbeta", "logit", y, list(grid),
    list(cut1 = -1.5, cut2 = 1.5, phi = 8, phi_prior_shape = 0.01,
         phi_prior_rate = 0.01, update_phi = TRUE),
    matrix(c(-1.5, 1.5, 8), 1L, 3L))
})

test_that("the survival scores are correct for events and for censoring", {
  set.seed(14)

  y <- stats::rnorm(n_grid)
  # Alternate events and censored observations so both branches are exercised.
  opts <- list(event = rep(c(1, 0), length.out = n_grid), sigma_hat = 1,
               update_sigma = TRUE)

  for (link in c("weibull", "loglogistic", "lognormal")) {
    expect_score_matches_difference("aft", link, y, list(grid), opts,
                                    matrix(1, 1L, 1L), check_info = TRUE)
  }
})

test_that("multi-predictor families are correct in each component", {
  set.seed(15)

  second <- matrix(stats::rnorm(n_grid) * 0.4, nrow = 1L)

  expect_score_matches_difference("location_scale", "identity",
                                  stats::rnorm(n_grid), list(grid, second),
                                  component = 0L, check_info = TRUE)
  # Both parameterizations. Reference coding carries num_cat - 1 predictors and
  # the symmetric coding one per category, so the same two-column eta describes
  # a three-category response under the first and a two-category one under the
  # second.
  expect_score_matches_difference("multinomial", "logit",
                                  sample(0:2, n_grid, replace = TRUE),
                                  list(grid, second),
                                  list(num_cat = 3L, symmetric = FALSE),
                                  component = 1L, check_info = TRUE)
  expect_score_matches_difference("multinomial", "logit",
                                  sample(0:1, n_grid, replace = TRUE),
                                  list(grid, second),
                                  list(num_cat = 2L, symmetric = TRUE),
                                  component = 1L, check_info = TRUE)

  counts <- stats::rpois(n_grid, 2)
  for (component in 0:1) {
    expect_score_matches_difference("zip", "log", counts,
                                    list(grid, second), component = component)
    expect_score_matches_difference(
      "zinb", "log", counts, list(grid, second),
      list(theta = 2, theta_prior_shape = 0.01, theta_prior_rate = 0.01,
           update_theta = TRUE), matrix(2, 1L, 1L), component = component)
  }
})
