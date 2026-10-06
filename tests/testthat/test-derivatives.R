# The sampler proposes leaf values from a Gaussian fitted by Fisher scoring, so
# a wrong analytic score sends the proposal to the wrong place. These tests
# check every family's analytic score against a central difference of its own
# log density, which is an independent route to the same quantity.
#
# The information is checked only for the families that report the observed
# second derivative. Several deliberately report the *expected* information
# instead, because the observed version can be negative -- probit and cloglog
# binomial, negative binomial, the zero-inflated mixtures, and and both beta
# families. For those, disagreement is the intended behavior, so
# only the score is compared. The gamma family is not among them: its response
# is strictly positive, so its observed curvature cannot go negative, and it
# reports the true one.

derivs <- function(family, link, y, eta, opts = list(),
                   aux = matrix(0, 1L, 0L), component = 0L,
                   by_difference = FALSE) {
  weights <- rep(1, length(y))
  .bartisan_derivs(y, weights, eta, family, link, opts, aux,
                   as.integer(component), by_difference, FALSE)
}

expect_score_matches_difference <- function(family, link, y, eta, opts = list(),
                                            aux = matrix(0, 1L, 0L),
                                           component = 0L,
                                           check_info = FALSE,
                                           tolerance = 1e-6) {
  analytic <- derivs(family, link, y, eta, opts, aux, component, FALSE)
  numeric <- derivs(family, link, y, eta, opts, aux, component, TRUE)

  scale <- max(1, max(abs(numeric[["d1"]])))
  testthat::expect_lt(max(abs(analytic[["d1"]] - numeric[["d1"]])) / scale,
                      tolerance)

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

# The two beta families take their score from a per-sweep table of the two
# digamma combinations, interpolated linearly, so the score matches a central
# difference of the exact log density only to the table's accuracy. On these
# grids that is 5.5e-6 of the score's scale for beta and 3.5e-6 for ordered
# beta, and a table read one entry off would be 1.8e-2, so 2e-5 separates the
# two. The responses are a fixed grid: drawn at random, the error crossed 1e-6
# for more than half of 500 draws.
y_grid <- seq(0.05, 0.95, length.out = n_grid)

test_that("the beta score is correct, tabulated derivatives and all", {
  expect_score_matches_difference(
    "beta", "logit", y_grid, list(grid),
    list(phi = 8, phi_prior_shape = 0.01, phi_prior_rate = 0.01,
         update_phi = TRUE),
    matrix(8, 1L, 1L), tolerance = 2e-5)
})

test_that("the ordered beta score is correct at both endpoints and inside", {
  y <- y_grid
  y[1:6] <- 0
  y[7:12] <- 1

  expect_score_matches_difference(
    "ordbeta", "logit", y, list(grid),
    list(cut1 = -1.5, cut2 = 1.5, phi = 8, phi_prior_shape = 0.01,
         phi_prior_rate = 0.01, update_phi = TRUE),
    matrix(c(-1.5, 1.5, 8), 1L, 3L), tolerance = 2e-5)
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

  expect_score_matches_difference("gaussian_ls", "identity",
                                  stats::rnorm(n_grid), list(grid, second),
                                  component = 0L, check_info = TRUE)
  # The gamma location-scale reports the observed curvature in the mean, as
  # `Gamma("log")` does, and the *expected* information in the log dispersion,
  # where the observed one can go negative away from the mode. So the score is
  # checked in both components and the information only in the first.
  expect_score_matches_difference("Gamma_ls", "log",
                                  stats::rgamma(n_grid, shape = 2, rate = 1),
                                  list(grid, second),
                                  component = 0L, check_info = TRUE)
  expect_score_matches_difference("Gamma_ls", "log",
                                  stats::rgamma(n_grid, shape = 2, rate = 1),
                                  list(grid, second),
                                  component = 1L)
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

# Where a stable form would lose its precision, the families hand over: the
# ordinal family to a difference of its log density once a category's
# probability underflows, and the ordered beta to the exact digamma combinations
# beyond the reach of its table, which ends at a predictor of 8.
test_that("the scores stay right where the stable forms hand over", {
  # Under the cloglog link the top category's probability is exp(-exp(c - eta)),
  # zero in double precision by eta = -9 though its log is finite, and the score
  # there is exp(c - eta).
  eta <- matrix(c(-9, -10), nrow = 1L)
  opts <- list(num_cat = 4L, cuts = c(0, 1, 2.5), update_cuts = TRUE)
  aux <- matrix(c(0, 1, 2.5), 1L, 3L)

  expect_score_matches_difference("ordinal", "cloglog", c(3, 3), list(eta),
                                  opts, aux)
  expect_equal(as.vector(derivs("ordinal", "cloglog", c(3, 3), list(eta), opts,
                                aux)[["d1"]]),
               exp(2.5 - as.vector(eta)), tolerance = 1e-6)

  expect_score_matches_difference(
    "ordbeta", "logit", rep(c(0.3, 0.7), 2L),
    list(matrix(c(-12, -9, 9, 12), nrow = 1L)),
    list(cut1 = -1.5, cut2 = 1.5, phi = 8, phi_prior_shape = 0.01,
         phi_prior_rate = 0.01, update_phi = TRUE),
    matrix(c(-1.5, 1.5, 8), 1L, 3L), tolerance = 2e-5)
})

# Past a predictor of about 37, expit() rounds to one, and three things built on
# it went with it: the logistic density in the ordinal score was p * (1 - p) and
# came out zero, so the score of the top category at -40 read 0 where it is 1;
# and both beta families formed their second shape as phi minus the first, which
# was zero there, so the log density was -Inf where it is finite and the score
# was NaN. Each is checked against something computed independently in R.
test_that("the logistic and beta families stay right far into the tails", {
  # The ordinal logit scores of the extreme categories have closed forms:
  # expit(c - eta) for the top one and -expit(eta - c) for the bottom one.
  cuts <- c(0, 1, 2.5)
  eta <- matrix(c(-40, -100, 40, 100), nrow = 1L)
  y <- c(3, 3, 0, 0)
  got <- derivs("ordinal", "logit", y, list(eta),
                list(num_cat = 4L, cuts = cuts, update_cuts = TRUE),
                matrix(cuts, 1L, 3L))
  expect_equal(as.vector(got[["d1"]]),
               c(stats::plogis(2.5 - eta[1:2]), -stats::plogis(eta[3:4] - 0)),
               tolerance = 1e-10)

  # The beta log density is R's, and finite; the score matches its difference.
  phi <- 8
  beta_eta <- matrix(c(-100, -40, 40, 100), nrow = 1L)
  beta_y <- c(0.3, 0.7, 0.3, 0.7)
  beta_opts <- list(phi = phi, phi_prior_shape = 0.01, phi_prior_rate = 0.01,
                    update_phi = TRUE)
  dens <- .bartisan_logdens(beta_y, rep(1, 4L), list(beta_eta), "beta",
                            "logit", beta_opts, matrix(phi, 1L, 1L))

  expect_true(all(is.finite(dens)))
  expect_equal(as.vector(dens),
               stats::dbeta(beta_y, phi * stats::plogis(beta_eta),
                            phi * stats::plogis(-beta_eta), log = TRUE),
               tolerance = 1e-10)
  expect_score_matches_difference("beta", "logit", beta_y, list(beta_eta),
                                  beta_opts, matrix(phi, 1L, 1L))

  # The ordered beta's interior density is the middle interval's probability
  # times the beta density. The interval is taken on whichever side keeps it a
  # difference of two small numbers rather than of two near one.
  ord_opts <- list(cut1 = -1.5, cut2 = 1.5, phi = phi, phi_prior_shape = 0.01,
                   phi_prior_rate = 0.01, update_phi = TRUE)
  ord_aux <- matrix(c(-1.5, 1.5, phi), 1L, 3L)
  dens <- .bartisan_logdens(beta_y, rep(1, 4L), list(beta_eta), "ordbeta",
                            "logit", ord_opts, ord_aux)
  e <- as.vector(beta_eta)
  span <- ifelse(e > 0, stats::plogis(1.5 - e) - stats::plogis(-1.5 - e),
                 stats::plogis(e + 1.5) - stats::plogis(e - 1.5))

  expect_true(all(is.finite(dens)))
  expect_equal(as.vector(dens),
               log(span) + stats::dbeta(beta_y, phi * stats::plogis(e),
                                        phi * stats::plogis(-e), log = TRUE),
               tolerance = 1e-10)
  expect_score_matches_difference("ordbeta", "logit", beta_y, list(beta_eta),
                                  ord_opts, ord_aux)
})
