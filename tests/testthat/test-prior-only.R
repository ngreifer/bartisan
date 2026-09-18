# `prior_only = TRUE` hands the engine a zero weight for every observation. The
# weight multiplies that observation's log density, its gradient and its
# curvature in one place, so the likelihood goes flat and the sampler targets
# the prior instead of the posterior.

test_that("a prior-only fit draws from the prior rather than the posterior", {
  skip_on_cran()
  skip_if_not_installed("rstantools")

  set.seed(2026)
  d <- sim_x(n = 250L, p = 3L, seed = 110L)
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(2 * d$x1 - 1))

  ctrl <- quick_control(num_trees = 10L, num_burn = 60L, num_draws = 200L)

  prior <- bartisan(y ~ ., d, family = stats::binomial(), control = ctrl,
                    prior_only = TRUE)
  post <- bartisan(y ~ ., d, family = stats::binomial(), control = ctrl)

  expect_s3_class(prior, "bartisan_fit")
  expect_true(isTRUE(prior[["prior_only"]]))
  expect_false(isTRUE(post[["prior_only"]]))

  # Wider uncertainty per observation, because the prior has learned nothing.
  # This is the signature to test, not the spread of the averaged prediction:
  # `predict()` averages over draws, and under the prior the forest wanders
  # symmetrically about its anchor, so the averaged function comes out flatter
  # than the posterior's rather than wider.
  uncertainty <- function(f) {
    mean(apply(stats::predict(f, type = "response", draws = TRUE), 2L,
               stats::sd))
  }

  expect_gt(uncertainty(prior), 1.5 * uncertainty(post))

  # The mean function is nearly constant across observations, being an average
  # over draws of a function the prior gives no reason to bend one way.
  flatness <- function(f) {
    stats::sd(stats::predict(f, type = "response"))
  }

  expect_lt(flatness(prior), flatness(post))

  # And it has not been fitted: the posterior tracks x1 and the prior does not.
  cor_with_x1 <- function(f) {
    abs(stats::cor(stats::predict(f, type = "response"), d$x1))
  }

  expect_gt(cor_with_x1(post), cor_with_x1(prior))
})

# The weight is a convenience and not the only way to flatten a likelihood. An
# update written against the response directly is flattened by telling it the
# observation carries no weight, which is what `dpm()` and the multinomial
# probit's latent utilities now do. The test of whether that worked is that
# permuting the response changes nothing: the intercept anchor takes the mean
# and the sd of `y`, and a permutation preserves both, so any movement in the
# draws is the data reaching them by some other route.
prior_only_eta <- function(y, family, data, seed = 7L) {
  set.seed(seed)

  fit <- bartisan(y ~ ., cbind(data, y = y), family = family,
                  prior_only = TRUE,
                  control = quick_control(num_trees = 10L, num_burn = 60L,
                                          num_draws = 120L))

  as.numeric(unlist(stats::predict(fit, type = "link", draws = TRUE)))
}

test_that("a permuted response leaves a prior-only fit where it was", {
  skip_on_cran()

  d <- sim_x(n = 150L, p = 2L, seed = 113L)
  y <- stats::rnorm(nrow(d), d$x1)

  set.seed(3)
  permuted <- sample(y)

  # Weighted throughout, so the two runs agree to within the last bits of the
  # anchor, which is a mean over the response and sums in a different order once
  # the response has been permuted.
  expect_equal(prior_only_eta(y, stats::gaussian(), d),
               prior_only_eta(permuted, stats::gaussian(), d),
               tolerance = 1e-8)

  expect_equal(prior_only_eta(y, gaussian_ls(), d),
               prior_only_eta(permuted, gaussian_ls(), d),
               tolerance = 1e-8)

  # The mixture consumes random numbers as it goes, so the two streams drift
  # apart by a hair and exact equality is the wrong bar. The bar that matters is
  # that the drift is nothing beside the spread of the prior itself: before the
  # atoms were told to ignore weightless residuals this gap was larger than the
  # prior's own standard deviation.
  a <- prior_only_eta(y, dpm(), d)
  b <- prior_only_eta(permuted, dpm(), d)

  expect_true(all(is.finite(c(a, b))))
  expect_lt(max(abs(a - b)), 0.25 * stats::sd(a))
})

test_that("every family can be drawn from its prior", {
  skip_on_cran()

  d <- sim_x(n = 120L, p = 2L, seed = 111L)
  d$y <- stats::rnorm(nrow(d))
  d$pos <- stats::rexp(nrow(d))
  d$ord <- factor(sample(1:3, nrow(d), replace = TRUE), ordered = TRUE)
  d$cat <- factor(sample(letters[1:3], nrow(d), replace = TRUE))
  d$unit <- stats::runif(nrow(d))
  d$unit[seq_len(12L)] <- 0
  d$unit[13:24] <- 1

  ctrl <- quick_control(num_trees = 5L, num_burn = 30L, num_draws = 50L)

  # The two cutpoint families were refused until the cutpoints carried the
  # induced-Dirichlet prior. With it, a flat likelihood leaves a proper density
  # to draw from and the replicates stay on the scale the response is measured
  # on rather than piling at one end of it.
  for (spec in list(list(ord ~ x1, ordinal()), list(unit ~ x1, ordbeta()))) {
    fit <- suppressWarnings(bartisan(spec[[1L]], d, family = spec[[2L]],
                                     control = ctrl, prior_only = TRUE))
    expect_s3_class(fit, "bartisan_fit")

    rep <- rstantools::posterior_predict(fit)
    expect_true(all(is.finite(rep)))
    # Every category is reachable under the prior rather than one absorbing it.
    expect_gt(length(unique(round(as.vector(rep), 6))), 1L)
  }

  # Everything else is reachable, whether the weight carries the whole
  # likelihood or an auxiliary update had to be told about it separately.
  for (family in list(stats::gaussian(), dpm(), gaussian_ls())) {
    expect_s3_class(bartisan(y ~ x1, d, family = family, control = ctrl,
                             prior_only = TRUE),
                    "bartisan_fit")
  }

  expect_s3_class(bartisan(pos ~ x1, d, family = Gamma_ls(), control = ctrl,
                           prior_only = TRUE),
                  "bartisan_fit")

  # The multinomial probit drew each latent utility with a variance of one over
  # the weight, which at zero weight left the covariance draw working on
  # infinities and warning that its matrix was no longer symmetric.
  #
  # Only that warning is the subject. A flat likelihood leaves the leaf scale
  # unidentified, so it can settle far above its prior median and warn about it
  # legitimately, and a bare `expect_no_warning()` fails on that instead
  # whenever the draws happen to run high.
  expect_no_warning(bartisan(cat ~ x1, d, family = multinomial("probit"),
                             control = ctrl, prior_only = TRUE),
                    message = "symmetric")
})

test_that("nothing that scores a fit against data runs on a prior-only one", {
  skip_on_cran()
  skip_if_not_installed("loo")
  skip_if_not_installed("rstantools")

  d <- sim_x(n = 120L, p = 2L, seed = 112L)
  d$y <- stats::rnorm(nrow(d), d$x1)

  prior <- bartisan(y ~ ., d, family = stats::gaussian(), prior_only = TRUE,
                    control = quick_control(num_trees = 5L, num_burn = 30L,
                                            num_draws = 60L))

  expect_error(loo::loo(prior), "never shown any")
  expect_error(loo::waic(prior), "never shown any")
  expect_error(loo::kfold(prior), "never shown any")

  # `model_performance()` reaches `loo()` on its own, so it drops those columns
  # rather than surfacing an error from a call that never mentioned loo.
  skip_if_not_installed("performance")

  perf <- suppressMessages(performance::model_performance(prior))

  expect_s3_class(perf, "data.frame")
  expect_false(any(c("ELPD", "LOOIC", "WAIC") %in% names(perf)))
})
