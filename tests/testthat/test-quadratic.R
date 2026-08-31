# The closed forms a quadratic target allows. A leaf value enters the additive
# predictor linearly, so when the family's log density is quadratic in the
# predictor the log target is quadratic in the leaf value -- and the leaf prior
# is Gaussian, so it is too. One pass over a node then determines the whole
# function exactly, which is where a conjugate sampler's advantage comes from.
#
# The point of these tests is that "exactly" is meant literally: the closed
# forms must reproduce what the general Fisher-scoring path computes, not merely
# come close.

test_that("the closed forms reproduce the general path", {
  d <- sim_x(n = 200, seed = 101)
  d$y <- 2 * d$x1 + stats::rnorm(200, sd = 0.4)

  settings <- list(
    list(label = "gaussian, soft", family = gaussian(), gate = "smoothstep",
         augment = FALSE, y = d$y),
    list(label = "gaussian, hard", family = gaussian(), gate = "hard",
         augment = FALSE, y = d$y),
    list(label = "probit augmented, soft", family = binomial("probit"),
         gate = "smoothstep", augment = TRUE,
         y = stats::rbinom(200, 1, stats::plogis(2 * d$x1 - 1))),
    list(label = "logit augmented, hard", family = binomial("logit"),
         gate = "hard", augment = TRUE,
         y = stats::rbinom(200, 1, stats::plogis(2 * d$x1 - 1))),
    list(label = "location-scale, soft", family = location_scale(),
         gate = "smoothstep", augment = FALSE, y = d$y))

  for (s in settings) {
    dd <- d
    dd$y <- s$y

    ctrl <- function(exact) {
      bartisan_control(num_trees = 10, num_burn = 0L, num_draws = 20L,
                       gate = s$gate, augment = s$augment, verbose = FALSE,
                      exact_quadratic = exact)
    }

    set.seed(5)
    general <- bartisan(y ~ ., data = dd, family = s$family,
                        control = ctrl(FALSE))
    set.seed(5)
    closed <- bartisan(y ~ ., data = dd, family = s$family, control = ctrl(TRUE))

    for (h in seq_along(general[["eta"]])) {
      # Agreement to rounding error, not to a tolerance: the closed form is an
      # identity, so the only difference is the order the arithmetic runs in.
      expect_lt(max(abs(general[["eta"]][[h]] - closed[["eta"]][[h]])), 1e-10,
                label = paste(s$label, "predictor", h))
    }

    expect_lt(max(abs(general[["loglik"]] - closed[["loglik"]])), 1e-8,
              label = paste(s$label, "log likelihood"))
  }
})

test_that("a family that is not quadratic is untouched by the shortcut", {
  d <- sim_x(n = 150, seed = 102)
  d$y <- stats::rpois(150, exp(1 + d$x1))

  ctrl <- function(exact) {
    bartisan_control(num_trees = 10, num_burn = 0L, num_draws = 30L,
                     verbose = FALSE, exact_quadratic = exact)
  }

  # The Poisson log density is not quadratic in the predictor, so there is
  # nothing for the shortcut to do and the draws must be bitwise identical.
  set.seed(5)
  a <- bartisan(y ~ ., data = d, family = poisson(), control = ctrl(FALSE))
  set.seed(5)
  b <- bartisan(y ~ ., data = d, family = poisson(), control = ctrl(TRUE))

  expect_identical(a[["eta"]], b[["eta"]])
  expect_identical(a[["loglik"]], b[["loglik"]])
})

test_that("is_quadratic is declared per predictor, not per family", {
  # A location-scale model is quadratic in its mean and not in its log standard
  # deviation, so the shortcut has to apply to the first forest and not the
  # second. Both are exercised above; here the check is that the fit is still
  # right, since a shortcut wrongly applied to the scale would distort it.
  skip_on_cran()

  set.seed(103)
  n <- 500
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  mean_true <- 2 * d$x1
  log_sd_true <- -1 + d$x2
  d$y <- mean_true + stats::rnorm(n, sd = exp(log_sd_true))

  fit <- bartisan(y ~ ., data = d, family = location_scale(),
                  control = bartisan_control(num_trees = 20, num_burn = 300,
                                             num_draws = 300, verbose = FALSE))

  eta <- predict(fit, type = "link")
  expect_gt(stats::cor(eta[, "mean"], mean_true), 0.9)
  expect_gt(stats::cor(eta[, "log_sd"], log_sd_true), 0.6)
})
