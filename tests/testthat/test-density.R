# The conditional density reuses the family's C++ log density, so these tests
# check it against R's own d* functions. Agreement to machine precision is the
# point: it confirms both the density itself and the normalizing constants that
# the sampler is free to drop.

test_that("the Gaussian density matches dnorm", {
  d <- sim_x(seed = 41)
  n <- nrow(d)
  d$y <- 2 * d$x1 + stats::rnorm(n, sd = 0.5)

  fit <- bartisan(y ~ ., data = d, family = gaussian(),
                  control = quick_control())

  dens <- predict(fit, type = "density", draws = TRUE)
  eta <- predict(fit, type = "link", draws = TRUE)
  sigma <- fit[["aux"]][, "sigma"]

  manual <- vapply(seq_len(n), function(i) {
    stats::dnorm(d$y[i], eta[, i], sigma)
  }, numeric(nrow(eta)))

  expect_equal(dens, manual, ignore_attr = TRUE)
})

test_that("the Poisson density matches dpois, factorial term included", {
  d <- sim_x(seed = 42)
  n <- nrow(d)
  d$y <- stats::rpois(n, exp(1 + d$x1))

  fit <- bartisan(y ~ ., data = d, family = poisson(),
                  control = quick_control())

  dens <- predict(fit, type = "density", draws = TRUE)
  eta <- predict(fit, type = "link", draws = TRUE)

  manual <- vapply(seq_len(n), function(i) {
    stats::dpois(d$y[i], exp(eta[, i]))
  }, numeric(nrow(eta)))

  expect_equal(dens, manual, ignore_attr = TRUE)
})

test_that("the binomial density includes the binomial coefficient", {
  d <- sim_x(seed = 43)
  n <- nrow(d)
  trials <- stats::rpois(n, 8) + 2L
  d$succ <- stats::rbinom(n, trials, stats::plogis(d$x1))
  d$fail <- trials - d$succ

  fit <- bartisan(cbind(succ, fail) ~ x1 + x2 + x3, data = d,
                  family = binomial(), control = quick_control())

  dens <- predict(fit, type = "density", draws = TRUE)
  eta <- predict(fit, type = "link", draws = TRUE)

  manual <- vapply(seq_len(n), function(i) {
    stats::dbinom(d$succ[i], trials[i], stats::plogis(eta[, i]))
  }, numeric(nrow(eta)))

  # Without the coefficient these would differ by choose(trials, succ), so this
  # is a real check on log_norm_const() and not just on the kernel.
  expect_equal(dens, manual, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(max(trials), 1)))
})

test_that("the gamma and negative binomial densities match theirs", {
  d <- sim_x(seed = 44)
  n <- nrow(d)

  d$y <- stats::rgamma(n, 4, rate = 4 / exp(1 + d$x1))
  fit <- bartisan(y ~ ., data = d, family = stats::Gamma("log"),
                  control = quick_control())
  dens <- predict(fit, type = "density", draws = TRUE)
  eta <- predict(fit, type = "link", draws = TRUE)
  shape <- fit[["aux"]][, "shape"]
  manual <- vapply(seq_len(n), function(i) {
    stats::dgamma(d$y[i], shape = shape, rate = shape / exp(eta[, i]))
  }, numeric(nrow(eta)))
  expect_equal(dens, manual, ignore_attr = TRUE)

  d$y <- stats::rnbinom(n, mu = exp(1 + d$x1), size = 3)
  fit <- bartisan(y ~ ., data = d, family = negbin(),
                  control = quick_control())
  dens <- predict(fit, type = "density", draws = TRUE)
  eta <- predict(fit, type = "link", draws = TRUE)
  theta <- fit[["aux"]][, "theta"]
  manual <- vapply(seq_len(n), function(i) {
    stats::dnbinom(d$y[i], mu = exp(eta[, i]), size = theta)
  }, numeric(nrow(eta)))
  expect_equal(dens, manual, ignore_attr = TRUE)
})

test_that("for a categorical response the density is the observed category's probability", {
  d <- sim_x(seed = 45)
  n <- nrow(d)
  d$y <- factor(sample(1:3, n, replace = TRUE), levels = 1:3, ordered = TRUE)

  fit <- bartisan(y ~ ., data = d, family = ordinal(), control = quick_control())

  probs <- predict(fit, type = "prob")
  picked <- probs[cbind(seq_len(n), as.integer(d$y))]

  expect_equal(predict(fit, type = "density"), picked, ignore_attr = TRUE)
})

test_that("the density on new data agrees with the density on the same rows", {
  d <- sim_x(seed = 46)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- bartisan(y ~ ., data = d, control = quick_control())

  expect_equal(predict(fit, newdata = d[1:10, ], type = "density"),
               predict(fit, type = "density")[1:10])

  # `log` is a pure transform of the returned value: the average over draws is
  # taken on the density scale either way.
  expect_equal(predict(fit, type = "density", log = TRUE),
               log(predict(fit, type = "density")))

  # And the summary really is the pointwise predictive density, not the mean of
  # the log density, which Jensen's inequality puts strictly below it.
  per_draw <- predict(fit, type = "density", draws = TRUE, log = TRUE)
  expect_gt(min(predict(fit, type = "density", log = TRUE) -
                  colMeans(per_draw)), 0)
})

test_that("the density needs the outcome in new data", {
  d <- sim_x(seed = 47)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- bartisan(y ~ ., data = d, control = quick_control())

  expect_error(predict(fit, newdata = d[, c("x1", "x2", "x3")],
                       type = "density"),
               "must contain the outcome")
})

test_that("a censored survival time contributes a survival probability", {
  d <- sim_x(seed = 48)
  n <- nrow(d)
  log_t <- 2 * d$x1 + stats::rnorm(n)
  cens <- stats::quantile(log_t, 0.7)
  d$time <- exp(pmin(log_t, cens))
  d$status <- as.numeric(log_t <= cens)

  fit <- bartisan(survival::Surv(time, status) ~ x1 + x2 + x3, data = d,
                  family = lognormal_aft(), control = quick_control())

  dens <- predict(fit, type = "density", draws = TRUE)
  eta <- predict(fit, type = "link", draws = TRUE)
  sigma <- fit[["aux"]][, "sigma"]

  # Censored rows contribute S(t), which is a probability, so it is bounded by
  # one; uncensored rows contribute a density, which need not be.
  censored <- d$status == 0
  expect_true(all(dens[, censored] <= 1))

  manual <- vapply(which(censored), function(i) {
    stats::pnorm(log(d$time[i]), eta[, i], sigma, lower.tail = FALSE)
  }, numeric(nrow(eta)))
  expect_equal(dens[, censored], manual, ignore_attr = TRUE)
})

test_that("a held-out log score prefers the informative model", {
  skip_on_cran()

  set.seed(49)
  n <- 200
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- stats::rpois(n, exp(0.5 + 2 * d$x1))

  train <- d[1:140, ]
  test <- d[141:n, ]

  ctrl <- bartisan_control(num_trees = 10, num_burn = 200, num_save = 200,
                           verbose = FALSE)

  informative <- bartisan(y ~ ., data = train, family = poisson(),
                          control = ctrl)
  scrambled_data <- train
  scrambled_data$x1 <- sample(scrambled_data$x1)
  scrambled <- bartisan(y ~ ., data = scrambled_data, family = poisson(),
                        control = ctrl)

  score <- function(fit) {
    sum(predict(fit, newdata = test, type = "density", log = TRUE))
  }

  expect_gt(score(informative), score(scrambled))
})

test_that("an undefined density warns rather than returning NaN silently", {
  # A link whose inverse does not cover the line, pushed off the training range
  # so that the forest extrapolates. This is the only way a *saved* draw can be
  # out of support: `bartisan()` rejects such proposals while sampling, so the
  # fit itself is valid and the failure appears only at a new `x`.
  d <- data.frame(x1 = stats::runif(400, -1, 1),
                  x2 = stats::runif(400, -1, 1))
  set.seed(881)
  d$y <- stats::rpois(nrow(d), pmax(2 + 3 * d$x1, 0.05))

  fit <- suppressMessages(
    bartisan(y ~ ., d, family = stats::poisson("identity"),
             control = quick_control(num_trees = 20L, num_burn = 200L,
                                     num_save = 200L)))

  # In sample the predictor stays where the sampler kept it valid.
  expect_silent(in_sample <- predict(fit, type = "density"))
  expect_false(anyNA(in_sample))

  nd <- data.frame(x1 = seq(-3, 3, length.out = 40), x2 = 0, y = 1)

  expect_warning(out <- predict(fit, newdata = nd, type = "density"),
                 "undefined")

  # The value stays NaN rather than becoming zero: zero would assert the outcome
  # is impossible, and would turn into -Inf in a log score.
  expect_true(anyNA(out))
  expect_false(any(out[!is.na(out)] == 0))

  # The warning reports the amplification -- a small fraction of draws makes a
  # much larger fraction of the returned values NaN, because the draws are
  # averaged before the log is taken.
  drawn <- suppressWarnings(predict(fit, newdata = nd, type = "density",
                                    draws = TRUE))
  expect_lt(mean(is.na(drawn)), mean(is.na(out)))
})
