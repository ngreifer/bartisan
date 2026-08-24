test_that("the zero-inflated constructors resolve", {
  expect_identical(as_genbart_family(zi_poisson())[["family"]], "zip")
  expect_identical(as_genbart_family(zi_negbin())[["family"]], "zinb")
  expect_identical(as_genbart_family("zi_poisson")[["family"]], "zip")
  expect_identical(as_genbart_family(ordbeta())[["family"]], "ordbeta")

  expect_error(zi_poisson("logit"), "should be")
  expect_error(ordbeta("probit"), "should be")
  expect_error(zi_negbin(theta = -1), "must be")
  expect_error(ordbeta(phi = 0), "must be")
})

test_that("zero-inflated models carry a forest for each component", {
  d <- sim_x(seed = 51)
  n <- nrow(d)
  structural <- stats::rbinom(n, 1, 0.3)
  d$y <- ifelse(structural == 1, 0, stats::rpois(n, exp(1 + d$x1)))

  fit <- genbart(y ~ ., data = d, family = zi_poisson(),
                 control = quick_control())

  expect_identical(fit[["num_forest"]], 2L)
  expect_identical(names(fit[["eta"]]), c("count", "zero"))
  expect_null(fit[["aux"]])
  expect_predictor_invariant(fit, d)

  # The fitted mean is the count mean scaled down by the inflation probability,
  # so it must sit below the count mean itself.
  link <- predict(fit, type = "link")
  expect_true(all(predict(fit, type = "response") < exp(link[, "count"]) + 1e-8))

  nb <- genbart(y ~ ., data = d, family = zi_negbin(),
                control = quick_control())
  expect_identical(colnames(nb[["aux"]]), "theta")
  expect_true(all(nb[["aux"]][, "theta"] > 0))
  expect_predictor_invariant(nb, d)
})

test_that("a fixed dispersion is held fixed", {
  d <- sim_x(seed = 52)
  d$y <- stats::rnbinom(nrow(d), mu = exp(1 + d$x1), size = 2)

  fit <- genbart(y ~ ., data = d, family = zi_negbin(theta = 2.5),
                 control = quick_control())

  expect_true(all(fit[["aux"]][, "theta"] == 2.5))
})

test_that("ordered beta fits a response with mass at both endpoints", {
  d <- sim_x(seed = 53)
  n <- nrow(d)
  eta <- 2 * (d$x1 - 0.5)
  p_one <- stats::plogis(eta - 1.5)
  p_zero <- 1 - stats::plogis(eta + 1.5)
  u <- stats::runif(n)
  y <- stats::rbeta(n, stats::plogis(eta) * 8, (1 - stats::plogis(eta)) * 8)
  y[u <= p_zero] <- 0
  y[u >= 1 - p_one] <- 1
  d$y <- y

  fit <- genbart(y ~ ., data = d, family = ordbeta(),
                 control = quick_control())

  expect_identical(fit[["num_forest"]], 1L)
  expect_identical(colnames(fit[["aux"]]), c("cut1", "cut2", "phi"))
  # Both cutpoints are drawn, and their order is enforced by construction.
  expect_true(all(fit[["aux"]][, "cut1"] < fit[["aux"]][, "cut2"]))
  expect_true(all(fit[["aux"]][, "phi"] > 0))
  expect_predictor_invariant(fit, d)

  # The fitted mean respects the bounds of the response.
  fitted <- predict(fit, type = "response")
  expect_true(all(fitted >= 0 & fitted <= 1))
})

test_that("ordered beta rejects responses it cannot model", {
  d <- sim_x(seed = 54)
  n <- nrow(d)

  d$y <- stats::runif(n, -1, 1)
  expect_error(genbart(y ~ ., data = d, family = ordbeta(),
                       control = quick_control()),
               "between 0 and 1")

  d$y <- stats::rbinom(n, 1, 0.5)
  expect_error(genbart(y ~ ., data = d, family = ordbeta(),
                       control = quick_control()),
               "needs some responses strictly")
})

test_that("zero-inflated recovery: both components are found", {
  skip_on_cran()

  set.seed(55)
  n <- 600
  X <- matrix(stats::runif(n * 4), nrow = n)
  colnames(X) <- paste0("x", 1:4)
  d <- as.data.frame(X)

  log_mu <- 1 + 1.5 * (X[, 1] - 0.5)
  p_zero <- stats::plogis(-0.5 + 2 * (X[, 2] - 0.5))
  structural <- stats::rbinom(n, 1, p_zero)
  d$y <- ifelse(structural == 1, 0, stats::rpois(n, exp(log_mu)))

  fit <- genbart(y ~ ., data = d, family = zi_poisson(),
                 control = genbart_control(num_trees = 20, num_burn = 400,
                                           num_save = 400, verbose = FALSE))

  link <- predict(fit, type = "link")
  expect_gt(stats::cor(link[, "count"], log_mu), 0.85)
  expect_gt(stats::cor(link[, "zero"], stats::qlogis(p_zero)), 0.7)

  # Each component's forest should lean on its own predictor.
  usage <- summary(fit)[["usage"]]
  expect_identical(rownames(usage[["count"]])[1L], "x1")
  expect_identical(rownames(usage[["zero"]])[1L], "x2")
})

test_that("ordered beta recovery: cutpoints and precision are found", {
  skip_on_cran()

  set.seed(56)
  n <- 600
  X <- matrix(stats::runif(n * 4), nrow = n)
  colnames(X) <- paste0("x", 1:4)
  d <- as.data.frame(X)

  eta <- 2 * (X[, 1] - 0.5) + 1.5 * (X[, 2] - 0.5)
  cut1 <- -1.5
  cut2 <- 1.5
  phi <- 8
  mu <- stats::plogis(eta)

  p_one <- stats::plogis(eta - cut2)
  p_zero <- 1 - stats::plogis(eta - cut1)
  u <- stats::runif(n)
  y <- stats::rbeta(n, mu * phi, (1 - mu) * phi)
  y[u <= p_zero] <- 0
  y[u >= 1 - p_one] <- 1
  d$y <- y

  fit <- genbart(y ~ ., data = d, family = ordbeta(),
                 control = genbart_control(num_trees = 20, num_burn = 400,
                                           num_save = 400, verbose = FALSE))

  expect_gt(stats::cor(predict(fit, type = "link"), eta), 0.9)
  expect_equal(mean(fit[["aux"]][, "cut1"]), cut1, tolerance = 0.35)
  expect_equal(mean(fit[["aux"]][, "cut2"]), cut2, tolerance = 0.35)
  expect_equal(mean(fit[["aux"]][, "phi"]), phi, tolerance = 0.35)

  # The endpoint masses are reproduced.
  expect_equal(mean(predict(fit, type = "response")), mean(y), tolerance = 0.05)
})
