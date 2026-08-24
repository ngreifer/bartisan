# The point of these is not that marginaleffects works -- it is that the draws
# reach it correctly. Every estimand it computes is the same transformation
# applied to every draw, so if the "posterior_draws" attribute is transposed, or
# stacked in the wrong order for a categorical family, the numbers come out
# plausible and wrong.

skip_if_no_me <- function() {
  testthat::skip_if_not_installed("marginaleffects")
}

test_that("the base generics report what insight will ask them for", {
  d <- sim_x(n = 80, seed = 91)
  d$y <- d$x1 + stats::rnorm(nrow(d))

  fit <- genbart(y ~ x1 + x2, d, control = quick_control())

  expect_s3_class(stats::formula(fit), "formula")
  expect_s3_class(stats::terms(fit), "terms")
  expect_identical(stats::nobs(fit), 80L)
  expect_s3_class(stats::model.frame(fit), "data.frame")
  # The model frame is the data the fit saw, which is what a counterfactual grid
  # is built from.
  expect_identical(nrow(stats::model.frame(fit)), 80L)
  expect_true(all(c("y", "x1", "x2") %in% names(stats::model.frame(fit))))
})

test_that("get_predict hands over draws in the orientation marginaleffects expects", {
  skip_if_no_me()

  d <- sim_x(n = 60, seed = 92)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))
  fit <- genbart(y ~ x1 + x2, d, control = quick_control())

  out <- marginaleffects::get_predict(fit)

  expect_true(all(c("rowid", "group", "estimate") %in% names(out)))
  expect_identical(nrow(out), 60L)

  draws <- attr(out, "posterior_draws")
  # Observations by draws, not the other way round. Getting this backwards is
  # the failure that produces plausible nonsense.
  expect_identical(dim(draws), c(60L, 30L))

  # And the point estimate is the row mean of those draws, so the two agree
  # about which margin is which.
  expect_equal(out$estimate, rowMeans(draws), ignore_attr = TRUE)

  # The draws are the same numbers predict() gives, transposed.
  direct <- stats::predict(fit, newdata = d, type = "response", draws = TRUE)
  expect_equal(draws, t(direct), ignore_attr = TRUE)
})

test_that("a categorical family stacks rows and draws the same way", {
  skip_if_no_me()

  d <- sim_x(n = 60, seed = 93)
  set.seed(1093)
  z <- 2 * d$x1 - d$x2 + stats::rnorm(nrow(d))
  d$y <- ordered(rowSums(outer(z, stats::quantile(z, c(1, 2) / 3), ">")) + 1L)

  fit <- genbart(y ~ x1 + x2, d, family = ordinal("probit"),
                 control = quick_control())

  out <- marginaleffects::get_predict(fit)
  draws <- attr(out, "posterior_draws")

  # Three categories, so three blocks of 60 rows.
  expect_identical(nrow(out), 180L)
  expect_identical(dim(draws), c(180L, 30L))
  expect_identical(as.character(out$group[c(1, 61, 121)]), c("1", "2", "3"))

  # Row order and draw order agree, which is the thing that would silently break.
  expect_equal(out$estimate, rowMeans(draws), ignore_attr = TRUE)

  # Probabilities over the categories sum to one for every row of every draw.
  by_cat <- array(draws, dim = c(60L, 3L, 30L))
  expect_equal(apply(by_cat, c(1L, 3L), sum), matrix(1, 60L, 30L),
               tolerance = 1e-10)
})

test_that("the estimands run, and their intervals come from the posterior", {
  skip_if_no_me()

  set.seed(94)
  n <- 300
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  z = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- 2 * d$x1 - d$x2 + (d$z == "b") + stats::rnorm(n, sd = 0.5)

  fit <- genbart(y ~ x1 + x2 + z, d, num_trees = 20L, num_burn = 200L,
                 num_save = 200L)

  p <- marginaleffects::avg_predictions(fit)
  expect_identical(nrow(p), 1L)
  expect_true(p$conf.low < p$estimate && p$estimate < p$conf.high)

  cmp <- marginaleffects::avg_comparisons(fit, variables = "z")
  expect_identical(nrow(cmp), 1L)
  # The factor contrast is a whole unit in truth, and shrinkage pulls it in, so
  # this checks the sign and the order of magnitude rather than the value.
  expect_gt(cmp$estimate, 0.3)
  expect_lt(cmp$estimate, 1.5)

  # A grouped prediction splits by the factor and orders as the factor does.
  by_z <- marginaleffects::avg_predictions(fit, by = "z")
  expect_identical(nrow(by_z), 2L)
  expect_lt(by_z$estimate[1L], by_z$estimate[2L])
})

test_that("a multi-predictor family works on the response scale and not the link one", {
  skip_if_no_me()

  d <- sim_x(n = 80, seed = 95)
  d$y <- stats::rnorm(nrow(d))

  fit <- genbart(y ~ x1 + x2, d, family = location_scale(),
                 control = quick_control())

  # The mean is one number per observation, so it goes through.
  out <- marginaleffects::get_predict(fit, type = "response")
  expect_identical(nrow(out), 80L)
  expect_identical(dim(attr(out, "posterior_draws")), c(80L, 30L))

  # The link is not, and saying so beats answering about the location and
  # leaving the scale unmentioned.
  expect_error(marginaleffects::get_predict(fit, type = "link"),
               "not a single quantity")

  # A zero-inflated family is the same shape of case.
  d$count <- stats::rpois(nrow(d), 2) * stats::rbinom(nrow(d), 1, 0.7)
  zi <- genbart(count ~ x1 + x2, d, family = zi_poisson(),
                control = quick_control())
  expect_identical(nrow(marginaleffects::get_predict(zi, type = "response")),
                   80L)
  expect_error(marginaleffects::get_predict(zi, type = "link"),
               "not a single quantity")
})

test_that("a multinomial fit gives one group per category", {
  skip_if_no_me()

  d <- sim_x(n = 60, seed = 97)
  set.seed(1097)
  d$y <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))

  fit <- genbart(y ~ x1 + x2, d, family = multinomial(),
                 control = quick_control())

  out <- marginaleffects::get_predict(fit)
  expect_identical(nrow(out), 180L)
  expect_identical(as.character(unique(out$group)), c("a", "b", "c"))
  expect_equal(out$estimate, rowMeans(attr(out, "posterior_draws")),
               ignore_attr = TRUE)
})

test_that("slopes are stable under a linear predictor transform and not under the quantile one", {
  skip_if_no_me()
  skip_on_cran()

  # The claim documented in ?genbart-marginaleffects: the default transform is an
  # empirical distribution function, so the fit is a step function of the
  # original predictor and its difference quotient diverges as the step shrinks.
  set.seed(96)
  n <- 600
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- 2 * d$x1 - d$x2 + stats::rnorm(n, sd = 0.5)

  slope_at <- function(transform, eps) {
    fit <- genbart(y ~ x1 + x2, d, x_transform = transform, num_trees = 20L,
                   num_burn = 200L, num_save = 200L)
    marginaleffects::avg_slopes(fit, variables = "x2", eps = eps)$estimate
  }

  linear_small <- slope_at("range", 1e-4)
  linear_big <- slope_at("range", 0.05)

  # Stable across a 500-fold change in the step, and of the right sign.
  expect_lt(abs(linear_small - linear_big), 0.25)
  expect_lt(linear_big, 0)

  # The quantile transform is not: at a small step the quotient is inflated by
  # roughly an order of magnitude.
  quantile_small <- slope_at("quantile", 1e-4)
  expect_gt(abs(quantile_small), 3 * abs(linear_small))
})
