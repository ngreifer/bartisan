test_that("genbart_control validates its arguments", {
  expect_s3_class(genbart_control(), "genbart_control")

  expect_error(genbart_control(num_trees = 0), "must be")
  expect_error(genbart_control(num_trees = 2.5), "whole number")
  expect_error(genbart_control(soft = "yes"), "must be")
  expect_error(genbart_control(bandwidth = -1), "must be")
  expect_error(genbart_control(gamma = 1.5), "must be")
  expect_error(genbart_control(sigma_mu = -1), "must be positive")
  expect_error(genbart_control(sigma_mu_ramp = 2), "must be between")
  expect_error(genbart_control(x_transform = "nope"), "should be one of")
})

test_that("genbart rejects a control that did not come from genbart_control", {
  d <- sim_x(seed = 31)
  d$y <- stats::rnorm(nrow(d))

  expect_error(genbart(y ~ ., data = d, control = list(num_trees = 5)),
               "must be the result of")
})

test_that("the leaf scale defaults to the documented value", {
  d <- sim_x(seed = 32)
  d$y <- 5 * d$x1 + stats::rnorm(nrow(d))

  ctrl <- quick_control(num_trees = 8L, k = 2, update_sigma_mu = FALSE,
                        sigma_mu_ramp = 0)
  fit <- genbart(y ~ ., data = d, control = ctrl)

  expected <- 3 * stats::sd(d$y) / (2 * sqrt(8))
  expect_equal(unname(fit[["sigma_mu"]][1L, 1L]), expected, tolerance = 1e-8)
})

test_that("both predictor transforms map into the unit interval", {
  d <- sim_x(seed = 33)
  d$x2 <- d$x2 * 1000 + 500
  d$y <- stats::rnorm(nrow(d))

  for (type in c("quantile", "range")) {
    fit <- genbart(y ~ ., data = d, control = quick_control(x_transform = type))
    mapped <- vapply(seq_along(fit[["unit_maps"]]), function(j) {
      z <- fit[["unit_maps"]][[j]](d[[j]])
      all(z >= 0 & z <= 1)
    }, logical(1L))
    expect_true(all(mapped))
  }
})

test_that("a binary predictor maps to zero and one under either transform", {
  for (type in c("quantile", "range")) {
    f <- make_unit_map(c(0, 1, 1, 0), type)
    expect_identical(f(c(0, 1)), c(0, 1))
  }
})

test_that("a constant predictor maps to a single interior value", {
  f <- make_unit_map(rep(3, 10), "quantile")
  expect_identical(f(c(3, 3)), c(0.5, 0.5))
})

test_that("the sparsity prior is disabled when there is only one group", {
  d <- sim_x(p = 1, seed = 34)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- genbart(y ~ x1, data = d, control = quick_control())

  expect_identical(ncol(fit[["counts"]][["eta"]]), 1L)
  expect_s3_class(fit, "genbart")
})

test_that("turning off the ramp leaves the leaf scale free to move", {
  d <- sim_x(seed = 35)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- genbart(y ~ ., data = d,
                 control = quick_control(sigma_mu_ramp = 0,
                                         update_sigma_mu = TRUE))

  expect_true(stats::sd(fit[["sigma_mu"]][, 1L]) > 0)
})

test_that("control arguments can be passed through genbart()'s dots", {
  d <- sim_x(seed = 41)
  d$y <- stats::rnorm(nrow(d))

  # A setting passed in `...` reaches the control list, and the ones passed to
  # genbart_control() alongside it survive.
  fit <- genbart(y ~ ., data = d, control = quick_control(),
                 num_trees = 3, soft = FALSE)

  expect_identical(fit[["control"]][["num_trees"]], 3L)
  expect_false(fit[["control"]][["soft"]])
  expect_identical(fit[["control"]][["num_burn"]], 30L)

  # `...` wins over the same setting given in `control`.
  fit2 <- genbart(y ~ ., data = d, control = quick_control(num_trees = 9),
                  num_trees = 4)

  expect_identical(fit2[["control"]][["num_trees"]], 4L)

  # Passing a setting through `...` is the same call as passing it through
  # `control`: same draws from the same seed.
  set.seed(1)
  a <- genbart(y ~ ., data = d, control = quick_control(), num_trees = 3,
               soft = FALSE)
  set.seed(1)
  b <- genbart(y ~ ., data = d,
               control = quick_control(num_trees = 3, soft = FALSE))

  expect_equal(a[["eta"]], b[["eta"]])

  # An argument whose default is NULL is kept when passed explicitly, rather
  # than being dropped from the merged list.
  fit3 <- genbart(y ~ ., data = d, control = quick_control(),
                  sigma_mu = NULL)
  expect_null(fit3[["control"]][["sigma_mu"]])
})

test_that("unknown dots arguments to genbart() are an error", {
  d <- sim_x(seed = 42)
  d$y <- stats::rnorm(nrow(d))

  expect_error(genbart(y ~ ., data = d, control = quick_control(),
                       num_tree = 3),
               "not an argument")
  expect_error(genbart(y ~ ., data = d, control = quick_control(),
                       nonsense = 1, other = 2),
               "not arguments")

  # Validation still happens in genbart_control(), rather than being skipped
  # for settings that arrive this way.
  expect_error(genbart(y ~ ., data = d, control = quick_control(),
                       num_trees = -1),
               "num_trees")
})
