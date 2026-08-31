test_that("variable_importance() reports usage and separates signal from noise", {
  d <- sim_x(n = 300, p = 4, seed = 771)
  set.seed(7711)
  d$y <- 2 * d$x1 + sin(3 * d$x2) + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(num_trees = 20L, num_burn = 200L,
                                          num_draws = 200L, sparsity = TRUE))

  vi <- variable_importance(fit)

  expect_s3_class(vi, "data.frame")
  expect_identical(names(vi), c("variable", "splits", "splits_lower",
                               "splits_upper", "prop_used"))
  expect_setequal(vi$variable, c("x1", "x2", "x3", "x4"))

  # Sorted by prop_used then splits, both decreasing.
  expect_false(is.unsorted(rev(vi$prop_used)))

  # The two predictors in the truth are used in far more draws than the two
  # that are not. The gap is what makes this usable as a selection rule.
  used <- stats::setNames(vi$prop_used, vi$variable)
  expect_gt(min(used[c("x1", "x2")]), max(used[c("x3", "x4")]))

  expect_true(all(vi$splits_lower <= vi$splits))
  expect_true(all(vi$splits <= vi$splits_upper))
  expect_true(all(vi$prop_used >= 0 & vi$prop_used <= 1))

  # The same numbers summary() prints.
  usage <- summary(fit)$usage$eta
  expect_equal(sort(vi$prop_used), sort(unname(usage[, "prop_used"])))
})

test_that("variable_importance() names the forest when there is more than one", {
  d <- sim_x(n = 200, p = 3, seed = 772)
  set.seed(7721)
  d$y <- d$x1 + stats::rnorm(nrow(d), sd = exp(-1 + d$x2))

  fit <- bartisan(y ~ ., d, family = location_scale(),
                  control = quick_control(num_trees = 10L))

  vi <- variable_importance(fit)

  expect_true("predictor" %in% names(vi))
  expect_length(unique(vi$predictor), length(fit$counts))
  expect_identical(nrow(vi), length(fit$counts) * 3L)
})

test_that("variable_importance() rejects things that are not fits", {
  expect_error(variable_importance(1), "must be a fit")
})
