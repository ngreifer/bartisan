# Partial dependence: the fitted surface averaged over the sample at each value
# of one or two predictors.

test_that("one numeric predictor gives a curve on its own grid", {
  d <- sim_x(n = 120L, p = 3L)
  d$y <- 2 * sin(pi * d$x1) + d$x2 + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ ., d, family = stats::gaussian(), control = quick_control())))

  pd <- partial_dependence(fit, ~ x1, grid = 10L)

  expect_s3_class(pd, "bartisan_partial")
  expect_identical(nrow(pd), 10L)
  expect_true(all(c("x1", "estimate", "lower", "upper") %in% names(pd)))
  expect_true(all(pd[["lower"]] <= pd[["estimate"]]))
  expect_true(all(pd[["estimate"]] <= pd[["upper"]]))
  expect_identical(attr(pd, "variables"), "x1")

  # The grid spans the predictor rather than inventing values beyond it.
  expect_equal(range(pd[["x1"]]), range(d[["x1"]]), tolerance = 1e-8)
})

test_that("a factor is evaluated at its levels whatever the grid says", {
  d <- sim_x(n = 120L, p = 2L)
  d$g <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  d$y <- d$x1 + as.integer(d$g) + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ x1 + x2 + g, d, family = stats::gaussian(),
             control = quick_control())))

  pd <- partial_dependence(fit, ~ g, grid = 25L)

  expect_identical(nrow(pd), 3L)
  expect_setequal(pd[["g"]], c("a", "b", "c"))
})

test_that("two predictors give the product of their grids", {
  d <- sim_x(n = 100L, p = 2L)
  d$g <- factor(sample(c("a", "b"), nrow(d), replace = TRUE))
  d$y <- d$x1 + as.integer(d$g) * d$x2 + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ x1 + x2 + g, d, family = stats::gaussian(),
             control = quick_control())))

  pd <- partial_dependence(fit, ~ x1 + g, grid = 6L)

  expect_identical(nrow(pd), 12L)
  expect_identical(attr(pd, "variables"), c("x1", "g"))
})

test_that("the average is taken within a draw, so the interval is on the average", {
  d <- sim_x(n = 100L, p = 2L)
  d$y <- d$x1 + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ ., d, family = stats::gaussian(),
             control = quick_control(num_draws = 100L))))

  pd <- partial_dependence(fit, ~ x1, grid = 3L)

  # Recomputed by hand at one grid point: set the column, predict every unit,
  # average across units within each draw, then summarize those averages. An
  # interval built the other way round, over units, would be far wider.
  nd <- fit[["model"]]
  nd[["x1"]] <- pd[["x1"]][2L]
  draws <- stats::predict(fit, newdata = nd, type = "response", draws = TRUE)
  by_hand <- post_summary(rowMeans(draws), level = 0.95)

  expect_equal(pd[["estimate"]][2L], by_hand[["mean"]], tolerance = 1e-8)
  expect_equal(pd[["lower"]][2L], by_hand[["lower"]], tolerance = 1e-8)
})

test_that("values overrides the grid for the predictors it names", {
  d <- sim_x(n = 80L, p = 2L)
  d$y <- d$x1 + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ ., d, family = stats::gaussian(), control = quick_control())))

  pd <- partial_dependence(fit, ~ x1, values = list(x1 = c(0.1, 0.9)))

  expect_identical(nrow(pd), 2L)
  expect_equal(pd[["x1"]], c(0.1, 0.9))
})

test_that("three predictors and unknown ones are refused", {
  d <- sim_x(n = 80L, p = 3L)
  d$y <- d$x1 + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ ., d, family = stats::gaussian(), control = quick_control())))

  expect_error(partial_dependence(fit, ~ x1 + x2 + x3),
               "one or two predictors")
  expect_error(partial_dependence(fit, ~ nope), "does not have")
})

test_that("the plot argument and the plot method draw the same thing", {
  skip_if_not_installed("ggplot2")

  d <- sim_x(n = 80L, p = 2L)
  d$g <- factor(sample(c("a", "b"), nrow(d), replace = TRUE))
  d$y <- d$x1 + as.integer(d$g) + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ x1 + x2 + g, d, family = stats::gaussian(),
             control = quick_control())))

  for (v in list(~ x1, ~ g, ~ x1 + g)) {
    pd <- partial_dependence(fit, v, grid = 5L)
    from_method <- plot(pd)
    from_arg <- partial_dependence(fit, v, grid = 5L, plot = TRUE)

    expect_s3_class(from_arg, "ggplot")
    expect_equal(from_arg[["data"]], from_method[["data"]])
    expect_equal(from_arg[["labels"]], from_method[["labels"]])
  }

  # `plot()` on the fit is the same drawing, which is what makes it sugar.
  expect_equal(plot(fit, ~ x1, grid = 5L)[["data"]],
               plot(partial_dependence(fit, ~ x1, grid = 5L))[["data"]])

  # And it says what it needs when told nothing.
  expect_error(plot(fit), "told which predictors")
})

# The grid loop predicts once per point over the whole sample, so it goes to
# workers when a plan has any. The streams are drawn before the branch, which is
# what keeps a family whose prediction simulates from depending on whether a
# plan was set.
test_that("the grid is the same in parallel as in sequence", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  d <- sim_x(n = 200L, p = 3L, seed = 97L)
  d$g <- factor(sample(c("a", "b"), nrow(d), replace = TRUE))
  d$y <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1))

  fit <- bartisan(y ~ ., d, family = stats::binomial(),
                  control = quick_control(num_trees = 5L, num_burn = 60L,
                                          num_draws = 100L))

  future::plan(future::sequential)
  one <- partial_dependence(fit, ~ x1, grid = 6L)

  started <- tryCatch({
    future::plan(future::multisession, workers = 2L)
    isTRUE(future::nbrOfWorkers() >= 2L)
  }, error = function(e) FALSE)

  skip_if_not(started, "no second worker available")

  many <- partial_dependence(fit, ~ x1, grid = 6L)

  expect_identical(nrow(one), 6L)
  expect_equal(as.data.frame(many), as.data.frame(one))

  # Two predictors give one row per combination, and the same must hold there.
  pair_seq <- local({
    future::plan(future::sequential)
    partial_dependence(fit, ~ x1 + g, grid = 4L)
  })

  future::plan(future::multisession, workers = 2L)
  pair_par <- partial_dependence(fit, ~ x1 + g, grid = 4L)

  expect_identical(nrow(pair_seq), 8L)
  expect_equal(as.data.frame(pair_par), as.data.frame(pair_seq))
})
