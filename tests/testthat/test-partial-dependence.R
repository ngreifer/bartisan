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

test_that("a predictor named through `$` counts as one variable", {
  # `all.vars(~ d$x1)` gives two names, so a single predictor read as two both
  # failed to match a column and spent half the two-variable budget.
  expect_identical(pd_variables(~ x1), "x1")
  expect_identical(pd_variables(~ d$x1), "d$x1")
  expect_identical(pd_variables(~ d$x1 + d$x2), c("d$x1", "d$x2"))

  # Three is still three, however they are written.
  expect_error(pd_variables(~ d$x1 + d$x2 + d$x3), "one or two predictors")
})

# Tree skipping: the predictor is a sum over trees, so a tree that never splits
# on a plotted column returns the same value at every grid point and is
# evaluated once. The optimization is invisible if it is right, so what these
# check is that it agrees with the path it replaces, on every shape of fit where
# the two could come apart.
#
# `pd_tree_mask()` returning `NULL` is the fallback, so mocking it out is how the
# ordinary path is reached for the comparison.
expect_same_as_unmasked <- function(fit, variables, ..., tolerance = 1e-10) {
  fast <- partial_dependence(fit, variables, ...)

  slow <- testthat::with_mocked_bindings(
    partial_dependence(fit, variables, ...),
    pd_tree_mask = function(object, vars) NULL
  )

  testthat::expect_equal(fast$estimate, slow$estimate, tolerance = tolerance)
  testthat::expect_equal(fast$lower, slow$lower, tolerance = tolerance)
  testthat::expect_equal(fast$upper, slow$upper, tolerance = tolerance)
  invisible(fast)
}

test_that("skipping trees does not change the curve", {
  d <- sim_x(n = 150L, p = 4L)
  d$g <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  # `x1` has to carry enough signal that the forest actually splits on it, or
  # the comparison below holds between two identical flat curves.
  d$y <- 6 * sin(pi * d$x1) + 3 * d$x2 + as.numeric(d$g) + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ ., d, family = stats::gaussian(),
             control = quick_control(num_trees = 20L, num_draws = 60L))))

  # A numeric predictor, a factor, two at once, and the link scale.
  expect_same_as_unmasked(fit, ~ x1, grid = 8L)
  expect_same_as_unmasked(fit, ~ g)
  expect_same_as_unmasked(fit, ~ x1 + g, grid = 5L)
  expect_same_as_unmasked(fit, ~ x2, grid = 6L, type = "link")

  # Some trees must actually have been skipped, or the comparison is vacuous.
  uses <- pd_tree_mask(fit, "x1")
  expect_false(is_null(uses))
  expect_true(any(uses))
  expect_false(all(uses))
})

test_that("the mask marks a term that mentions the variable, however written", {
  d <- sim_x(n = 120L, p = 3L)
  d$y <- 4 * d$x1 + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ log(x1 + 1) + x2 + x3, d, family = stats::gaussian(),
             control = quick_control(num_trees = 20L, num_draws = 60L))))

  # The term is `log(x1 + 1)` and the variable is `x1`, so a mask keyed on the
  # label alone would miss it and hold the curve flat with no sign of doing so.
  expect_true(any(pd_tree_mask(fit, "x1")))

  # The model frame holds `log(x1 + 1)` rather than `x1`, so the grid over the
  # untransformed variable needs the original data to average over.
  expect_same_as_unmasked(fit, ~ x1, newdata = d, grid = 6L)
})

test_that("a grouping factor falls back rather than freezing its intercepts", {
  d <- sim_x(n = 120L, p = 2L)
  d$g <- factor(rep(letters[1:4], length.out = nrow(d)))
  d$y <- d$x1 + as.numeric(d$g) + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ x1 + x2 + (1 | g), d, family = stats::gaussian(),
             control = quick_control())))

  # A `(1 | g)` intercept is a drawn parameter rather than a tree, so no choice
  # of trees holds it fixed while the grid moves `g`.
  expect_null(pd_tree_mask(fit, "g"))

  # A predictor that is not the grouping factor is still optimized, and the
  # random part rides along in the constant base.
  expect_false(is_null(pd_tree_mask(fit, "x1")))
  expect_same_as_unmasked(fit, ~ x1, grid = 6L)
})

test_that("a varying coefficient is combined after the trees are summed", {
  d <- sim_x(n = 120L, p = 3L)
  d$z <- stats::runif(nrow(d))
  d$y <- d$x1 + d$z * (1 + d$x2) + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ x1 + x2 + vc(z), d, family = stats::gaussian(),
             control = quick_control())))

  # `z` multiplies its own forest, so the part that does not move still lands in
  # the prediction multiplied by a value the grid changes. Combining after the
  # sum is what keeps that right.
  expect_same_as_unmasked(fit, ~ z, grid = 6L)
  expect_same_as_unmasked(fit, ~ x2, grid = 6L)
})

test_that("a bcf propensity score is treated as moving with the grid", {
  d <- sim_x(n = 200L, p = 3L)
  d$z <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1 - 0.5))
  d$y <- d$x1 + 0.5 * d$z + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bcf(y ~ x1 + x2 + x3, treat = ~ z, data = d, family = stats::gaussian(),
        num_trees = 10L, num_burn = 30L, num_draws = 40L, verbose = FALSE)))

  # The score is a function of the covariates and is rebuilt from the data at
  # every grid point, so trees splitting on it cannot go in the fixed base. It
  # is marked as moving whether or not the plotted variable feeds it, which is
  # the cheap side of that trade.
  labels <- fit[["term_labels"]]
  score_group <- match(colnames(fit[["bcf"]][["propensity"]]), labels)
  expect_false(is.na(score_group))

  # A variable that names no term at all still leaves the score group moving, so
  # its mask is the one the score's own columns give. `any()` of that mask would
  # instead be asking whether a 10-tree chain happened to split on the score,
  # which it need not.
  by_score <- pd_tree_mask(fit, colnames(fit[["bcf"]][["propensity"]]))
  nothing <- pd_tree_mask(fit, "nothing_at_all")

  expect_false(is.null(nothing))
  expect_identical(nothing, by_score)

  expect_same_as_unmasked(fit, ~ x1, grid = 5L)
  expect_same_as_unmasked(fit, ~ x2, grid = 5L)
})

test_that("a variable no tree splits on gives a flat curve at no cost", {
  d <- sim_x(n = 150L, p = 4L)
  d$g <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  # The factor dominates a short chain, so the forest never splits on `x1`.
  d$y <- 2 * sin(pi * d$x1) + d$x2 + as.numeric(d$g) + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ ., d, family = stats::gaussian(), control = quick_control())))

  uses <- pd_tree_mask(fit, "x1")

  skip_if(any(uses), "this chain did split on x1, so there is nothing to check")

  # Every tree lands in the base, the grid costs one forest evaluation in total,
  # and the curve is exactly flat rather than flat to rounding.
  pd <- partial_dependence(fit, ~ x1, grid = 6L)
  expect_identical(diff(range(pd[["estimate"]])), 0)
  expect_same_as_unmasked(fit, ~ x1, grid = 6L)
})

test_that("an argument the fast path does not know hands the grid to predict()", {
  d <- sim_x(n = 120L, p = 3L)
  d$y <- 4 * d$x1 + stats::rnorm(nrow(d))

  fit <- suppressMessages(suppressWarnings(
    bartisan(y ~ ., d, family = stats::gaussian(),
             control = quick_control(num_trees = 20L, num_draws = 60L))))

  taken <- function(...) {
    !identical(body(pd_predictor(fit, d, "x1", "response", list(...))),
               body(pd_predictor(fit, d, "x1", "nonesuch", list())))
  }

  # The names it reads itself keep the fast path.
  expect_true(taken(iterations = 1:10))
  expect_true(taken(offset = NULL, log = FALSE))
  expect_true(taken())

  # Anything else, and anything positional, goes back to `predict()`, which is
  # what knows where such an argument belongs.
  expect_false(taken(nonesuch = 1))
  expect_false(taken(iterations = 1:10, nonesuch = 1))
  expect_false(taken(1))

  # And the fallback still gives the same curve, which is the point of taking it.
  expect_equal(partial_dependence(fit, ~ x1, grid = 5L, iterations = 1:20)$estimate,
               partial_dependence(fit, ~ x1, grid = 5L, iterations = 1:20,
                                  nonesuch = 1)$estimate)
})
