test_that("variable_importance() reports usage and separates signal from noise", {
  d <- sim_x(n = 300, p = 4, seed = 771)
  set.seed(7711)
  d$y <- 2 * d$x1 + sin(3 * d$x2) + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(num_trees = 20L, num_burn = 200L,
                                          num_draws = 200L, sparsity = TRUE))

  vi <- variable_importance(fit)

  expect_s3_class(vi, "bartisan_importance")
  expect_s3_class(vi, "data.frame")
  expect_identical(names(vi), c("variable", "prop_used", "prop_splits",
                               "splits", "splits_lower", "splits_upper"))
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

  fit <- bartisan(y ~ ., d, family = gaussian_ls(),
                  control = quick_control(num_trees = 10L))

  vi <- variable_importance(fit)

  expect_true("predictor" %in% names(vi))
  expect_length(unique(vi$predictor), length(fit$counts))
  expect_identical(nrow(vi), length(fit$counts) * 3L)
})

test_that("variable_importance() rejects things that are not fits", {
  expect_error(variable_importance(1), "must be a fit")
  expect_error(variable_importance(1, level = 1), "must be a fit")
})

test_that("prop_splits is a share, and so survives a change of forest size", {
  d <- sim_x(n = 300, p = 4, seed = 773)
  set.seed(7731)
  d$y <- 2 * d$x1 + sin(3 * d$x2) + stats::rnorm(nrow(d), sd = 0.3)

  small <- bartisan(y ~ ., d, family = stats::gaussian(),
                    control = quick_control(num_trees = 10L, num_burn = 200L,
                                            num_draws = 200L, sparsity = TRUE))
  large <- bartisan(y ~ ., d, family = stats::gaussian(),
                    control = quick_control(num_trees = 40L, num_burn = 200L,
                                            num_draws = 200L, sparsity = TRUE))

  a <- variable_importance(small)
  b <- variable_importance(large)

  # The shares add to one within each fit, which is the property that makes
  # them comparable between fits.
  expect_equal(sum(a$prop_splits), 1)
  expect_equal(sum(b$prop_splits), 1)

  # Four times the trees means several times the rules, which is what makes the
  # raw count useless for comparing the two fits.
  expect_gt(sum(b$splits) / sum(a$splits), 2)

  # The share does not scale that way, and the ranking is what carries over.
  # It is not invariant: measured here, the leading predictor's share moved from
  # .51 to .64 when the forest went from 10 trees to 40, so the column removes
  # the dependence on `num_trees` without removing the forest from the answer.
  expect_true(all(a$prop_splits >= 0 & a$prop_splits <= 1))
  expect_identical(a$variable[1L], b$variable[1L])
  expect_identical(a$variable[2L], b$variable[2L])
})

test_that("draws = TRUE hands back the counts themselves", {
  d <- sim_x(n = 200, p = 3, seed = 774)
  set.seed(7741)
  d$y <- d$x1 + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(num_trees = 10L, num_draws = 40L))

  counts <- variable_importance(fit, draws = TRUE)

  expect_identical(counts, fit$counts$eta)
  expect_identical(dim(counts), c(40L, 3L))

  # The summary has to be a summary of exactly these numbers.
  vi <- variable_importance(fit)
  o <- match(vi$variable, colnames(counts))
  expect_equal(vi$splits, unname(colMeans(counts)[o]))
  expect_equal(vi$prop_used, unname(colMeans(counts > 0)[o]))

  # Several forests give one matrix each, the way `coef()` does.
  ls_fit <- bartisan(y ~ ., d, family = gaussian_ls(),
                     control = quick_control(num_trees = 10L))
  expect_type(variable_importance(ls_fit, draws = TRUE), "list")
  expect_length(variable_importance(ls_fit, draws = TRUE), length(ls_fit$counts))

  expect_error(variable_importance(fit, draws = TRUE, plot = TRUE),
               "different things")
})

test_that("the print method says which reading the fit supports", {
  d <- sim_x(n = 200, p = 3, seed = 775)
  set.seed(7751)
  d$y <- d$x1 + stats::rnorm(nrow(d), sd = 0.3)

  ctrl <- function(sparsity) {
    quick_control(num_trees = 10L, num_draws = 60L, sparsity = sparsity)
  }

  on_fit <- bartisan(y ~ ., d, family = stats::gaussian(), control = ctrl(TRUE))
  off_fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                      control = ctrl(FALSE))

  expect_output(print(variable_importance(on_fit)), "Variable importance")

  # The interval columns are in the object and out of the printed table, and
  # `sparsity = FALSE` is called out because it is what makes `prop_used`
  # unreadable as a selection rule.
  # The table goes to stdout and the notes below it are cli bullets, which are
  # conditions rather than printed output, so the two are captured separately.
  shown <- capture.output(print(variable_importance(on_fit)))
  notes <- capture_messages(print(variable_importance(on_fit)))

  expect_false(any(grepl("splits_lower", shown)))
  expect_true(any(grepl("95% interval", notes)))
  expect_false(any(grepl("sparsity = FALSE", notes)))

  notes_off <- capture_messages(print(variable_importance(off_fit)))
  expect_true(any(grepl("sparsity = FALSE", notes_off)))

  # And a subset still prints, which it cannot do by carrying the attributes.
  expect_output(print(subset(variable_importance(on_fit), prop_used > 0)),
                "Variable importance")
})

test_that("plot = TRUE draws the table it would otherwise return", {
  skip_if_not_installed("ggplot2")

  d <- sim_x(n = 200, p = 3, seed = 776)
  set.seed(7761)
  d$y <- d$x1 + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(num_trees = 10L))

  expect_s3_class(variable_importance(fit, plot = TRUE), "ggplot")

  # One panel per forest when there is more than one.
  ls_fit <- bartisan(y ~ ., d, family = gaussian_ls(),
                     control = quick_control(num_trees = 10L))
  expect_s3_class(variable_importance(ls_fit, plot = TRUE), "ggplot")
})

test_that("the plot method and the plot argument are one drawing", {
  skip_if_not_installed("ggplot2")

  d <- sim_x(n = 200, p = 4, seed = 777)
  set.seed(7771)
  d$y <- 2 * d$x1 + sin(3 * d$x2) + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(num_trees = 10L, num_draws = 100L))

  imp <- variable_importance(fit)
  by_arg <- variable_importance(fit, plot = TRUE)
  by_method <- plot(imp)

  expect_s3_class(by_method, "ggplot")

  # `plot = TRUE` calls the method, so these are the same drawing rather than
  # two of them; the ggplot objects themselves carry environments that will not
  # compare equal, so the comparison is on the pieces.
  expect_equal(by_method$data, by_arg$data)
  expect_identical(by_method$labels, by_arg$labels)

  # Subsetting first is how a wide model is made readable, so the method has to
  # accept what `head()` and `subset()` return.
  expect_identical(nrow(plot(head(imp, 2L))$data), 2L)
  expect_identical(nrow(plot(subset(imp, prop_used > 0))$data),
                   nrow(subset(imp, prop_used > 0)))
})
