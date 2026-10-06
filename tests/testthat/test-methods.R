test_that("variable_importance() reports usage and separates signal from noise", {
  skip_on_cran()

  d <- sim_x(n = 300, p = 4, seed = 771)
  set.seed(7711)
  d$y <- 2 * d$x1 + sin(3 * d$x2) + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(num_trees = 20L, num_burn = 200L,
                                          num_draws = 200L, sparsity = TRUE))

  vi <- variable_importance(fit)

  # The two predictors in the truth are used in far more draws than the two
  # that are not. The gap is what makes this usable as a selection rule.
  used <- stats::setNames(vi$prop_used, vi$variable)
  expect_gt(min(used[c("x1", "x2")]), max(used[c("x3", "x4")]))
})

test_that("variable_importance() returns a sorted table that summary() agrees with", {
  d <- sim_x(n = 100, p = 4, seed = 771)
  set.seed(7711)
  d$y <- 2 * d$x1 + sin(3 * d$x2) + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(sparsity = TRUE))

  vi <- variable_importance(fit)

  expect_s3_class(vi, "bartisan_importance")
  expect_s3_class(vi, "data.frame")
  expect_identical(names(vi), c("variable", "prop_used", "prop_splits",
                               "splits", "splits_lower", "splits_upper"))
  expect_setequal(vi$variable, c("x1", "x2", "x3", "x4"))

  # Sorted by prop_used then splits, both decreasing.
  expect_false(is.unsorted(rev(vi$prop_used)))

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
  expect_error(variable_importance(1), "must inherit from class")
  expect_error(variable_importance(1, level = 1), "must inherit from class")
})

test_that("prop_splits is a share, and so survives a change of forest size", {
  skip_on_cran()

  d <- sim_x(n = 300, p = 4, seed = 773)
  set.seed(7731)
  d$y <- 2 * d$x1 + sin(3 * d$x2) + stats::rnorm(nrow(d), sd = 0.3)

  # Each fit is seeded rather than started from wherever the one before it left
  # the stream, so neither is a property of how many draws the other took.
  set.seed(7732)
  small <- bartisan(y ~ ., d, family = stats::gaussian(),
                    control = quick_control(num_trees = 10L, num_burn = 200L,
                                            num_draws = 200L, sparsity = TRUE))
  set.seed(7733)
  large <- bartisan(y ~ ., d, family = stats::gaussian(),
                    control = quick_control(num_trees = 40L, num_burn = 200L,
                                            num_draws = 200L, sparsity = TRUE))

  a <- variable_importance(small)
  b <- variable_importance(large)

  # Four times the trees means several times the rules, which is what makes the
  # raw count useless for comparing the two fits.
  expect_gt(sum(b$splits) / sum(a$splits), 2)

  # The share does not scale that way, and the ranking is what carries over.
  # It is not invariant: measured here, the leading predictor's share moved from
  # .51 to .64 when the forest went from 10 trees to 40, so the column removes
  # the dependence on `num_trees` without removing the forest from the answer.
  #
  # Which of the two leads is not the claim, and asserting it was asking for a
  # coin flip: `y` depends on both `x1` and `x2`, and holding this data fixed
  # while the sampler's stream moved, the leading predictor matched across the
  # two forest sizes in 14 of 25 draws. The pair carries over in 25 of 25.
  expect_setequal(a$variable[1:2], b$variable[1:2])
})

test_that("prop_splits is a share of the forest's rules", {
  d <- sim_x(n = 100, p = 4, seed = 773)
  set.seed(7731)
  d$y <- 2 * d$x1 + sin(3 * d$x2) + stats::rnorm(nrow(d), sd = 0.3)

  set.seed(7732)
  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(sparsity = TRUE))

  vi <- variable_importance(fit)

  # The shares add to one within each fit, which is the property that makes
  # them comparable between fits.
  expect_equal(sum(vi$prop_splits), 1)
  expect_true(all(vi$prop_splits >= 0 & vi$prop_splits <= 1))
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
  # unreadable as a selection rule. Both the table and the notes below it go to
  # stdout, so the header row is what says which columns were printed; the note
  # naming the two interval columns is itself one of the lines captured.
  shown <- capture.output(print(variable_importance(on_fit)))
  header <- grep("prop_used", shown, value = TRUE)

  expect_length(header, 1L)
  expect_false(any(grepl("splits_lower", header)))

  expect_match(printed_text(variable_importance(on_fit)), "95% interval")
  expect_no_match(printed_text(variable_importance(on_fit)),
                  "sparsity = FALSE", fixed = TRUE)

  expect_match(printed_text(variable_importance(off_fit)), "sparsity = FALSE",
               fixed = TRUE)

  # And a subset still prints, which it cannot do by carrying the attributes.
  expect_output(print(subset(variable_importance(on_fit), prop_used > 0)),
                "Variable importance")
})

test_that("plot() draws the importance table", {
  skip_if_not_installed("ggplot2")

  d <- sim_x(n = 200, p = 3, seed = 776)
  set.seed(7761)
  d$y <- d$x1 + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(num_trees = 10L))

  expect_s3_class(plot(variable_importance(fit)), "ggplot")

  # One panel per forest when there is more than one.
  ls_fit <- bartisan(y ~ ., d, family = gaussian_ls(),
                     control = quick_control(num_trees = 10L))
  expect_s3_class(plot(variable_importance(ls_fit)), "ggplot")

  # Drawing is the method's job alone; there is no argument that does it.
  expect_error(variable_importance(fit, plot = TRUE), "unused argument")
})

test_that("the plot method draws the table, and a subset of it", {
  skip_if_not_installed("ggplot2")

  d <- sim_x(n = 200, p = 4, seed = 777)
  set.seed(7771)
  d$y <- 2 * d$x1 + sin(3 * d$x2) + stats::rnorm(nrow(d), sd = 0.3)

  fit <- bartisan(y ~ ., d, family = stats::gaussian(),
                  control = quick_control(num_trees = 10L, num_draws = 100L))

  imp <- variable_importance(fit)
  expect_s3_class(plot(imp), "ggplot")

  # Subsetting first is how a wide model is made readable, so the method has to
  # accept what `head()` and `subset()` return.
  expect_identical(nrow(plot(head(imp, 2L))$data), 2L)
  expect_identical(nrow(plot(subset(imp, prop_used > 0))$data),
                   nrow(subset(imp, prop_used > 0)))
})

# The parts of a fit that `print()` and `summary()` describe only when the model
# has them: a random part, forests of different sizes, and a scale for each
# grouping factor in each additive predictor.
test_that("print() and summary() report the group intercepts and their scales", {
  d <- sim_x(n = 100L, p = 2L, seed = 778L)
  d$g <- factor(rep(letters[1:4], length.out = nrow(d)))
  set.seed(7781)
  d$y <- d$x1 + stats::rnorm(nrow(d))

  ctrl <- quick_control(num_burn = 10L, num_draws = 10L)

  fit <- bartisan(y ~ x1 + x2 + (1 | g), d, family = stats::gaussian(),
                  control = ctrl)

  expect_match(printed_text(fit), "Random intercepts: g (4 levels)",
               fixed = TRUE)

  # One row per grouping factor, summarizing the draws of its scale.
  scales <- summary(fit)[["tau"]]
  expect_identical(rownames(scales), "g")
  expect_equal(unname(scales[, "mean"]), mean(fit[["tau"]][[1L]][, "g"]))
  expect_match(printed_text(summary(fit)), "Random-effect scales",
               fixed = TRUE)

  # With two additive predictors each has its own scale, and the rows say which
  # predictor they belong to.
  ls_fit <- bartisan(y ~ x1 + x2 + (1 | g), d, family = gaussian_ls(),
                     control = quick_control(num_trees = c(mean = 5L, log_sd = 3L),
                                             num_burn = 10L, num_draws = 10L))

  expect_identical(rownames(summary(ls_fit)[["tau"]]),
                   c("g [mean]", "g [log_sd]"))
  expect_match(printed_text(ls_fit), "2 forests of 5 and 3 trees", fixed = TRUE)
})

test_that("the family line names the family the caller asked for", {
  d <- sim_x(n = 90L, p = 2L, seed = 779L)
  set.seed(7791)
  d$b <- stats::rbinom(nrow(d), 1L, stats::plogis(d$x1))
  d$m <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  d$y <- stats::rnorm(nrow(d))

  ctrl <- quick_control(num_burn = 10L, num_draws = 10L)

  # A link the engine does not carry is applied on the R side, and says so.
  cauchit <- bartisan(b ~ x1 + x2, d, family = stats::binomial("cauchit"),
                      control = ctrl)
  expect_match(printed_text(cauchit),
               'Family: "binomial" with the "cauchit" link (supplied from R)',
               fixed = TRUE)

  # The engine's "mnp" is the multinomial probit the caller wrote.
  mnp <- bartisan(m ~ x1 + x2, d, family = multinomial("probit"),
                  control = ctrl)
  expect_match(printed_text(mnp),
               'Family: "multinomial" with the "probit" link', fixed = TRUE)

  # A custom family reports the name it was given rather than "custom", which
  # with the identity link would say nothing.
  normal <- custom_family(
    function(y, eta, ...) stats::dnorm(y, eta[, 1L], 1, log = TRUE),
    num_predictors = 1L, name = "unit normal")
  custom <- bartisan(y ~ x1 + x2, d, family = normal, control = ctrl)
  expect_match(printed_text(custom), 'Family: "unit normal" (supplied from R)',
               fixed = TRUE)
})

# `summary()` leaves out a section it has nothing for rather than failing on it,
# so the rest of the summary still prints.
test_that("a summary prints around what the fit does not have", {
  d <- sim_x(n = 60L, p = 2L, seed = 780L)
  set.seed(7801)
  d$y <- d$x1 + stats::rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + x2, d, family = stats::gaussian(),
                  control = quick_control(num_burn = 10L, num_draws = 10L))

  no_loglik <- fit
  no_loglik[["loglik"]][] <- NA_real_
  s <- summary(no_loglik)

  expect_null(s[["convergence"]])
  expect_no_match(printed_text(s), "Log likelihood: R-hat", fixed = TRUE)
  expect_match(printed_text(s), "Convergence and mixing", fixed = TRUE)

  no_counts <- fit
  no_counts[["counts"]] <- NULL
  s <- summary(no_counts)

  expect_null(s[["importance"]])
  expect_match(printed_text(s), "Convergence and mixing", fixed = TRUE)
})

# A proportional hazards baseline has one rate per time bin, which is too many
# rows to read, so the summary prints the first and last six.
test_that("a long block of nuisance parameters is shown by its ends", {
  d <- sim_x(n = 80L, p = 2L, seed = 781L)
  set.seed(7811)
  d$time <- stats::rexp(nrow(d), exp(-d$x1))
  d$status <- stats::rbinom(nrow(d), 1L, 0.8)

  fit <- bartisan(cbind(time, status) ~ x1 + x2, d, family = ph(num_bins = 12L),
                  control = quick_control(num_burn = 10L, num_draws = 10L))

  aux <- summary(fit)[["aux"]]
  shown <- printed_text(summary(fit))
  hidden <- rownames(aux)[-c(1:6, nrow(aux) - 5:0)]

  expect_gt(nrow(aux), 12L)
  expect_match(shown, sprintf("%d more, omitted", length(hidden)), fixed = TRUE)

  for (nm in rownames(aux)) {
    if (nm %in% hidden) {
      expect_no_match(shown, nm, fixed = TRUE)
    }
    else {
      expect_match(shown, nm, fixed = TRUE)
    }
  }
})
