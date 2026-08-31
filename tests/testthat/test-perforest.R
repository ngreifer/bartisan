# Arguments that refer to different forests of one model: a formula per forest,
# and every prior setting keyed by forest name or given positionally.

sim_ls <- function(n = 300, seed = 1) {
  set.seed(seed)
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                  x3 = stats::rnorm(n))
  d$y <- d$x1 + exp(0.5 * d$x2) * stats::rnorm(n)
  d
}

ls_control <- function(...) {
  quick_control(num_trees = 8, num_burn = 60, num_draws = 60, ...)
}

test_that("a formula per forest holds each forest to its own predictors", {
  d <- sim_ls()

  fit <- bartisan(list(y ~ x1 + x2, ~ x2 + x3), data = d,
                  family = location_scale(), control = ls_control())

  # The frame is the union, so both forests have every predictor in their data.
  expect_setequal(attr(stats::terms(fit), "term.labels"),
                  c("x1", "x2", "x3"))

  # And each splits only on the terms its own formula named.
  mean_splits <- colMeans(fit[["counts"]][["mean"]])
  sd_splits <- colMeans(fit[["counts"]][["log_sd"]])

  expect_identical(unname(mean_splits[["x3"]]), 0)
  expect_identical(unname(sd_splits[["x1"]]), 0)
  expect_gt(mean_splits[["x1"]], 0)
  expect_gt(sd_splits[["x3"]], 0)
})

test_that("the list of formulas can be named, and names reorder it", {
  d <- sim_ls()

  a <- bartisan(list(y ~ x1, ~ x3), data = d, family = location_scale(),
                control = ls_control())
  b <- bartisan(list(log_sd = ~ x3, mean = y ~ x1), data = d,
                family = location_scale(), control = ls_control())

  expect_identical(unname(colMeans(a[["counts"]][["log_sd"]])[["x1"]]), 0)
  expect_identical(unname(colMeans(b[["counts"]][["log_sd"]])[["x1"]]), 0)
  expect_identical(unname(colMeans(b[["counts"]][["mean"]])[["x3"]]), 0)
})

test_that("only the first formula needs a response, and a second must agree", {
  d <- sim_ls()
  d$z <- d$y

  expect_no_error(bartisan(list(y ~ x1, y ~ x2), data = d,
                           family = location_scale(), control = ls_control()))

  expect_error(bartisan(list(y ~ x1, z ~ x2), data = d,
                        family = location_scale(), control = ls_control()),
               "different response")
  expect_error(bartisan(list(~ x1, ~ x2), data = d,
                        family = location_scale(), control = ls_control()),
               "two-sided")
})

test_that("a per-forest formula expands . without taking in the response", {
  d <- sim_ls()

  # `~ .` on a later formula has no left-hand side of its own, so the response
  # has to be put back before the terms are taken or the outcome becomes a
  # predictor of itself.
  fit <- bartisan(list(y ~ x1, ~ .), data = d, family = location_scale(),
                  control = ls_control())

  expect_false("y" %in% attr(stats::terms(fit), "term.labels"))
  expect_setequal(attr(stats::terms(fit), "term.labels"), c("x1", "x2", "x3"))
})

test_that("an interaction written either way is one predictor", {
  d <- sim_ls()

  fit <- bartisan(list(y ~ x1 * x2, ~ x2:x1), data = d,
                  family = location_scale(), control = ls_control())

  expect_setequal(attr(stats::terms(fit), "term.labels"),
                  c("x1", "x2", "x1:x2"))
})

test_that("prior settings are per forest, named or positional", {
  d <- sim_ls()

  named <- bartisan(y ~ ., data = d, family = location_scale(),
                    control = ls_control(num_trees = c(mean = 8, log_sd = 3),
                                         sparsity = c(mean = TRUE,
                                                      log_sd = FALSE)))
  positional <- bartisan(y ~ ., data = d, family = location_scale(),
                         control = ls_control(num_trees = c(8, 3),
                                              sparsity = c(TRUE, FALSE)))

  expect_identical(named[["num_trees"]], c(8L, 3L))
  expect_identical(positional[["num_trees"]], c(8L, 3L))
  expect_identical(unname(named[["control"]][["update_s"]]), c(TRUE, FALSE))

  # A forest a named argument does not mention keeps the default rather than
  # borrowing the value given for the other forest.
  partial <- bartisan(y ~ ., data = d, family = location_scale(),
                      control = ls_control(k = c(log_sd = 8)))

  expect_identical(partial[["control"]][["k"]], c(log_sd = 8))
  expect_gt(mean(partial[["sigma_mu"]][, 1L]),
            mean(partial[["sigma_mu"]][, 2L]))
})

test_that("naming some entries and not others matches R's own rule", {
  d <- sim_ls()

  # Named entries take the forest they name; the rest fill what is left in
  # order. `list(y ~ x1, log_sd = ~ x2)` is the case worth having.
  fit <- bartisan(list(y ~ x1, log_sd = ~ x2), data = d,
                  family = location_scale(),
                  control = ls_control(num_trees = c(8, log_sd = 3)))

  expect_identical(fit[["num_trees"]], c(8L, 3L))
  expect_identical(unname(colMeans(fit[["counts"]][["mean"]])[["x2"]]), 0)
  expect_identical(unname(colMeans(fit[["counts"]][["log_sd"]])[["x1"]]), 0)

  # More unnamed entries than forests left to take them is an error.
  expect_error(bartisan(y ~ ., data = d, family = location_scale(),
                        control = ls_control(num_trees = c(log_sd = 3, 8, 9))),
               "unnamed entr")
})

test_that("per-forest arguments are checked against the family's forests", {
  d <- sim_ls()

  expect_error(bartisan(y ~ ., data = d, family = location_scale(),
                        control = ls_control(num_trees = c(mu = 8))),
               "does not have")
  expect_error(bartisan(y ~ ., data = d, family = location_scale(),
                        control = ls_control(num_trees = c(8, 3, 2))),
               "3 values")
  expect_error(bartisan(list(y ~ x1, ~ x2, ~ x3), data = d,
                        family = location_scale(), control = ls_control()),
               "3 formulas")
})

test_that("the multinomial families take one value for all their forests", {
  set.seed(4)
  n <- 200
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  d$y <- factor(sample(c("a", "b", "c"), n, replace = TRUE))

  # Their forests are the levels of one parameter, so splitting a setting across
  # them says nothing a caller would mean.
  expect_error(bartisan(y ~ ., data = d, family = multinomial(),
                        control = ls_control(num_trees = c(8, 3, 2))),
               "one value for this family")
  expect_error(bartisan(list(y ~ x1, ~ x2, ~ x1), data = d,
                        family = multinomial(), control = ls_control()),
               "single formula for this family")

  expect_no_error(bartisan(y ~ ., data = d, family = multinomial(),
                           control = ls_control(num_trees = 8)))
})

test_that("a single formula and scalar settings reach the engine unchanged", {
  d <- sim_ls()

  # The whole point of the recycling rule: the ordinary call is the same call it
  # was, down to the draws.
  set.seed(3)
  a <- bartisan(y ~ x1 + x2 + x3, data = d, family = gaussian(),
                control = ls_control())
  set.seed(3)
  b <- bartisan(list(y ~ x1 + x2 + x3), data = d, family = gaussian(),
                control = ls_control())

  expect_equal(a[["eta"]][[1L]], b[["eta"]][[1L]])
})

test_that("a forest left with nothing to split on is an error", {
  d <- sim_ls()

  expect_error(bartisan(list(y ~ x1 + x2, ~ x3), data = d,
                        family = location_scale(),
                        control = ls_control(split_prior = c(x3 = 0))),
               "nothing in common")
})
