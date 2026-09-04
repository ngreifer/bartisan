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
                  family = gaussian_ls(), control = ls_control())

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

  a <- bartisan(list(y ~ x1, ~ x3), data = d, family = gaussian_ls(),
                control = ls_control())
  b <- bartisan(list(log_sd = ~ x3, mean = y ~ x1), data = d,
                family = gaussian_ls(), control = ls_control())

  expect_identical(unname(colMeans(a[["counts"]][["log_sd"]])[["x1"]]), 0)
  expect_identical(unname(colMeans(b[["counts"]][["log_sd"]])[["x1"]]), 0)
  expect_identical(unname(colMeans(b[["counts"]][["mean"]])[["x3"]]), 0)
})

test_that("only the first formula needs a response, and a second must agree", {
  d <- sim_ls()
  d$z <- d$y

  expect_no_error(bartisan(list(y ~ x1, y ~ x2), data = d,
                           family = gaussian_ls(), control = ls_control()))

  expect_error(bartisan(list(y ~ x1, z ~ x2), data = d,
                        family = gaussian_ls(), control = ls_control()),
               "different response")
  expect_error(bartisan(list(~ x1, ~ x2), data = d,
                        family = gaussian_ls(), control = ls_control()),
               "two-sided")
})

test_that("a per-forest formula expands . without taking in the response", {
  d <- sim_ls()

  # `~ .` on a later formula has no left-hand side of its own, so the response
  # has to be put back before the terms are taken or the outcome becomes a
  # predictor of itself.
  fit <- bartisan(list(y ~ x1, ~ .), data = d, family = gaussian_ls(),
                  control = ls_control())

  expect_false("y" %in% attr(stats::terms(fit), "term.labels"))
  expect_setequal(attr(stats::terms(fit), "term.labels"), c("x1", "x2", "x3"))
})

test_that("an interaction written either way is one predictor", {
  d <- sim_ls()

  fit <- bartisan(list(y ~ x1 * x2, ~ x2:x1), data = d,
                  family = gaussian_ls(), control = ls_control())

  expect_setequal(attr(stats::terms(fit), "term.labels"),
                  c("x1", "x2", "x1:x2"))
})

test_that("prior settings are per forest, named or positional", {
  d <- sim_ls()

  named <- bartisan(y ~ ., data = d, family = gaussian_ls(),
                    control = ls_control(num_trees = c(mean = 8, log_sd = 3),
                                         sparsity = c(mean = TRUE,
                                                      log_sd = FALSE)))
  positional <- bartisan(y ~ ., data = d, family = gaussian_ls(),
                         control = ls_control(num_trees = c(8, 3),
                                              sparsity = c(TRUE, FALSE)))

  expect_identical(named[["num_trees"]], c(8L, 3L))
  expect_identical(positional[["num_trees"]], c(8L, 3L))
  expect_identical(unname(named[["control"]][["update_s"]]), c(TRUE, FALSE))

  # A forest a named argument does not mention keeps the default rather than
  # borrowing the value given for the other forest.
  partial <- bartisan(y ~ ., data = d, family = gaussian_ls(),
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
                  family = gaussian_ls(),
                  control = ls_control(num_trees = c(8, log_sd = 3)))

  expect_identical(fit[["num_trees"]], c(8L, 3L))
  expect_identical(unname(colMeans(fit[["counts"]][["mean"]])[["x2"]]), 0)
  expect_identical(unname(colMeans(fit[["counts"]][["log_sd"]])[["x1"]]), 0)

  # More unnamed entries than forests left to take them is an error.
  expect_error(bartisan(y ~ ., data = d, family = gaussian_ls(),
                        control = ls_control(num_trees = c(log_sd = 3, 8, 9))),
               "unnamed entr")
})

test_that("per-forest arguments are checked against the family's forests", {
  d <- sim_ls()

  expect_error(bartisan(y ~ ., data = d, family = gaussian_ls(),
                        control = ls_control(num_trees = c(mu = 8))),
               "does not have")
  expect_error(bartisan(y ~ ., data = d, family = gaussian_ls(),
                        control = ls_control(num_trees = c(8, 3, 2))),
               "3 values")
  expect_error(bartisan(list(y ~ x1, ~ x2, ~ x3), data = d,
                        family = gaussian_ls(), control = ls_control()),
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

test_that("a forest whose predictors are all zero-weighted is an error", {
  d <- sim_ls()

  # The forest's formula names `x3`, so this is a `split_prior` that has zeroed
  # out the only thing it could have split on, which is a mistake rather than a
  # statement. An intercept-only formula is the way to say the parameter is
  # constant; see the test below.
  expect_error(bartisan(list(y ~ x1 + x2, ~ x3), data = d,
                        family = gaussian_ls(),
                        control = ls_control(split_prior = c(x3 = 0))),
               "weight of zero")
})

test_that("an intercept-only forest is a constant parameter, not an error", {
  d <- sim_ls()

  set.seed(3)
  fit <- bartisan(list(y ~ x1 + x2 + x3, ~ 1), data = d,
                  family = gaussian_ls(), control = ls_control())

  # Every tree in that forest is a stump, so it never splits and its predictor
  # is one drawn scalar rather than a function of the covariates.
  expect_equal(sum(fit[["counts"]][["log_sd"]]), 0)
  expect_equal(stats::sd(colMeans(fit[["eta"]][[2L]])), 0)

  # The mean forest is unaffected and still splits.
  expect_gt(sum(fit[["counts"]][["mean"]]), 0)
})

test_that("every family taking several formulas accepts an intercept-only one", {
  # The point of `~ 1` is that it works wherever more than one formula does, so
  # that a nuisance parameter and an empty forest are close to the same thing.
  # The multinomial families are excluded because their forests are the levels
  # of one vector-valued parameter rather than separate components, which is
  # what `joint_forests()` marks and what refuses a formula list outright.
  set.seed(12)
  n <- 300L
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$num <- 2 * d$x1 + stats::rnorm(n, sd = 0.5)
  d$pos <- stats::rgamma(n, shape = 4, rate = 4 / exp(1 + d$x1))
  d$cnt <- stats::rpois(n, exp(0.5 + d$x1)) * stats::rbinom(n, 1L, 0.7)

  cases <- list(
    list(y = "num", family = gaussian_ls()),
    list(y = "pos", family = Gamma_ls()),
    list(y = "cnt", family = zi_poisson()),
    list(y = "cnt", family = zi_negbin()),
    list(y = "num", family = custom_family(
      function(y, eta) stats::dnorm(y, eta[, 1L], exp(eta[, 2L]), log = TRUE),
      num_predictors = 2L))
  )

  for (case in cases) {
    f <- stats::reformulate(c("x1", "x2"), response = case[["y"]])
    label <- case[["family"]][["family"]]

    set.seed(2)
    fit <- bartisan(list(f, ~ 1), data = d, family = case[["family"]],
                    control = quick_control())

    # The second forest never splits, and its predictor is one number.
    expect_equal(sum(fit[["counts"]][[2L]]), 0, info = label)
    expect_equal(stats::sd(colMeans(fit[["eta"]][[2L]])), 0, info = label)

    # The first forest is untouched.
    expect_gt(sum(fit[["counts"]][[1L]]), 0)
  }
})

test_that("an intercept-only inflation forest is a constant zero probability", {
  # `zi_poisson()` with `~ 1` on its second forest is the ordinary zero-inflated
  # Poisson, with one structural-zero probability rather than a forest of them,
  # so a known constant probability should come back.
  skip_on_cran()

  set.seed(21)
  n <- 1200L
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  lambda <- exp(0.4 + 1.5 * d$x1)
  d$count <- ifelse(stats::rbinom(n, 1L, 0.3) == 1L, 0L,
                    stats::rpois(n, lambda))

  set.seed(1)
  fit <- bartisan(list(count ~ x1 + x2, ~ 1), data = d,
                  family = zi_poisson(),
                  control = bartisan_control(num_trees = 50L, num_burn = 300L,
                                             num_draws = 300L, gate = "hard"))

  expect_equal(mean(stats::plogis(fit[["eta"]][[2L]])), 0.3, tolerance = 0.08)
  expect_lt(sqrt(mean((colMeans(fit[["eta"]][[1L]]) - log(lambda))^2)), 0.2)
})
