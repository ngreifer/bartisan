# Missing predictor values, handled by missingness incorporated in attributes:
# the splitting rule itself says where a missing value goes, and that choice is
# drawn from its prior with the variable and the cutpoint.

with_missing <- function(n = 120, seed = 51, rate = 0.25, columns = "x1") {
  d <- sim_x(n = n, seed = seed)
  d$y <- 2 * d$x1 + stats::rnorm(n, sd = 0.4)

  for (nm in columns) {
    d[[nm]][stats::runif(n) < rate] <- NA
  }

  d
}

test_that("na.action decides whether missing predictors are kept", {
  d <- with_missing()
  complete <- sum(!is.na(d$x1))

  # na.omit drops the rows, which is what every other modelling function does
  # and is no longer the default here.
  dropped <- genbart(y ~ ., data = d, na.action = stats::na.omit,
                     control = quick_control())
  expect_identical(dropped[["n"]], complete)
  expect_false(any(dropped[["has_na"]]))

  # The default keeps them, and so does asking for na.pass by name.
  by_default <- genbart(y ~ ., data = d, control = quick_control())
  expect_identical(by_default[["n"]], nrow(d))
  expect_identical(names(which(by_default[["has_na"]])), "x1")

  kept <- genbart(y ~ ., data = d, na.action = stats::na.pass,
                  control = quick_control())
  expect_identical(kept[["n"]], nrow(d))
  expect_identical(names(which(kept[["has_na"]])), "x1")

  # Nothing missing anywhere means nothing to record, whichever na.action.
  clean <- sim_x(n = 60, seed = 52)
  clean$y <- stats::rnorm(60)
  fit <- genbart(y ~ ., data = clean, na.action = stats::na.pass,
                 control = quick_control())
  expect_false(any(fit[["has_na"]]))
})

test_that("a fit with missing predictors predicts and stays consistent", {
  d <- with_missing(columns = c("x1", "x2"))

  fit <- genbart(y ~ ., data = d, na.action = stats::na.pass,
                 control = quick_control())

  expect_false(anyNA(predict(fit)))
  expect_false(anyNA(predict(fit, type = "link")))

  # The invariant that matters: replaying the stored forests over the same data,
  # missing values and all, has to reproduce the predictor the sampler recorded.
  # That is what checks the missing-value branch is encoded and decoded right.
  expect_predictor_invariant(fit, d)
})

test_that("hard rules handle missing values too", {
  d <- with_missing(seed = 53)

  fit <- genbart(y ~ ., data = d, na.action = stats::na.pass,
                 control = quick_control(soft = FALSE))

  expect_false(anyNA(predict(fit)))
  expect_predictor_invariant(fit, d)
})

test_that("a missing factor level is a level the rules can split on", {
  d <- sim_x(n = 120, seed = 54)
  d$g <- factor(sample(c("a", "b", "c"), 120, replace = TRUE))
  d$g[stats::runif(120) < 0.2] <- NA
  d$y <- 2 * d$x1 + stats::rnorm(120, sd = 0.4)

  fit <- genbart(y ~ ., data = d, na.action = stats::na.pass,
                 control = quick_control())

  # A missing factor value leaves every one of its indicator columns missing.
  expect_setequal(names(which(fit[["has_na"]])), c("ga", "gb", "gc"))
  expect_predictor_invariant(fit, d)
})

test_that("a column constant where observed is kept for its missingness", {
  d <- sim_x(n = 100, seed = 55)
  d$flat <- 1
  d$flat[stats::runif(100) < 0.3] <- NA
  d$y <- stats::rnorm(100)

  # Constant columns are dropped, but this one carries information in whether it
  # is there at all, so it survives.
  fit <- genbart(y ~ ., data = d, na.action = stats::na.pass,
                 control = quick_control())
  expect_true("flat" %in% names(fit[["has_na"]]))
  expect_true(fit[["has_na"]][["flat"]])

  # With nothing missing it is constant and goes.
  d$flat <- 1
  expect_warning(genbart(y ~ ., data = d, na.action = stats::na.pass,
                         control = quick_control()),
                 "constant predictor column")
})

test_that("rows the model cannot use are dropped and named as such", {
  d <- with_missing(seed = 56)
  d$y[1:5] <- NA

  expect_warning(fit <- genbart(y ~ ., data = d, na.action = stats::na.pass,
                                control = quick_control()),
                 "missing response")
  expect_identical(fit[["n"]], nrow(d) - 5L)

  # A missing weight is the same situation.
  d2 <- with_missing(seed = 57)
  w <- rep(1, nrow(d2))
  w[1:3] <- NA
  expect_warning(fit <- genbart(y ~ ., data = d2, weights = w,
                                na.action = stats::na.pass,
                                control = quick_control()),
                 "missing response")
  expect_identical(fit[["n"]], nrow(d2) - 3L)
})

test_that("new data may only be missing where the fitting data was", {
  d <- with_missing(seed = 58)

  fit <- genbart(y ~ ., data = d, na.action = stats::na.pass,
                 control = quick_control())

  expect_false(anyNA(predict(fit, newdata = d)))

  # x2 had no missing values, so no rule in the forest carries an answer for one.
  d2 <- d
  d2$x2[1] <- NA
  expect_error(predict(fit, newdata = d2), "has missing values in")

  # And a fit with no missing values at all rejects any missing new data.
  clean <- sim_x(n = 60, seed = 59)
  clean$y <- stats::rnorm(60)
  plain <- genbart(y ~ ., data = clean, control = quick_control())
  clean$x1[1] <- NA
  expect_error(predict(plain, newdata = clean), "has missing values in")
})

test_that("informative missingness is recovered", {
  skip_on_cran()

  # The response depends on whether x1 is there, and on nothing else about it.
  # No imputation scheme can find this; a rule that splits on missingness can.
  set.seed(60)
  n <- 400
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  x3 = stats::runif(n))
  gone <- stats::runif(n) < 0.35
  d$y <- 2 * gone + stats::rnorm(n, sd = 0.3)
  d$x1[gone] <- NA

  fit <- genbart(y ~ ., data = d, na.action = stats::na.pass,
                 control = genbart_control(num_trees = 20, num_burn = 300,
                                           num_save = 300, verbose = FALSE))

  p <- predict(fit)
  expect_equal(mean(p[gone]), 2, tolerance = 0.15)
  expect_equal(mean(p[!gone]), 0, tolerance = 0.15)

  # And the sampler should be spending its splits on x1, since it is the only
  # variable that carries anything.
  used <- colMeans(fit[["counts"]][["eta"]])
  expect_gt(used[["x1"]], 5 * max(used[c("x2", "x3")]))
})

test_that("a flat likelihood still reproduces the tree prior with missing data", {
  skip_on_cran()

  # The sharpest check that the missing-value branch of a rule is drawn from its
  # prior and therefore cancels from every acceptance ratio. With the likelihood
  # flat the target is the tree prior, which says nothing about missing values,
  # so the distribution of tree sizes must be the same as it is with complete
  # data. If the extra draw were mishandled the tree size would drift.
  leaf_pmf <- function(gamma = 0.95, beta = 2, max_leaves = 12L,
                       max_depth = 60L) {
    p <- c(1, rep(0, max_leaves - 1L))
    for (d in max_depth:0) {
      rho <- gamma * (1 + d)^(-beta)
      conv <- vapply(2:max_leaves, function(L) {
        sum(p[1:(L - 1L)] * p[(L - 1L):1])
      }, numeric(1))
      p <- c(1 - rho, rho * conv)
    }
    p
  }

  set.seed(61)

  n <- 120
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  x3 = stats::runif(n))
  d$y <- stats::rnorm(n)
  # Missing in every predictor, so every rule the sampler draws needs the extra
  # choice and none of them can dodge it.
  for (nm in c("x1", "x2", "x3")) {
    d[[nm]][stats::runif(n) < 0.3] <- NA
  }

  fit <- genbart(y ~ ., data = d, weights = rep(1e-10, n),
                 na.action = stats::na.pass,
                 control = genbart_control(num_trees = 1, num_burn = 2000,
                                           num_save = 20000, soft = FALSE,
                                           sigma_mu = 0.4,
                                           update_sigma_mu = FALSE,
                                           update_s = FALSE,
                                           update_alpha = FALSE,
                                           sigma_mu_ramp = 0, verbose = FALSE))
  sampled <- rowSums(fit[["counts"]][["eta"]]) + 1L

  reference <- leaf_pmf()
  expect_equal(mean(sampled), sum(reference * seq_along(reference)),
               tolerance = 0.03)
})
