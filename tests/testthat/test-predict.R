test_that("link and response predictions have the documented shapes", {
  d <- sim_x(seed = 21)
  n <- nrow(d)
  d$y <- stats::rbinom(n, 1, stats::plogis(2 * d$x1 - 1))

  fit <- genbart(y ~ ., data = d, family = binomial(),
                 control = quick_control())

  expect_length(predict(fit), n)
  expect_identical(dim(predict(fit, draws = TRUE)), c(30L, n))

  p <- predict(fit, type = "response")
  expect_length(p, n)
  expect_true(all(p > 0 & p < 1))

  probs <- predict(fit, type = "prob")
  expect_identical(dim(probs), c(n, 2L))
  expect_equal(unname(rowSums(probs)), rep(1, n))

  expect_identical(dim(predict(fit, type = "prob", draws = TRUE)),
                   c(30L, n, 2L))

  cls <- predict(fit, type = "class")
  expect_s3_class(cls, "factor")
  expect_length(cls, n)
})

test_that("multi-predictor families return one column per predictor", {
  d <- sim_x(seed = 22)
  n <- nrow(d)
  d$y <- factor(sample(c("a", "b", "c"), n, replace = TRUE))

  # The default symmetric coding carries one forest per category.
  fit <- genbart(y ~ ., data = d, family = multinomial(),
                 control = quick_control())

  link <- predict(fit, type = "link")
  expect_identical(dim(link), c(n, 3L))
  expect_identical(colnames(link), c("a", "b", "c"))

  expect_type(predict(fit, type = "link", draws = TRUE), "list")
  expect_length(predict(fit, type = "link", draws = TRUE), 3L)

  # For a categorical response the response scale is the category
  # probabilities, so it agrees with type = "prob".
  expect_identical(dim(predict(fit, type = "prob")), c(n, 3L))
  expect_equal(predict(fit), predict(fit, type = "prob"))
  expect_identical(dim(predict(fit, draws = TRUE)), c(30L, n, 3L))
})

test_that("reference coding drops the reference category's predictor", {
  d <- sim_x(seed = 22)
  n <- nrow(d)
  d$y <- factor(sample(c("a", "b", "c"), n, replace = TRUE))

  fit <- genbart(y ~ ., data = d, family = multinomial(reference = "b"),
                 control = quick_control())

  link <- predict(fit, type = "link")
  expect_identical(dim(link), c(n, 2L))
  # The reference is moved to the front, so the remaining levels keep their
  # original order behind it.
  expect_identical(colnames(link), c("a", "c"))

  probs <- predict(fit, type = "prob")
  expect_identical(colnames(probs), c("b", "a", "c"))
  expect_equal(unname(rowSums(probs)), rep(1, n))
})

test_that("predictions on new data reproduce predictions on the fitting data", {
  d <- sim_x(seed = 23)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  for (soft in c(TRUE, FALSE)) {
    fit <- genbart(y ~ ., data = d, control = quick_control(soft = soft))
    expect_predictor_invariant(fit, d)

    # Predicting the same rows in a different order permutes the answer and
    # nothing else.
    perm <- rev(seq_len(nrow(d)))
    expect_equal(predict(fit, newdata = d[perm, ]),
                 predict(fit)[perm], tolerance = 1e-6)
  }
})

test_that("iterations selects a subset of draws", {
  d <- sim_x(seed = 24)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- genbart(y ~ ., data = d, control = quick_control())

  expect_identical(nrow(predict(fit, draws = TRUE, iterations = 1:5)), 5L)
  expect_equal(predict(fit, draws = TRUE, iterations = 3)[1, ],
               fit[["eta"]][["eta"]][3, ])

  expect_error(predict(fit, iterations = 0), "must lie between")
  expect_error(predict(fit, iterations = 999), "must lie between")
})

test_that("new data must supply the predictors the model was fit with", {
  d <- sim_x(seed = 25)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- genbart(y ~ ., data = d, control = quick_control())

  expect_error(predict(fit, newdata = d[, c("x1", "x2")]), "x3")

  nd <- d
  is.na(nd$x1[1L]) <- TRUE
  expect_error(predict(fit, newdata = nd), "missing values")
})

test_that("factor levels unseen in fitting are rejected, known ones are kept", {
  d <- sim_x(seed = 26)
  d$g <- factor(sample(c("a", "b"), nrow(d), replace = TRUE))
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- genbart(y ~ ., data = d, control = quick_control())

  nd <- d[1:5, ]
  nd$g <- factor(c("a", "b", "a", "b", "a"), levels = c("a", "b"))
  expect_length(predict(fit, newdata = nd), 5L)

  bad <- d[1:5, ]
  bad$g <- factor(rep("z", 5))
  expect_error(predict(fit, newdata = bad))
})

test_that("prob and class are refused for families without categories", {
  d <- sim_x(seed = 27)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- genbart(y ~ ., data = d, control = quick_control())

  expect_error(predict(fit, type = "prob"), "available only for")
  expect_error(predict(fit, type = "class"), "available only for")
})

test_that("survival predictions are on the time scale and positive", {
  d <- sim_x(seed = 28)
  n <- nrow(d)
  log_t <- 2 * d$x1 + stats::rnorm(n)
  d$time <- exp(log_t)
  d$status <- 1

  fit <- genbart(survival::Surv(time, status) ~ x1 + x2 + x3, data = d,
                 family = lognormal_aft(), control = quick_control())

  med <- predict(fit, type = "response")
  expect_length(med, n)
  expect_true(all(med > 0))

  # For a log-normal accelerated failure time model the error has median zero,
  # so the median survival time is the exponentiated predictor.
  expect_equal(med, colMeans(exp(fit[["eta"]][["eta"]])), tolerance = 1e-8)
})

test_that("ordinal class predictions come back ordered", {
  d <- sim_x(seed = 29)
  d$y <- factor(sample(1:3, nrow(d), replace = TRUE), levels = 1:3,
                ordered = TRUE)

  fit <- genbart(y ~ ., data = d, family = ordinal(), control = quick_control())

  expect_s3_class(predict(fit, type = "class"), "ordered")
})
