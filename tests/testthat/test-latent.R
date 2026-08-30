# `type = "mean"` and `type = "stdlv"`, both of which follow WeightIt's ordinal
# predictions. Each is a transformation of quantities the model already reports,
# so what the tests check is that the transformation is the stated one -- against
# WeightIt where it can be, and against arithmetic where it cannot.

test_that("the mean prediction is the probabilities weighted by the labels", {
  d <- sim_x(n = 200, seed = 121)
  set.seed(1121)
  z <- 2 * d$x1 - d$x2 + stats::rnorm(nrow(d))
  d$y <- ordered(rowSums(outer(z, stats::quantile(z, c(1, 2) / 3), ">")) + 1L,
                 levels = 1:3, labels = c("1", "2", "4"))

  fit <- bartisan(y ~ ., d, family = ordinal("probit"),
                  control = quick_control(gate = "hard"))

  p <- stats::predict(fit, newdata = d, type = "prob")
  expected <- drop(p %*% c(1, 2, 4))

  expect_equal(stats::predict(fit, newdata = d, type = "mean"), expected,
               ignore_attr = TRUE)

  # Bounded by the extreme labels, since it is a weighted average of them.
  m <- stats::predict(fit, newdata = d, type = "mean")
  expect_true(all(m >= 1 & m <= 4))

  # `values` overrides the labels and is used as given.
  v <- c("1" = 0, "2" = 10, "4" = 100)
  expect_equal(stats::predict(fit, newdata = d, type = "mean", values = v),
               drop(p %*% c(0, 10, 100)), ignore_attr = TRUE)

  # Draws come back as draws by observations, and their mean is the point
  # estimate -- so the two orientations agree.
  dm <- stats::predict(fit, newdata = d, type = "mean", draws = TRUE)
  expect_identical(dim(dm), c(30L, 200L))
  expect_equal(colMeans(dm), stats::predict(fit, newdata = d, type = "mean"),
               ignore_attr = TRUE)
})

test_that("the mean prediction refuses labels it cannot read, and says what to do", {
  d <- sim_x(n = 120, seed = 122)
  set.seed(1122)
  z <- d$x1 - d$x2 + stats::rnorm(nrow(d))
  d$y <- ordered(rowSums(outer(z, stats::quantile(z, c(1, 2) / 3), ">")) + 1L,
                 levels = 1:3, labels = c("low", "mid", "high"))

  fit <- bartisan(y ~ ., d, family = ordinal("probit"), control = quick_control())

  expect_error(stats::predict(fit, type = "mean"), "cannot be read")

  # And `values` is the way through.
  out <- stats::predict(fit, type = "mean",
                        values = c(low = 1, mid = 2, high = 3))
  expect_true(all(out >= 1 & out <= 3))

  # A `values` that does not name every level is refused rather than recycled.
  expect_error(stats::predict(fit, type = "mean", values = c(low = 1, mid = 2)),
               "named for every response level")
  expect_error(stats::predict(fit, type = "mean", values = c(1, 2, 3)),
               "named for every response level")

  # And it is ignored, with a warning, for any other type.
  expect_warning(stats::predict(fit, type = "prob",
                                values = c(low = 1, mid = 2, high = 3)),
                 "ignored")
})

test_that("a binomial mean prediction with 0/1 labels is the fitted probability", {
  d <- sim_x(n = 150, seed = 123)
  set.seed(1123)
  d$y <- factor(stats::rbinom(nrow(d), 1, stats::plogis(2 * d$x1 - 1)),
                levels = 0:1)

  fit <- bartisan(y ~ ., d, family = binomial(), control = quick_control())

  # E[Y] is P(Y = 1) when the labels are 0 and 1, so the two types must agree
  # exactly rather than approximately.
  expect_equal(stats::predict(fit, type = "mean"),
               stats::predict(fit, type = "response"), ignore_attr = TRUE)
})

test_that("the standardized latent variable has the scale it claims", {
  skip_on_cran()

  set.seed(124)
  n <- 1500
  d <- data.frame(x1 = stats::runif(n, -1, 2), x2 = stats::runif(n, 0, 3))
  lp <- 1.4 * d$x1 - 0.7 * d$x2
  cuts <- c(-1, 0.5, 2)

  errors <- list(probit = function(k) stats::rnorm(k),
                 logit = function(k) stats::rlogis(k),
                 cloglog = function(k) log(-log1p(-stats::runif(k))))
  variances <- c(probit = 1, logit = pi^2 / 3, cloglog = pi^2 / 6)
  means <- c(probit = 0, logit = 0, cloglog = digamma(1))

  for (lk in names(errors)) {
    d$y <- ordered(rowSums(outer(lp + errors[[lk]](n), cuts, ">")) + 1L)

    fit <- bartisan(y ~ x1 + x2, d, family = ordinal(lk), gate = "hard",
                    num_trees = 30, num_burn = 300, num_save = 300)

    eta <- stats::predict(fit, newdata = d, type = "link", draws = TRUE)
    got <- stats::predict(fit, newdata = d, type = "stdlv", draws = TRUE)

    # The definition, rebuilt from the predictor and the link's own constants.
    spread <- apply(fit[["eta"]][[1L]], 1L, stats::var)
    scale <- sqrt(spread + variances[[lk]])
    expected <- (eta - means[[lk]]) / scale

    expect_equal(got, expected, info = lk)

    # The divisor is a property of the fit, not of what is being predicted, so
    # predicting a subset does not change the scale.
    half <- stats::predict(fit, newdata = d[1:100, ], type = "stdlv",
                           draws = TRUE)
    expect_equal(half, got[, 1:100], info = lk)
  }
})

test_that("the standardized latent variable matches WeightIt on scale", {
  skip_on_cran()
  skip_if_not_installed("WeightIt")

  # A linear truth, so WeightIt's ordinal model is correctly specified and the
  # two should agree about the scale. They differ by a constant, because they
  # identify the latent's location differently -- WeightIt drops the intercept
  # column, bartisan centers the predictor -- and a standardized quantity is used
  # for differences, which that constant leaves alone.
  set.seed(125)
  n <- 3000
  d <- data.frame(x1 = stats::runif(n, -1, 2), x2 = stats::runif(n, 0, 3))
  lp <- 1.4 * d$x1 - 0.7 * d$x2

  for (lk in c("probit", "logit", "cloglog")) {
    err <- switch(lk,
                  probit = stats::rnorm(n),
                  logit = stats::rlogis(n),
                  cloglog = log(-log1p(-stats::runif(n))))
    d$y <- ordered(rowSums(outer(lp + err, c(-1, 0.5, 2), ">")) + 1L)

    fit <- bartisan(y ~ x1 + x2, d, family = ordinal(lk), gate = "hard",
                    num_trees = 50, num_burn = 300, num_save = 300)
    wi <- WeightIt::ordinal_weightit(y ~ x1 + x2, data = d, link = lk)

    ours <- stats::predict(fit, type = "stdlv")
    theirs <- stats::predict(wi, type = "stdlv")

    # Same scale, and the difference really is a constant rather than noise.
    expect_equal(stats::sd(ours) / stats::sd(theirs), 1, tolerance = 0.05,
                 info = lk)
    expect_gt(stats::cor(ours, theirs), 0.98)
    expect_lt(stats::sd(ours - theirs), 0.1)
  }
})

test_that("the standardized latent variable is refused where there is no latent", {
  d <- sim_x(n = 100, seed = 126)
  set.seed(1126)
  d$y <- stats::rnorm(nrow(d))

  gauss <- bartisan(y ~ ., d, control = quick_control())
  expect_error(stats::predict(gauss, type = "stdlv"), "threshold crossing")

  d$ym <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
  multi <- bartisan(ym ~ x1 + x2, d, family = multinomial(),
                    control = quick_control())
  expect_error(stats::predict(multi, type = "stdlv"), "threshold crossing")

  # A binomial fit has a latent variable, but only under a link that names its
  # distribution -- a Cauchy error has no variance to divide by.
  d$yb <- stats::rbinom(nrow(d), 1L, 0.5)
  cauchit <- bartisan(yb ~ x1 + x2, d, family = binomial("cauchit"),
                      control = quick_control())
  expect_error(stats::predict(cauchit, type = "stdlv"), "known latent")

  # And "mean" is refused for a family with no categories at all.
  expect_error(stats::predict(gauss, type = "mean"), "available only")
})

test_that("a binomial fit has a standardized latent variable too", {
  d <- sim_x(n = 400, seed = 131)
  set.seed(1131)
  signal <- 1.5 * d$x1 - d$x2

  # The scale is the checkable part: the standard deviation of the standardized
  # latent has to be the predictor's own, divided by the standard deviation of
  # the latent it indexes.
  error_variance <- c(probit = 1, logit = pi^2 / 3, cloglog = pi^2 / 6)

  for (link in names(error_variance)) {
    d$y <- stats::rbinom(nrow(d), 1L,
                         stats::binomial(link)$linkinv(signal))
    fit <- bartisan(y ~ x1 + x2, d,
                    control = quick_control(num_trees = 20L, num_burn = 150L,
                                            num_save = 150L),
                   family = stats::binomial(link))

    standardized <- stats::predict(fit, type = "stdlv")
    divisor <- sqrt(mean(apply(fit[["eta"]][[1L]], 1L, stats::var)) +
                      error_variance[[link]])

    expect_equal(stats::sd(standardized),
                 stats::sd(stats::predict(fit, type = "link")) / divisor,
                 tolerance = 0.01)
  }
})

test_that("a symmetric link puts the same error in both families", {
  # The scale and location of the standardized latent come from one function, so
  # this checks it directly rather than through two chains that would only agree
  # up to Monte Carlo error.
  for (link in c("probit", "logit")) {
    expect_identical(latent_error("binomial", link),
                     latent_error("ordinal", link))
  }

  # The complementary log-log error is not symmetric, so the two families add it
  # with opposite signs. Same variance, mirrored mean.
  binary <- latent_error("binomial", "cloglog")
  ordered <- latent_error("ordinal", "cloglog")

  expect_identical(binary[2L], ordered[2L])
  expect_identical(binary[1L], -ordered[1L])
  expect_equal(ordered[1L], digamma(1))
})

test_that("a binomial and a two-category ordinal probit fit agree", {
  d <- sim_x(n = 300, seed = 133)
  set.seed(1133)
  d$y <- stats::rbinom(nrow(d), 1L, stats::pnorm(1.5 * d$x1 - d$x2))
  d$ordered <- factor(d$y, levels = c(0L, 1L), ordered = TRUE)

  # Two ways of writing the same model. They are fitted by different samplers,
  # so they agree up to Monte Carlo error rather than exactly.
  chain <- quick_control(num_trees = 20L, num_burn = 150L, num_save = 150L)
  binary <- bartisan(y ~ x1 + x2, d, family = stats::binomial("probit"),
                     control = chain)
  ordinal_fit <- bartisan(ordered ~ x1 + x2, d, family = ordinal("probit"),
                          control = chain)

  from_binary <- stats::predict(binary, type = "stdlv")
  from_ordinal <- stats::predict(ordinal_fit, type = "stdlv")

  expect_gt(stats::cor(from_binary, from_ordinal), 0.98)
  expect_equal(stats::sd(from_binary), stats::sd(from_ordinal),
               tolerance = 0.1)
  expect_lt(abs(mean(from_binary) - mean(from_ordinal)), 0.05)
})

test_that("the complementary log-log error enters the two families with opposite signs", {
  d <- sim_x(n = 400, seed = 135)
  set.seed(1135)
  signal <- 1.5 * d$x1 - d$x2
  d$y <- stats::rbinom(nrow(d), 1L,
                       stats::binomial("cloglog")$linkinv(signal))

  fit <- bartisan(y ~ x1 + x2, d, family = stats::binomial("cloglog"),
                  control = quick_control(num_trees = 20L, num_burn = 150L,
                                          num_save = 150L))

  # A binomial cloglog model says Y = 1 exactly when a smallest extreme value
  # variate falls below the index, so the latent is `eta - e` and the error it
  # adds has mean `+gamma`. Derived here rather than taken from the code.
  euler <- -digamma(1)
  divisor <- sqrt(mean(apply(fit[["eta"]][[1L]], 1L, stats::var)) + pi^2 / 6)
  predictor <- stats::predict(fit, type = "link")

  expect_equal(mean(stats::predict(fit, type = "stdlv")),
               (mean(predictor) - euler) / divisor,
               tolerance = 0.01)

  # The ordinal convention is the other sign, and it is a long way off, so this
  # is not a test that would pass either way.
  expect_gt(abs(mean(stats::predict(fit, type = "stdlv")) -
                  (mean(predictor) + euler) / divisor),
            0.5)
})
