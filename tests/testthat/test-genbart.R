test_that("every family fits and returns coherent dimensions", {
  d <- sim_x()
  n <- nrow(d)
  eta <- 2 * d$x1 - 1

  cases <- list(
    gaussian = list(y = eta + stats::rnorm(n), family = gaussian(), h = 1L),
    binomial = list(y = stats::rbinom(n, 1, stats::plogis(eta)),
                    family = binomial(), h = 1L),
    probit = list(y = stats::rbinom(n, 1, stats::pnorm(eta)),
                  family = binomial("probit"), h = 1L),
    cloglog = list(y = stats::rbinom(n, 1, 1 - exp(-exp(eta))),
                   family = binomial("cloglog"), h = 1L),
    poisson = list(y = stats::rpois(n, exp(eta)), family = poisson(), h = 1L),
    negbin = list(y = stats::rnbinom(n, mu = exp(eta), size = 2),
                  family = negbin(), h = 1L),
    Gamma = list(y = stats::rgamma(n, 2, rate = 2 / exp(eta)),
                 family = Gamma("log"), h = 1L),
    location_scale = list(y = eta + stats::rnorm(n), family = location_scale(),
                          h = 2L)
  )

  for (nm in names(cases)) {
    case <- cases[[nm]]
    dd <- d
    dd$y <- case$y

    fit <- genbart(y ~ ., data = dd, family = case$family,
                   control = quick_control())

    expect_s3_class(fit, "genbart")
    expect_identical(fit[["num_forest"]], case$h)
    expect_length(fit[["eta"]], case$h)
    expect_identical(dim(fit[["eta"]][[1L]]), c(30L, n))
    expect_identical(dim(fit[["counts"]][[1L]]), c(30L, 3L))
    expect_true(all(is.finite(fit[["loglik"]])))
    expect_predictor_invariant(fit, dd)
  }
})

test_that("ordinal and multinomial fit and carry their levels", {
  d <- sim_x(seed = 2)
  n <- nrow(d)

  d$y <- factor(sample(1:4, n, replace = TRUE), levels = 1:4, ordered = TRUE)
  fit <- genbart(y ~ ., data = d, family = ordinal(), control = quick_control())

  expect_identical(fit[["num_cat"]], 4L)
  expect_identical(fit[["num_forest"]], 1L)
  # One cutpoint per boundary. With three or more categories all of them are
  # free and the predictor is centered instead; see ?genbart-families.
  expect_identical(colnames(fit[["aux"]]), paste0("cut", 1:3))
  expect_true(all(abs(rowMeans(fit[["eta"]][[1L]])) < 1e-8))
  expect_gt(stats::sd(fit[["aux"]][, "cut1"]), 0)
  # Cutpoints stay ordered, which the sampler enforces by construction.
  expect_true(all(fit[["aux"]][, "cut2"] < fit[["aux"]][, "cut3"]))
  expect_predictor_invariant(fit, d)

  d$y <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  fit <- genbart(y ~ ., data = d, family = multinomial(),
                 control = quick_control())

  expect_identical(fit[["num_forest"]], 3L)
  expect_identical(names(fit[["eta"]]), c("a", "b", "c"))
  expect_predictor_invariant(fit, d)

  fit <- genbart(y ~ ., data = d, family = multinomial(reference = "a"),
                 control = quick_control())

  expect_identical(fit[["num_forest"]], 2L)
  expect_identical(names(fit[["eta"]]), c("b", "c"))
  expect_predictor_invariant(fit, d)
})

test_that("a two-category ordinal model reduces to binary regression", {
  d <- sim_x(seed = 3)
  d$y <- factor(stats::rbinom(nrow(d), 1, 0.5), ordered = TRUE)

  fit <- genbart(y ~ ., data = d, family = ordinal(), control = quick_control())

  # With one boundary there is nothing to draw, so the single cutpoint stays
  # pinned and the model is a logistic regression.
  expect_identical(unname(fit[["aux"]][, "cut1"]), rep(0, 30))
})

test_that("accelerated failure time models fit right-censored data", {
  d <- sim_x(seed = 4)
  n <- nrow(d)
  log_t <- 2 * d$x1 + stats::rlogis(n)
  cens <- stats::quantile(log_t, 0.8)

  d$time <- exp(pmin(log_t, cens))
  d$status <- as.numeric(log_t <= cens)

  for (fam in list(weibull_aft(), loglogistic_aft(), lognormal_aft())) {
    fit <- genbart(survival::Surv(time, status) ~ x1 + x2 + x3, data = d,
                   family = fam, control = quick_control())

    expect_identical(colnames(fit[["aux"]]), "sigma")
    expect_true(all(fit[["aux"]][, "sigma"] > 0))
    expect_predictor_invariant(fit, d)
  }
})

test_that("hard and soft decision rules both fit", {
  d <- sim_x(seed = 5)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  hard <- genbart(y ~ ., data = d, control = quick_control(soft = FALSE))
  soft <- genbart(y ~ ., data = d, control = quick_control(soft = TRUE))

  expect_false(hard[["soft"]])
  expect_true(soft[["soft"]])

  # Bandwidths are irrelevant to hard rules but still recorded, and they are
  # drawn only for soft rules.
  expect_true(stats::sd(as.vector(soft[["bandwidth"]])) > 0)

  expect_predictor_invariant(hard, d)
  expect_predictor_invariant(soft, d)
})

test_that("factors become one group in the sparsity prior", {
  d <- sim_x(seed = 6)
  d$g <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- genbart(y ~ ., data = d, control = quick_control())

  # Four terms, so four groups, even though the factor spans three columns.
  expect_identical(fit[["group_names"]], c("x1", "x2", "x3", "g"))
  expect_identical(ncol(fit[["counts"]][[1L]]), 4L)
  expect_length(fit[["unit_maps"]], 6L)
  expect_predictor_invariant(fit, d)
})

test_that("weights, offsets and subsets are honored", {
  d <- sim_x(seed = 7)
  n <- nrow(d)
  d$y <- stats::rpois(n, exp(1 + d$x1))
  d$logn <- stats::runif(n, 0.5, 2)

  fit <- genbart(y ~ x1 + x2 + x3 + offset(log(logn)), data = d,
                 family = poisson(), control = quick_control())
  expect_s3_class(fit, "genbart")

  # An observation-level offset is not recoverable from the predictors, so
  # predicting new data must insist on being given one rather than quietly
  # returning predictions that omit it.
  expect_error(predict(fit, newdata = d), "requires .*offset")
  expect_length(predict(fit, newdata = d, offset = log(d$logn)), n)

  # Predictions on the fitting data still work, since the offset is already
  # baked into the stored predictor.
  expect_length(predict(fit), n)

  sub <- genbart(y ~ x1 + x2 + x3, data = d, subset = seq_len(40),
                 family = poisson(), control = quick_control())
  expect_identical(sub[["n"]], 40L)

  wt <- genbart(y ~ x1 + x2 + x3, data = d, weights = rep(2, n),
                family = poisson(), control = quick_control())
  expect_s3_class(wt, "genbart")
})

test_that("degenerate inputs are rejected", {
  d <- sim_x(seed = 8)
  d$y <- stats::rnorm(nrow(d))

  expect_error(genbart(y ~ 1, data = d, control = quick_control()),
               "at least one predictor")

  d$constant <- 1
  expect_warning(genbart(y ~ x1 + constant, data = d,
                         control = quick_control()),
                 "constant predictor")

  d2 <- d
  d2$y <- rep(1, nrow(d2))
  expect_error(genbart(y ~ x1, data = d2, control = quick_control()),
               "no variation")

  d3 <- d
  d3$y <- -1
  # Matched on a fragment that cannot span a line break in the message source:
  # testthat stops cli reflowing condition messages, so a regex crossing one
  # fails under R CMD check and passes everywhere else.
  expect_error(genbart(y ~ x1, data = d3, family = poisson(),
                       control = quick_control()),
               "requires a response of non-negative")
  expect_error(genbart(y ~ x1, data = d3, family = Gamma("log"),
                       control = quick_control()),
               "strictly positive")
})

test_that("print and summary report the fit", {
  d <- sim_x(seed = 9)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- genbart(y ~ ., data = d, control = quick_control())

  expect_output(print(fit), "Generalized BART")
  expect_output(print(fit), "soft decision rules")

  s <- summary(fit)
  expect_s3_class(s, "summary.genbart")
  expect_identical(rownames(s[["usage"]][["eta"]]), c("x1", "x2", "x3")[
    order(colMeans(fit[["counts"]][["eta"]] > 0), decreasing = TRUE)])
  expect_true(all(s[["usage"]][["eta"]][, "prop_used"] <= 1))
  expect_output(print(s), "Predictor usage")

  expect_error(summary(fit, level = 1.5), "must be between")
})
