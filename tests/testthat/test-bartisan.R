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
                 family = stats::Gamma("log"), h = 1L),
    location_scale = list(y = eta + stats::rnorm(n), family = location_scale(),
                          h = 2L)
  )

  for (nm in names(cases)) {
    case <- cases[[nm]]
    dd <- d
    dd$y <- case$y

    fit <- bartisan(y ~ ., data = dd, family = case$family,
                    control = quick_control())

    expect_s3_class(fit, "bartisan_fit")
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
  fit <- bartisan(y ~ ., data = d, family = ordinal(), control = quick_control())

  expect_identical(fit[["num_cat"]], 4L)
  expect_identical(fit[["num_forest"]], 1L)
  # One cutpoint per boundary. With three or more categories all of them are
  # free and the predictor is centered instead; see ?bartisan-families.
  expect_identical(colnames(fit[["aux"]]), paste0("cut", 1:3))
  expect_true(all(abs(rowMeans(fit[["eta"]][[1L]])) < 1e-8))
  expect_gt(stats::sd(fit[["aux"]][, "cut1"]), 0)
  # Cutpoints stay ordered, which the sampler enforces by construction.
  expect_true(all(fit[["aux"]][, "cut2"] < fit[["aux"]][, "cut3"]))
  expect_predictor_invariant(fit, d)

  d$y <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  fit <- bartisan(y ~ ., data = d, family = multinomial(),
                  control = quick_control())

  expect_identical(fit[["num_forest"]], 3L)
  expect_identical(names(fit[["eta"]]), c("a", "b", "c"))
  expect_predictor_invariant(fit, d)

  fit <- bartisan(y ~ ., data = d, family = multinomial(reference = "a"),
                  control = quick_control())

  expect_identical(fit[["num_forest"]], 2L)
  expect_identical(names(fit[["eta"]]), c("b", "c"))
  expect_predictor_invariant(fit, d)
})

test_that("a two-category ordinal model reduces to binary regression", {
  d <- sim_x(seed = 3)
  d$y <- factor(stats::rbinom(nrow(d), 1, 0.5), ordered = TRUE)

  fit <- bartisan(y ~ ., data = d, family = ordinal(), control = quick_control())

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
    fit <- bartisan(survival::Surv(time, status) ~ x1 + x2 + x3, data = d,
                    family = fam, control = quick_control())

    expect_identical(colnames(fit[["aux"]]), "sigma")
    expect_true(all(fit[["aux"]][, "sigma"] > 0))
    expect_predictor_invariant(fit, d)
  }
})

test_that("hard and soft decision rules both fit", {
  d <- sim_x(seed = 5)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  hard <- bartisan(y ~ ., data = d, control = quick_control(gate = "hard"))
  soft <- bartisan(y ~ ., data = d, control = quick_control(gate = "smoothstep"))

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

  fit <- bartisan(y ~ ., data = d, control = quick_control())

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

  fit <- bartisan(y ~ x1 + x2 + x3 + offset(log(logn)), data = d,
                  family = poisson(), control = quick_control())
  expect_s3_class(fit, "bartisan_fit")

  # An observation-level offset is not recoverable from the predictors, so
  # predicting new data must insist on being given one rather than quietly
  # returning predictions that omit it.
  expect_error(predict(fit, newdata = d), "requires .*offset")
  expect_length(predict(fit, newdata = d, offset = log(d$logn)), n)

  # Predictions on the fitting data still work, since the offset is already
  # baked into the stored predictor.
  expect_length(predict(fit), n)

  sub <- bartisan(y ~ x1 + x2 + x3, data = d, subset = seq_len(40),
                  family = poisson(), control = quick_control())
  expect_identical(sub[["n"]], 40L)

  wt <- bartisan(y ~ x1 + x2 + x3, data = d, weights = rep(2, n),
                 family = poisson(), control = quick_control())
  expect_s3_class(wt, "bartisan_fit")
})

test_that("degenerate inputs are rejected", {
  d <- sim_x(seed = 8)
  d$y <- stats::rnorm(nrow(d))

  expect_error(bartisan(y ~ 1, data = d, control = quick_control()),
               "at least one predictor")

  d$constant <- 1
  expect_warning(bartisan(y ~ x1 + constant, data = d,
                          control = quick_control()),
                 "constant predictor")

  d2 <- d
  d2$y <- rep(1, nrow(d2))
  expect_error(bartisan(y ~ x1, data = d2, control = quick_control()),
               "no variation")

  d3 <- d
  d3$y <- -1
  # Matched on a fragment that cannot span a line break in the message source:
  # testthat stops cli reflowing condition messages, so a regex crossing one
  # fails under R CMD check and passes everywhere else.
  expect_error(bartisan(y ~ x1, data = d3, family = poisson(),
                       control = quick_control()),
               "requires a response of non-negative")
  expect_error(bartisan(y ~ x1, data = d3, family = stats::Gamma("log"),
                       control = quick_control()),
               "strictly positive")
})

test_that("print and summary report the fit", {
  d <- sim_x(seed = 9)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- bartisan(y ~ ., data = d, control = quick_control())

  expect_output(print(fit), "Generalized BART")
  expect_output(print(fit), "soft decision rules")

  s <- summary(fit)
  expect_s3_class(s, "summary.bartisan_fit")
  expect_identical(rownames(s[["usage"]][["eta"]]), c("x1", "x2", "x3")[
    order(colMeans(fit[["counts"]][["eta"]] > 0), decreasing = TRUE)])
  expect_true(all(s[["usage"]][["eta"]][, "prop_used"] <= 1))
  expect_output(print(s), "Predictor usage")

  expect_error(summary(fit, level = 1.5), "must be between")
})

test_that("the residual scale prior reads the weights it was given", {
  # `residual_scale()` anchors the prior on sigma. It is a property of the data
  # the model is actually being fitted to, so the weights have to reach it.
  d <- sim_x(n = 300, seed = 141)
  set.seed(1141)
  y <- 2 * d$x1 + stats::rnorm(nrow(d), sd = 0.5)
  x <- as.matrix(d)

  unweighted <- residual_scale(y, x)

  # Only the relative sizes can matter, so a constant weight changes nothing --
  # including the unit weights the unweighted call implies.
  expect_identical(residual_scale(y, x, rep.int(1, length(y))), unweighted)
  expect_equal(residual_scale(y, x, rep.int(7, length(y))), unweighted)

  # A zero weight drops a row rather than shrinking the scale, so a fit that
  # zeroes out a noisy half reports the clean half's scale.
  noisy <- y
  half <- seq_len(150L)
  noisy[half] <- noisy[half] + stats::rnorm(150L, sd = 5)
  keep <- c(rep.int(0, 150L), rep.int(1, 150L))

  expect_equal(residual_scale(noisy, x, keep),
               residual_scale(noisy[-half], x[-half, ]))

  # And ignoring the weights would have been badly wrong here.
  expect_gt(residual_scale(noisy, x), 5 * residual_scale(noisy, x, keep))
})

test_that("the family is read off the response when none is named", {
  d <- sim_x(n = 150, seed = 231)
  set.seed(1231)

  # expect_message() returns the condition, not the value, so the fit is caught
  # by assignment inside it.
  infer <- function(y) {
    d$y <- y
    fit <- NULL
    # Matched on the family expression rather than the leading word, which
    # `arg` sentence-cases.
    expect_message(fit <- bartisan(y ~ x1 + x2, d, control = quick_control()),
                   "family = ")
    c(fit[["family"]][["family"]], fit[["family"]][["link"]])
  }

  expect_identical(infer(stats::rbinom(nrow(d), 1L, 0.5)),
                   c("binomial", "logit"))
  expect_identical(infer(stats::rbinom(nrow(d), 1L, 0.5) == 1L),
                   c("binomial", "logit"))
  expect_identical(infer(factor(sample(c("a", "b"), nrow(d), TRUE))),
                   c("binomial", "logit"))
  expect_identical(infer(factor(sample(c("a", "b", "c"), nrow(d), TRUE))),
                   c("multinomial", "logit"))
  expect_identical(infer(sample(c("a", "b", "c"), nrow(d), TRUE)),
                   c("multinomial", "logit"))
  expect_identical(infer(factor(sample(c("lo", "hi"), nrow(d), TRUE),
                                levels = c("lo", "hi"), ordered = TRUE)),
                   c("ordinal", "logit"))

  # A continuous numeric response gets the mixture, which matches `gaussian()`
  # when the errors are normal and beats it otherwise.
  expect_identical(infer(stats::rnorm(nrow(d))), c("dpm", "identity"))

  # A count is left to the numeric default rather than read as Poisson, and a
  # two-valued numeric whose values are not zero and one is not read as binomial:
  # guessing either would be a modeling decision rather than a reading of the
  # response's type.
  expect_identical(infer(stats::rpois(nrow(d), 30)), c("dpm", "identity"))
  expect_identical(infer(sample(c(1, 2), nrow(d), TRUE)),
                   c("dpm", "identity"))

  # The number of distinct values does not enter into it. A response with only a
  # handful of them probably wants `ordinal()`, but that is a modeling decision
  # and not something to switch families over: a threshold there would make the
  # default arbitrary and hard to predict.
  expect_identical(infer(sample(c(1, 2, 3), nrow(d), TRUE)),
                   c("dpm", "identity"))
  expect_identical(infer(sample(seq_len(9), nrow(d), TRUE)),
                   c("dpm", "identity"))
  expect_identical(infer(sample(seq_len(10), nrow(d), TRUE)),
                   c("dpm", "identity"))
})

test_that("weights with no family named is an error rather than a substitution", {
  d <- sim_x(n = 120, seed = 233)
  set.seed(1233)
  d$y <- stats::rnorm(nrow(d))
  w <- stats::runif(nrow(d), 0.5, 2)

  # The numeric default is `dpm()`, which refuses weights, and only the caller
  # can say whether they meant to drop the weights or to change family.
  expect_error(bartisan(y ~ x1 + x2, d, weights = w, control = quick_control()),
               "does not take prior weights")

  # Naming a family that takes them works, whichever family it is.
  expect_no_error(bartisan(y ~ x1 + x2, d, family = gaussian(), weights = w,
                          control = quick_control()))

  d$coarse <- sample(c(1, 2, 3), nrow(d), TRUE)
  expect_no_error(bartisan(coarse ~ x1 + x2, d, family = ordinal(), weights = w,
                          control = quick_control()))

  # A response whose default is not the mixture is unaffected: only the numeric
  # default runs into this.
  d$bin <- stats::rbinom(nrow(d), 1L, 0.4)
  expect_message(bartisan(bin ~ x1 + x2, d, weights = w,
                          control = quick_control()),
                 "family = ")
})

test_that("a survival response is read as an accelerated failure time model", {
  skip_if_not_installed("survival")

  d <- sim_x(n = 150, seed = 233)
  set.seed(1233)
  d$time <- stats::rexp(nrow(d), 0.5)
  d$event <- stats::rbinom(nrow(d), 1L, 0.8)

  # `dpm_aft()` rather than one of the parametric accelerated failure time
  # families: it was the most accurate of the six over the truths compared in
  # `vignette("survival")` and is cheaper to fit than most of them.
  expect_message(
    fit <- bartisan(survival::Surv(time, event) ~ x1 + x2, d,
                    control = quick_control()),
    "dpm_aft")
  expect_identical(fit[["family"]][["family"]], "dpm_aft")
  expect_identical(fit[["family"]][["link"]], "identity")

  # The same thing as a plain two-column matrix, which the survival families
  # also accept.
  expect_message(
    matrix_fit <- bartisan(cbind(time, event) ~ x1 + x2, d,
                           control = quick_control()),
    "dpm_aft")
  expect_identical(matrix_fit[["family"]][["family"]], "dpm_aft")

  # Naming a parametric family still gets one.
  expect_silent(
    named <- bartisan(survival::Surv(time, event) ~ x1 + x2, d,
                      family = weibull_aft(), control = quick_control()))
  expect_identical(named[["family"]][["link"]], "weibull")

  # `dpm_aft()` takes no prior weights, so a weighted fit that names no family
  # is an error naming the alternatives rather than a silent substitution --
  # the same rule `dpm()` carries for a numeric response.
  expect_error(
    bartisan(survival::Surv(time, event) ~ x1 + x2, d,
             weights = rep(1.5, nrow(d)), control = quick_control()),
    "does not take prior weights")

  # Two columns of counts are a binomial response, not a survival one: the
  # second column is not an event indicator.
  d$hits <- stats::rbinom(nrow(d), 8L, 0.4)
  d$misses <- 8L - d$hits
  expect_message(
    counted <- bartisan(cbind(hits, misses) ~ x1 + x2, d,
                        control = quick_control()),
    "binomial")
  expect_identical(counted[["family"]][["family"]], "binomial")
})

test_that("naming the family silences the message", {
  d <- sim_x(n = 100, seed = 235)
  set.seed(1235)
  d$y <- stats::rnorm(nrow(d))

  expect_silent(bartisan(y ~ x1 + x2, d, family = stats::gaussian(),
                        control = quick_control()))
  expect_silent(bartisan(y ~ x1 + x2, d, family = "gaussian",
                         control = quick_control()))
})
