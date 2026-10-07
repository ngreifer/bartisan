# The two routes out of the compiled families: a link supplied from R, composed
# onto the scale the engine's family works on, and a whole log density supplied
# from R. Both call back into the interpreter, so both go through the blocked
# evaluation path, and the first thing to establish is that that path computes
# the same thing as the per-observation one.

test_that("blocked evaluation reproduces the per-observation path", {
  d <- sim_x(n = 120, seed = 31)

  responses <- list(
    gaussian = list(y = 2 * d$x1 + stats::rnorm(120, sd = 0.4),
                    family = gaussian()),
    binomial = list(y = stats::rbinom(120, 1, 0.4), family = binomial()),
    poisson = list(y = stats::rpois(120, 2), family = poisson()),
    # Two additive predictors, so the block carries more than one row per
    # observation and the layout matters.
    gaussian_ls = list(y = 2 * d$x1 + stats::rnorm(120, sd = 0.4),
                          family = gaussian_ls()))

  for (nm in names(responses)) {
    dd <- d
    dd$y <- responses[[nm]][["y"]]

    # Short chains: the two paths differ only in the order the compiler is free
    # to associate a multiply-add, so they agree to rounding error, but a
    # long chain would eventually let a rounding difference flip an
    # accept/reject and diverge for a reason that is not a defect.
    set.seed(5)
    plain <- bartisan(y ~ ., data = dd, family = responses[[nm]][["family"]],
                      control = quick_control(num_burn = 0L, num_draws = 20L))
    set.seed(5)
    blocked <- bartisan(y ~ ., data = dd, family = responses[[nm]][["family"]],
                        control = quick_control(num_burn = 0L, num_draws = 20L,
                                                block_eval = TRUE))

    for (h in seq_along(plain[["eta"]])) {
      expect_equal(plain[["eta"]][[h]], blocked[["eta"]][[h]],
                   tolerance = 1e-10, info = nm)
    }

    expect_equal(plain[["loglik"]], blocked[["loglik"]], tolerance = 1e-10,
                 info = nm)
  }
})

cauchit_opts <- function(derivative = TRUE) {
  out <- list(link_theta = function(eta) stats::qlogis(stats::pcauchy(eta)))

  if (derivative) {
    out[["link_dtheta"]] <- function(eta) {
      p <- stats::pcauchy(eta)
      stats::dcauchy(eta) / (p * (1 - p))
    }
  }

  out
}

test_that("a composed link gives the derivatives its closed form does", {
  set.seed(32)
  n <- 40
  y <- stats::rbinom(n, 1, 0.4)
  eta <- matrix(stats::rnorm(3 * n, 0, 0.8), nrow = 3L)
  y_wide <- matrix(y, nrow = nrow(eta), ncol = n, byrow = TRUE)
  opts <- cauchit_opts()

  # The log likelihood of a cauchit binomial model is
  # y log F(eta) + (1 - y) log(1 - F(eta)) with F the Cauchy cdf, so the score
  # is (y - F) f / (F (1 - F)) and the expected information f^2 / (F (1 - F)).
  p <- stats::pcauchy(eta)
  f <- stats::dcauchy(eta)
  score <- (y_wide - p) * f / (p * (1 - p))

  # Both routes into the family: one observation at a time, and a whole draw
  # through the block methods, which is what the sampler uses.
  for (blocked in c(FALSE, TRUE)) {
    got <- .bartisan_derivs(y, rep(1, n), list(eta), "binomial", "logit", opts,
                            matrix(0, 3L, 0L), 0L, FALSE, blocked)

    expect_equal(got[["d1"]], score, tolerance = 1e-10, info = blocked)

    # The information drops the term in the second derivative of the link,
    # whose expectation is zero, so it is the expected information of the
    # composite and is never negative.
    expect_equal(got[["info"]], f^2 / (p * (1 - p)), tolerance = 1e-10,
                 info = blocked)
    expect_true(all(got[["info"]] >= 0))
  }

  # And the log density itself, which is what the acceptance ratio uses.
  got <- .bartisan_logdens(y, rep(1, n), list(eta), "binomial", "logit", opts,
                           matrix(0, 3L, 0L))
  expect_equal(got, y_wide * log(p) + (1 - y_wide) * log(1 - p),
               tolerance = 1e-12)
})

test_that("a link with no derivative falls back on differences", {
  set.seed(33)
  n <- 30
  y <- stats::rbinom(n, 1, 0.4)
  eta <- matrix(stats::rnorm(2 * n, 0, 0.6), nrow = 2L)

  for (blocked in c(FALSE, TRUE)) {
    with_deriv <- .bartisan_derivs(y, rep(1, n), list(eta), "binomial", "logit",
                                   cauchit_opts(TRUE), matrix(0, 2L, 0L), 0L,
                                  FALSE, blocked)
    without <- .bartisan_derivs(y, rep(1, n), list(eta), "binomial", "logit",
                                cauchit_opts(FALSE), matrix(0, 2L, 0L), 0L,
                               FALSE, blocked)

    expect_equal(with_deriv[["d1"]], without[["d1"]], tolerance = 1e-6)
    expect_equal(with_deriv[["info"]], without[["info"]], tolerance = 1e-6)
  }
})

# A link supplied from R is the caller's code, so what it returns is checked
# before the engine reads past the end of it.
test_that("a link supplied from R that returns the wrong thing is reported", {
  set.seed(34)
  n <- 20
  y <- stats::rbinom(n, 1, 0.4)
  eta <- matrix(stats::rnorm(n), nrow = 1L)

  logdens_with <- function(link) {
    .bartisan_logdens(y, rep(1, n), list(eta), "binomial", "logit",
                      list(link_theta = link), matrix(0, 1L, 0L))
  }

  expect_error(logdens_with(function(eta) 0), "returned 1 values")
  expect_error(logdens_with(function(eta) rep("a", length(eta))),
               "must return a numeric vector")
})

# The same mistake in a link object handed to a family is caught before the fit.
# The sampler calls the link one value at a time, where a function returning a
# single number whatever it is given looks right, so such a link once ran a whole
# fit and then failed in `colMeans()` on the fitted values.
test_that("a link built by the caller is checked before it is used", {
  d <- sim_x(n = 40L, p = 1L, seed = 35L)
  d$y <- stats::rbinom(nrow(d), 1L, 0.4)

  link_with <- function(linkinv) {
    structure(list(linkfun = stats::qlogis, linkinv = linkinv,
                   mu.eta = stats::dlogis, valideta = function(eta) TRUE,
                   name = "mine"),
              class = "link-glm")
  }

  fit_with <- function(linkinv) {
    bartisan(y ~ x1, d, family = stats::binomial(link = link_with(linkinv)),
             control = quick_control(num_burn = 5L, num_draws = 5L))
  }

  expect_error(fit_with(function(eta) stats::plogis(eta)[1L]),
               "returned 1 of type")
  expect_error(fit_with(function(eta) as.character(stats::plogis(eta))),
               "of type <character>")
  expect_s3_class(suppressMessages(fit_with(stats::plogis)), "bartisan_fit")
})

test_that("stats family objects carry an uncompiled link through a fit", {
  d <- sim_x(n = 100, seed = 34)
  truth <- 2 * d$x1 - 1
  d$y <- stats::rbinom(100, 1, stats::pcauchy(truth))

  fit <- bartisan(y ~ ., data = d, family = binomial("cauchit"),
                  control = quick_control())

  expect_identical(fit[["family"]][["link"]], "cauchit")
  expect_predictor_invariant(fit, d)

  # The response scale is the caller's inverse link, not the engine's. Compared
  # draw by draw, since averaging and applying the link do not commute.
  probs <- predict(fit, type = "response", draws = TRUE)
  expect_equal(probs,
               stats::pcauchy(predict(fit, type = "link", draws = TRUE)))
  expect_true(all(probs > 0 & probs < 1))

  # And a link the caller invents, supplied the way glm() takes one.
  invented <- structure(
    list(linkfun = function(mu) stats::qnorm(mu) / 2,
         linkinv = function(eta) stats::pnorm(2 * eta),
         mu.eta = function(eta) 2 * stats::dnorm(2 * eta),
         valideta = function(eta) TRUE, name = "half-probit"),
    class = "link-glm")

  fit <- bartisan(y ~ ., data = d, family = binomial(invented),
                  control = quick_control())

  expect_identical(fit[["family"]][["link"]], "half-probit")
  expect_equal(predict(fit, type = "response", draws = TRUE),
               stats::pnorm(2 * predict(fit, type = "link", draws = TRUE)))
  expect_predictor_invariant(fit, d)
})

test_that("a nuisance parameter is still drawn under a composed link", {
  d <- sim_x(n = 120, seed = 35)
  d$y <- exp(1 + d$x1) + stats::rnorm(120, sd = 0.3)

  fit <- bartisan(y ~ ., data = d, family = gaussian("log"),
                  control = quick_control())

  # The residual scale belongs to the wrapped Gaussian family and is drawn on
  # its own scale, which the wrapper has to leave alone.
  expect_identical(colnames(fit[["aux"]]), "sigma")
  expect_true(all(fit[["aux"]][, "sigma"] > 0))
  expect_predictor_invariant(fit, d)
})

test_that("custom_family reproduces the family it imitates", {
  set.seed(36)
  n <- 50
  y <- stats::rpois(n, 2)
  eta <- matrix(stats::rnorm(4 * n, 0.7, 0.4), nrow = 4L)

  # The Poisson log density without the term free of eta, which the sampler
  # never needs and a custom family is under no obligation to supply.
  opts <- list(num_predictors = 1L,
               logdens = function(y, eta) y * eta[, 1L] - exp(eta[, 1L]),
               derivatives = NULL, name = "poisson-by-hand")

  compiled <- .bartisan_logdens(y, rep(1, n), list(eta), "poisson", "log",
                                list(), matrix(0, 4L, 0L))
  supplied <- .bartisan_logdens(y, rep(1, n), list(eta), "custom", "identity",
                                opts, matrix(0, 4L, 0L))

  expect_equal(compiled, supplied - rep(lgamma(y + 1), each = 4L),
               tolerance = 1e-12)

  # The score of a Poisson model is y - exp(eta) and the information exp(eta);
  # the differences the custom family takes should land on both.
  a <- .bartisan_derivs(y, rep(1, n), list(eta), "poisson", "log", list(),
                        matrix(0, 4L, 0L), 0L, FALSE, FALSE)
  b <- .bartisan_derivs(y, rep(1, n), list(eta), "custom", "identity", opts,
                        matrix(0, 4L, 0L), 0L, FALSE, TRUE)

  expect_equal(a[["d1"]], b[["d1"]], tolerance = 1e-6)
  expect_equal(a[["info"]], b[["info"]], tolerance = 1e-4)
})

test_that("supplied derivatives are used and must be the right shape", {
  set.seed(37)
  n <- 40
  y <- stats::rpois(n, 2)
  eta <- matrix(stats::rnorm(2 * n, 0.5, 0.3), nrow = 2L)
  logdens <- function(y, eta) y * eta[, 1L] - exp(eta[, 1L])

  exact <- .bartisan_derivs(
    y, rep(1, n), list(eta), "custom", "identity",
    list(num_predictors = 1L, logdens = logdens, name = "p",
         derivatives = function(y, eta, h) {
           list(score = y - exp(eta[, h]), info = exp(eta[, h]))
         }),
    matrix(0, 2L, 0L), 0L, FALSE, TRUE)

  differenced <- .bartisan_derivs(
    y, rep(1, n), list(eta), "custom", "identity",
    list(num_predictors = 1L, logdens = logdens, derivatives = NULL,
         name = "p"),
    matrix(0, 2L, 0L), 0L, FALSE, TRUE)

  expect_equal(exact[["d1"]], differenced[["d1"]], tolerance = 1e-6)
  expect_equal(exact[["info"]], differenced[["info"]], tolerance = 1e-4)

  # The same for supplied derivatives that return the wrong number of values.
  expect_error(
    .bartisan_derivs(
      y, rep(1, n), list(eta), "custom", "identity",
      list(num_predictors = 1L, logdens = logdens, name = "p",
           derivatives = function(y, eta, h) {
             list(score = (y - exp(eta[, h]))[-1L], info = exp(eta[, h]))
           }),
      matrix(0, 2L, 0L), 0L, FALSE, TRUE),
    "one score and one information value per observation")

  # A function that returns the wrong number of values is a mistake worth a
  # message rather than a read past the end of a vector.
  expect_error(
    .bartisan_logdens(y, rep(1, n), list(eta), "custom", "identity",
                      list(num_predictors = 1L, name = "p", derivatives = NULL,
                           logdens = function(y, eta) 0),
                     matrix(0, 2L, 0L)),
    "returned 1 values")

  expect_error(
    .bartisan_logdens(y, rep(1, n), list(eta), "custom", "identity",
                      list(num_predictors = 1L, name = "p", derivatives = NULL,
                           logdens = function(y, eta) rep("a", length(y))),
                     matrix(0, 2L, 0L)),
    "must return a numeric vector")
})

test_that("a broken log density errors cleanly and leaves nothing behind", {
  d <- sim_x(n = 80, seed = 43)
  d$y <- stats::rpois(80, 2)

  broken <- custom_family(function(y, eta) stop("deliberate failure"),
                          start = log(mean(d$y)))

  # Three times, because the interesting failure would be the forest leaking on
  # each attempt; the fit is abandoned mid-sweep, so the trees have to be
  # released by a destructor rather than at the end of the run.
  for (i in 1:3) {
    expect_error(bartisan(y ~ ., data = d, family = broken,
                          control = quick_control()),
                 "deliberate failure")
  }

  # And the session is still usable.
  fit <- bartisan(y ~ ., data = d,
                  family = custom_family(
                   function(y, eta) y * eta[, 1L] - exp(eta[, 1L]),
                   start = log(mean(d$y))),
                 control = quick_control())

  expect_true(all(is.finite(predict(fit))))
})

test_that("supplied derivatives are used by both routes alike", {
  set.seed(44)
  n <- 60
  y <- stats::rpois(n, 2)
  eta <- matrix(stats::rnorm(3 * n, 0.5, 0.3), nrow = 3L)

  opts <- list(num_predictors = 1L, name = "p",
               logdens = function(y, eta) y * eta[, 1L] - exp(eta[, 1L]),
               derivatives = function(y, eta, h) {
                 list(score = y - exp(eta[, h]), info = exp(eta[, h]))
               })

  per_unit <- .bartisan_derivs(y, rep(1, n), list(eta), "custom", "identity",
                               opts, matrix(0, 3L, 0L), 0L, FALSE, FALSE)
  blocked <- .bartisan_derivs(y, rep(1, n), list(eta), "custom", "identity",
                              opts, matrix(0, 3L, 0L), 0L, FALSE, TRUE)

  expect_identical(per_unit, blocked)
})

test_that("custom_family fits, predicts and validates its arguments", {
  d <- sim_x(n = 120, seed = 38)
  d$y <- stats::rpois(120, exp(1 + d$x1))

  fit <- bartisan(y ~ ., data = d,
                  family = custom_family(
                   function(y, eta) y * eta[, 1L] - exp(eta[, 1L]),
                   start = log(mean(d$y))),
                 control = quick_control())

  expect_identical(fit[["family"]][["family"]], "custom")
  expect_identical(names(fit[["eta"]]), "eta")
  expect_predictor_invariant(fit, d)

  # There is no mean to report for a density the package cannot interpret, so
  # the response scale is the predictor itself.
  expect_equal(predict(fit, type = "response"), predict(fit, type = "link"))
  expect_error(predict(fit, type = "prob"), "`type` should be one of")

  expect_error(custom_family("not a function"), "must be a function")
  expect_error(custom_family(function(y, eta) y, num_predictors = 0))
  expect_error(custom_family(function(y, eta) y, derivatives = 1),
               "must be a function")
  expect_error(custom_family(function(y, eta) y, num_predictors = 2L,
                             start = c(1, 2, 3)),
               "one per additive predictor")
})

test_that("a custom family with two predictors recovers both surfaces", {
  skip_on_cran()

  set.seed(39)
  n <- 400
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  mean_true <- 2 * d$x1
  log_sd_true <- -1 + d$x2
  d$y <- mean_true + stats::rnorm(n, sd = exp(log_sd_true))

  fit <- bartisan(y ~ ., data = d,
                  family = custom_family(
                   function(y, eta) {
                     stats::dnorm(y, eta[, 1L], exp(eta[, 2L]), log = TRUE)
                   },
                   num_predictors = 2L, start = c(0, -1),
                   name = "location-scale by hand"),
                 control = bartisan_control(num_trees = 20, num_burn = 300,
                                            num_draws = 300, verbose = FALSE))

  eta <- predict(fit, type = "link")
  expect_identical(colnames(eta), c("eta1", "eta2"))
  expect_gt(stats::cor(eta[, 1L], mean_true), 0.9)
  expect_gt(stats::cor(eta[, 2L], log_sd_true), 0.6)
})

test_that("the reported density is right for both R-supplied routes", {
  d <- sim_x(n = 100, seed = 40)
  truth <- 2 * d$x1 - 1

  # A composed link. The reported density is the average over draws of the
  # density at each draw, so it has to be compared that way round.
  d$y <- stats::rbinom(100, 1, stats::pcauchy(truth))
  fit <- bartisan(y ~ ., data = d, family = binomial("cauchit"),
                  control = quick_control())

  p <- predict(fit, type = "response", draws = TRUE)
  y_wide <- matrix(d$y, nrow(p), ncol(p), byrow = TRUE)

  # The two routes reach the probability through different arithmetic -- the
  # engine composes the link inside the log density -- so they agree to rounding
  # error rather than exactly.
  expect_equal(predict(fit, newdata = d, type = "density"),
               colMeans(ifelse(y_wide == 1, p, 1 - p)),
               tolerance = 1e-10)

  # A composed link over a family with a nuisance parameter, which exercises the
  # wrapper's delegation of the parameter and of the eta-free terms.
  d$y <- exp(1 + d$x1) + stats::rnorm(100, sd = 0.3)
  fit <- bartisan(y ~ ., data = d, family = gaussian("log"),
                  control = quick_control())

  e <- predict(fit, type = "link", draws = TRUE)
  y_wide <- matrix(d$y, nrow(e), ncol(e), byrow = TRUE)

  expect_equal(predict(fit, newdata = d, type = "density"),
               colMeans(stats::dnorm(y_wide, exp(e), fit[["aux"]][, "sigma"])),
               tolerance = 1e-8)

  # A custom family. Its log density omits the term free of eta, so what comes
  # back is the Poisson density times y factorial, which is what was asked for.
  d$y <- stats::rpois(100, exp(truth + 2))
  fit <- bartisan(y ~ ., data = d,
                  family = custom_family(
                   function(y, eta) y * eta[, 1L] - exp(eta[, 1L]),
                   start = log(mean(d$y))),
                 control = quick_control())

  e <- predict(fit, type = "link", draws = TRUE)
  y_wide <- matrix(d$y, nrow(e), ncol(e), byrow = TRUE)

  expect_equal(predict(fit, newdata = d, type = "density"),
               colMeans(exp(y_wide * e - exp(e))), tolerance = 1e-8)
})

test_that("a separable response is flagged rather than passed off silently", {
  skip_on_cran()

  set.seed(42)
  d <- data.frame(x1 = stats::runif(200), x2 = stats::runif(200))
  d$y <- as.numeric(d$x1 > 0.5)

  expect_warning(bartisan(y ~ ., data = d, family = binomial(),
                         control = quick_control(num_burn = 200L,
                                                 num_draws = 200L)),
                 "close to separable")

  # And it does not fire on a response the predictors do not separate.
  d$y <- stats::rbinom(200, 1, stats::plogis(3 * (d$x1 - 0.5)))
  expect_no_warning(bartisan(y ~ ., data = d, family = binomial(),
                            control = quick_control(num_burn = 200L,
                                                    num_draws = 200L)))
})

test_that("a pinned leaf scale is never flagged as separable", {
  set.seed(42)
  d <- data.frame(x1 = stats::runif(100), x2 = stats::runif(100))
  d$y <- as.numeric(d$x1 > 0.5)

  # Pinning the leaf scale is the documented remedy, and it silences the
  # warning because there is no longer a drawn scale to run away.
  expect_no_warning(bartisan(y ~ ., data = d, family = binomial(),
                            control = quick_control(update_sigma_mu = FALSE)))
})

test_that("a fit holding R closures survives a round trip through a file", {
  d <- sim_x(n = 100, seed = 45)
  truth <- 2 * d$x1 - 1

  d$y <- stats::rbinom(100, 1, stats::pcauchy(truth))
  linked <- bartisan(y ~ ., data = d, family = binomial("cauchit"),
                     control = quick_control())

  d2 <- d
  d2$y <- stats::rpois(100, exp(truth + 2))
  supplied <- bartisan(y ~ ., data = d2,
                       family = custom_family(
                        function(y, eta) y * eta[, 1L] - exp(eta[, 1L]),
                        start = log(mean(d2$y))),
                      control = quick_control())

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(linked = linked, supplied = supplied), path)
  back <- readRDS(path)

  # The link and the log density are R functions carried on the fit object, so
  # everything that reaches back into them has to still work.
  expect_equal(predict(back[["linked"]], newdata = d, type = "response"),
               predict(linked, newdata = d, type = "response"))
  expect_equal(predict(back[["linked"]], newdata = d, type = "density"),
               predict(linked, newdata = d, type = "density"))
  expect_equal(predict(back[["supplied"]], newdata = d2, type = "density"),
               predict(supplied, newdata = d2, type = "density"))
})

# A nuisance parameter for a user-written likelihood is carried as an additive
# predictor whose forest the engine pins at depth zero: one tree that can never
# split, so the forest is a single scalar drawn by the same Laplace-plus-
# Metropolis step as any leaf. The caller never sees that -- `logdens` is handed
# `aux` as a plain vector and the draws come back in `fit$aux` -- and these tests
# are mostly about that separation holding.

test_that("a nuisance parameter is declared, drawn, and reported as one", {
  skip_on_cran()

  set.seed(151)
  n <- 800
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  x3 = stats::runif(n))
  truth <- 1 + 1.5 * sin(pi * d$x1)
  sigma <- 0.4
  d$y <- stats::rnorm(n, truth, sigma)

  # A Gaussian written out by hand, with its scale drawn rather than fixed in the
  # closure. Deliberately started at 3, seven times the truth, because the point
  # of drawing it is that the caller does not know it.
  by_hand <- custom_family(
    logdens = function(y, eta, aux) {
      stats::dnorm(y, eta[, 1], exp(aux[1]), log = TRUE)
    },
    aux_names = "log_sigma", aux_start = log(3), start = mean(d$y))

  fit <- bartisan(y ~ ., d, family = by_hand,
                  control = bartisan_control(num_burn = 500, num_draws = 1000,
                                             verbose = FALSE))

  # The interface keeps the two kinds apart: one additive predictor, one
  # nuisance parameter, and the pinned forest is not reported as either a
  # predictor or a leaf scale.
  expect_length(fit[["eta"]], 1L)
  expect_identical(colnames(fit[["aux"]]), "log_sigma")
  expect_identical(ncol(fit[["sigma_mu"]]), 1L)

  drawn <- exp(fit[["aux"]][, "log_sigma"])

  # It has to reach the truth from a start seven times too high, which is the
  # regression test for the proposal trust region: an undamped independence
  # proposal cannot be accepted from that far out and the chain sat at its
  # starting value forever.
  expect_equal(mean(drawn), sigma, tolerance = 0.1)
  expect_gt(stats::sd(drawn), 0)

  # And the posterior it reaches is the right one, not merely the right centre:
  # `gaussian()` draws the same parameter by a conjugate step, so the two should
  # agree on the spread as well, and on sigma / sqrt(2n) from theory.
  reference <- bartisan(y ~ ., d, family = gaussian(),
                        control = bartisan_control(num_burn = 500,
                                                   num_draws = 1000,
                                                 verbose = FALSE))
  expect_equal(stats::sd(drawn), stats::sd(reference[["aux"]][, "sigma"]),
               tolerance = 0.25)
  expect_equal(stats::sd(drawn), sigma / sqrt(2 * n), tolerance = 0.3)
})

test_that("several nuisance parameters are kept in order and named", {
  skip_on_cran()

  set.seed(152)
  n <- 400
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- stats::rnorm(n, 2 * d$x1, 0.5)

  # A t likelihood: a scale and a degrees of freedom, both drawn.
  student <- custom_family(
    logdens = function(y, eta, aux) {
      stats::dt((y - eta[, 1]) / exp(aux[1]), df = exp(aux[2]),
                log = TRUE) - aux[1]
    },
    aux_names = c("log_sigma", "log_df"),
    aux_start = c(log(0.5), log(10)), start = mean(d$y))

  fit <- bartisan(y ~ ., d, family = student,
                  control = bartisan_control(num_trees = 20, num_burn = 300,
                                             num_draws = 300, verbose = FALSE))

  expect_identical(colnames(fit[["aux"]]), c("log_sigma", "log_df"))
  expect_length(fit[["eta"]], 1L)
  expect_true(all(is.finite(fit[["aux"]])))
  expect_equal(mean(exp(fit[["aux"]][, "log_sigma"])), 0.5, tolerance = 0.2)

  # The errors are normal, so the degrees of freedom should not settle small.
  expect_gt(mean(exp(fit[["aux"]][, "log_df"])), 5)

  # `summary()` and the convergence table treat them as parameters.
  expect_output(print(summary(fit)), "log_sigma")
})

test_that("the nuisance parameters reach the density and derivative routes", {
  # These go through the block path, where a paired evaluation stacks two values
  # of one component in a single block -- so a nuisance column is constant only
  # in runs. Reading it once per block instead of once per run made the two
  # halves share a value, which silently flattened the target.
  set.seed(153)
  n <- 12
  y <- stats::rnorm(n, 1, 0.5)
  ld <- function(y, eta, aux) stats::dnorm(y, eta[, 1], exp(aux[1]), log = TRUE)
  opts <- list(num_predictors = 1L, num_aux = 1L, aux_names = "log_sigma",
               logdens = ld, derivatives = NULL, name = "custom")

  # Two draws with *different* nuisance values, which is what the runs are for.
  eta <- list(matrix(1, 2L, n), rbind(rep(log(0.5), n), rep(log(2), n)))
  aux <- matrix(c(log(0.5), log(2)), 2L, 1L)

  got <- .bartisan_logdens(y, rep(1, n), eta, "custom", "identity", opts, aux)
  want <- rbind(stats::dnorm(y, 1, 0.5, log = TRUE),
                stats::dnorm(y, 1, 2, log = TRUE))

  expect_equal(got, want, tolerance = 1e-12)

  # The score with respect to the nuisance, against its closed form.
  for (h in 0:1) {
    d <- .bartisan_derivs(y, rep(1, n), eta, "custom", "identity", opts, aux,
                          as.integer(h), FALSE, FALSE)
    s <- c(0.5, 2)
    r <- (y - 1) / rep(s, each = 1)
    want_score <- if (h == 0) {
      rbind((y - 1) / 0.5^2, (y - 1) / 2^2)
    } else {
      rbind(-1 + ((y - 1) / 0.5)^2, -1 + ((y - 1) / 2)^2)
    }
    expect_equal(d[["d1"]], want_score, tolerance = 1e-5,
                 label = paste("component", h))
  }
})

test_that("nuisance parameters are validated at construction", {
  # A two-argument log density cannot receive them, and saying so at
  # construction is better than a confusing error from inside the sampler.
  expect_error(custom_family(function(y, eta) 0, aux_names = "sigma"),
               "third argument")

  expect_error(custom_family(function(y, eta, aux) 0,
                             aux_names = c("a", "a")),
               "distinct")
  expect_error(custom_family(function(y, eta, aux) 0, aux_names = c("a", "")),
               "distinct")
  expect_error(custom_family(function(y, eta, aux) 0,
                             aux_names = c("a", "b"),
                             aux_start = c(1, 2, 3)),
               "one per nuisance parameter")

  # Naming none is the default and leaves the two-argument form working.
  plain <- custom_family(function(y, eta) 0)
  expect_identical(plain[["num_aux"]], 0L)
  expect_identical(plain[["aux_names"]], character())

  # Giving only starting values names them positionally.
  positional <- custom_family(function(y, eta, aux) 0, aux_start = c(0, 0))
  expect_identical(positional[["aux_names"]], c("aux1", "aux2"))

  # A single starting value is recycled across them.
  recycled <- custom_family(function(y, eta, aux) 0,
                            aux_names = c("a", "b"), aux_start = 1)
  expect_identical(recycled[["aux_start"]], c(1, 1))
})

test_that("a custom fit with nuisance parameters has a pointwise log likelihood", {
  d <- sim_x(n = 80, seed = 154)
  d$y <- stats::rnorm(80, d$x1, 0.5)

  by_hand <- custom_family(
    logdens = function(y, eta, aux) {
      stats::dnorm(y, eta[, 1], exp(aux[1]), log = TRUE)
    },
    aux_names = "log_sigma", aux_start = 0, start = mean(d$y))

  fit <- bartisan(y ~ ., d, family = by_hand, control = quick_control())

  # Each draw's density is the caller's function at that draw's predictor and
  # nuisance value. The engine takes a nuisance parameter as a pinned forest,
  # and these all failed while it was handed them as a matrix of values.
  eta <- fit[["eta"]][[1L]]
  sigma <- exp(fit[["aux"]][, "log_sigma"])
  want <- stats::dnorm(matrix(d$y, nrow(eta), ncol(eta), byrow = TRUE), eta,
                       sigma, log = TRUE)

  got <- predict(fit, type = "density", draws = TRUE, log = TRUE)
  expect_equal(got, want, tolerance = 1e-10, ignore_attr = TRUE)

  # New data takes the same route, with the predictor replayed from the trees.
  got_new <- predict(fit, newdata = d[1:10, ], type = "density", draws = TRUE,
                     log = TRUE)
  expect_equal(got_new, want[, 1:10], tolerance = 1e-6, ignore_attr = TRUE)

  skip_if_not_installed("loo")
  expect_equal(dim(suppressWarnings(loo::loo(fit))[["pointwise"]]), c(80L, 5L))
  expect_equal(dim(suppressWarnings(loo::waic(fit))[["pointwise"]]), c(80L, 3L))
})

test_that("a custom family's response reaches its density as given", {
  skip_if_not_installed("patrick")
  skip_if_not_installed("survival")

  d <- sim_x(n = 80, seed = 155)
  d$y1 <- stats::rnorm(80, d$x1)
  d$y2 <- stats::rnorm(80, d$x2, 0.5)
  d$g <- factor(sample(c("a", "b", "c"), 80, replace = TRUE))
  d$time <- stats::rexp(80, exp(d$x1 - 1))
  d$status <- stats::rbinom(80, 1, 0.8)

  # Each density checks the shape it is handed on every call, so a fit that
  # completes has shown that the shape survives the subsetting in the sampler,
  # paired evaluations included.
  bivariate <- function(y, eta) {
    stopifnot(is.matrix(y), ncol(y) == 2L, nrow(y) == nrow(eta))
    stats::dnorm(y[, 1], eta[, 1], log = TRUE) +
      stats::dnorm(y[, 2], eta[, 2], 0.5, log = TRUE)
  }

  categorical <- function(y, eta) {
    stopifnot(is.factor(y), identical(levels(y), c("a", "b", "c")),
              length(y) == nrow(eta))
    lp <- cbind(0, eta)
    lp[cbind(seq_along(y), as.integer(y))] - log(rowSums(exp(lp)))
  }

  exponential <- function(y, eta) {
    stopifnot(inherits(y, "Surv"), nrow(y) == nrow(eta))
    y[, "status"] * eta[, 1] - y[, "time"] * exp(eta[, 1])
  }

  patrick::with_parameters_test_that(
    "response:",
    {
      fit <- bartisan(formula, data = d,
                      family = custom_family(logdens,
                                             num_predictors = num_predictors),
                      control = quick_control())

      response <- stats::model.response(stats::model.frame(formula, d), "any")
      expect_identical(fit[["y"]], response)

      # The stored densities are the caller's function on the whole response,
      # one draw at a time.
      draw_eta <- function(s) {
        vapply(fit[["eta"]], function(e) e[s, ], numeric(nrow(d)))
      }
      want <- t(vapply(seq_len(nrow(fit[["eta"]][[1L]])),
                       function(s) logdens(response, draw_eta(s)),
                       numeric(nrow(d))))

      got <- predict(fit, type = "density", draws = TRUE, log = TRUE)
      expect_equal(got, want, tolerance = 1e-10, ignore_attr = TRUE)
    },
    patrick::cases(
      matrix = list(formula = cbind(y1, y2) ~ x1 + x2 + x3,
                    logdens = bivariate, num_predictors = 2L),
      factor = list(formula = g ~ x1 + x2 + x3, logdens = categorical,
                    num_predictors = 2L),
      Surv = list(formula = survival::Surv(time, status) ~ x1 + x2 + x3,
                  logdens = exponential, num_predictors = 1L)
    )
  )
})

test_that("a factor response is read on the fit's levels in new data", {
  d <- sim_x(n = 60, seed = 156)
  d$g <- factor(sample(c("a", "b", "c"), 60, replace = TRUE))

  categorical <- function(y, eta) {
    lp <- cbind(0, eta)
    lp[cbind(seq_along(y), as.integer(y))] - log(rowSums(exp(lp)))
  }

  fit <- bartisan(g ~ ., d,
                  family = custom_family(categorical, num_predictors = 2L),
                  control = quick_control())

  # Dropping the unused level would make "b" code 1 and "c" code 2, so a
  # density reading the codes would score each row against the wrong category.
  new <- d[d$g != "a", ][1:5, ]
  new$g <- droplevels(new$g)

  eta <- predict(fit, newdata = new, type = "link", draws = TRUE)
  relevelled <- factor(as.character(new$g), levels = c("a", "b", "c"))
  want <- t(vapply(seq_len(nrow(eta[[1L]])), function(s) {
    categorical(relevelled, cbind(eta[[1L]][s, ], eta[[2L]][s, ]))
  }, numeric(5L)))

  got <- predict(fit, newdata = new, type = "density", draws = TRUE,
                 log = TRUE)
  expect_equal(got, want, tolerance = 1e-10, ignore_attr = TRUE)

  new$g <- factor(c("b", "d", "c", "b", "c"))
  expect_error(predict(fit, newdata = new, type = "density"),
               "not fit with")
})

test_that("supplied derivatives receive the nuisance values", {
  set.seed(157)
  n <- 30
  y <- stats::rnorm(n, 1, 0.5)
  ld <- function(y, eta, aux) {
    stats::dnorm(y, eta[, 1], exp(aux[1]), log = TRUE)
  }
  dv <- function(y, eta, h, aux) {
    v <- exp(2 * aux[1])
    list(score = (y - eta[, h]) / v, info = rep(1 / v, length(y)))
  }

  eta <- list(matrix(stats::rnorm(2 * n, 1, 0.2), 2L, n),
              rbind(rep(log(0.5), n), rep(log(2), n)))
  aux <- matrix(c(log(0.5), log(2)), 2L, 1L)
  opts <- list(num_predictors = 1L, num_aux = 1L, aux_names = "log_sigma",
               logdens = ld, name = "custom")

  exact <- .bartisan_derivs(y, rep(1, n), eta, "custom", "identity",
                            c(opts, list(derivatives = dv)), aux, 0L, FALSE,
                            TRUE)
  differenced <- .bartisan_derivs(y, rep(1, n), eta, "custom", "identity",
                                  c(opts, list(derivatives = NULL)), aux, 0L,
                                  FALSE, TRUE)

  # Two draws with different scales, so a score computed without the nuisance
  # value would match at most one of them.
  expect_equal(exact[["d1"]], differenced[["d1"]], tolerance = 1e-6)
  expect_equal(exact[["info"]], differenced[["info"]], tolerance = 1e-4)

  # And through a fit, where the nuisance value has to be passed on by the R
  # side as well as by the engine.
  d <- sim_x(n = 60, seed = 157)
  d$y <- stats::rnorm(60, d$x1, 0.5)
  fit <- bartisan(y ~ ., d,
                  family = custom_family(ld, derivatives = dv,
                                         aux_names = "log_sigma",
                                         start = mean(d$y)),
                  control = quick_control())
  expect_true(all(is.finite(fit[["aux"]])))

  expect_error(custom_family(ld, aux_names = "log_sigma",
                             derivatives = function(y, eta, h) list()),
               "fourth argument")
})

test_that("a bivariate normal through custom_family recovers its correlation", {
  skip_on_cran()

  set.seed(158)
  n <- 400
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  x3 = stats::runif(n))
  z1 <- stats::rnorm(n)
  z2 <- 0.6 * z1 + sqrt(1 - 0.6^2) * stats::rnorm(n)
  d$y1 <- sin(2 * pi * d$x1) + 0.5 * z1
  d$y2 <- 2 * (d$x2 - 0.5)^2 + 0.3 * z2

  bvn <- custom_family(
    logdens = function(y, eta, aux) {
      s <- exp(aux[1:2])
      rho <- tanh(aux[3])
      z1 <- (y[, 1] - eta[, 1]) / s[1]
      z2 <- (y[, 2] - eta[, 2]) / s[2]
      -sum(log(s)) - log1p(-rho^2) / 2 -
        (z1^2 - 2 * rho * z1 * z2 + z2^2) / (2 * (1 - rho^2))
    },
    num_predictors = 2, start = c(mean(d$y1), mean(d$y2)),
    aux_names = c("log_sd1", "log_sd2", "atanh_rho"))

  fit <- bartisan(cbind(y1, y2) ~ x1 + x2 + x3, d, family = bvn,
                  control = bartisan_control(num_trees = 20, num_burn = 300,
                                             num_draws = 300, verbose = FALSE))

  aux <- fit[["aux"]]
  expect_equal(mean(tanh(aux[, "atanh_rho"])), 0.6, tolerance = 0.1)
  expect_equal(mean(exp(aux[, "log_sd1"])), 0.5, tolerance = 0.15)
  expect_equal(mean(exp(aux[, "log_sd2"])), 0.3, tolerance = 0.15)
  expect_gt(stats::cor(colMeans(fit[["eta"]][[1L]]), sin(2 * pi * d$x1)), 0.95)
})
