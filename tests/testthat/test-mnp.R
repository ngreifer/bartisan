# Multinomial probit. The model differs from the multinomial logit in one way
# that matters -- the latent variables are correlated -- and in one way that
# constrains everything else: its likelihood is a Gaussian orthant probability
# with no closed form, so the probabilities are simulated.

# Draw a categorical response from the multinomial probit model itself, so that
# the truth is the thing being fitted rather than an approximation to it.
sim_mnp <- function(latent, sigma, seed) {
  set.seed(seed)
  n <- nrow(latent)
  drawn <- latent + matrix(stats::rnorm(n * ncol(latent)), n) %*% chol(sigma)
  winner <- apply(drawn, 1L, which.max)
  code <- ifelse(apply(drawn, 1L, max) < 0, 0L, winner)
  factor(c("ref", "a", "b")[code + 1L], levels = c("ref", "a", "b"))
}

# The category probabilities the truth implies, by simulation.
true_probs <- function(latent, sigma, reps = 4000L, seed = 99L) {
  set.seed(seed)
  n <- nrow(latent)
  out <- matrix(0, n, ncol(latent) + 1L)

  for (r in seq_len(reps)) {
    drawn <- latent + matrix(stats::rnorm(n * ncol(latent)), n) %*% chol(sigma)
    code <- ifelse(apply(drawn, 1L, max) < 0, 0L, apply(drawn, 1L, which.max))
    at <- cbind(seq_len(n), code + 1L)
    out[at] <- out[at] + 1
  }

  colnames(out) <- c("ref", "a", "b")
  out / reps
}

test_that("two categories are binary probit, exactly as the trace constraint says", {
  d <- sim_x(n = 500, seed = 201)
  set.seed(1201)
  d$y <- factor(ifelse(stats::rbinom(nrow(d), 1L,
                                     stats::pnorm(1.5 * sin(pi * d$x1) - d$x2)) == 1L,
                       "yes", "no"),
                levels = c("no", "yes"))

  chain <- quick_control(num_trees = 30L, num_burn = 300L, num_draws = 300L)
  probit <- bartisan(y ~ ., d, family = multinomial("probit"), control = chain)
  binary <- bartisan(y ~ ., d, family = stats::binomial("probit"),
                     control = chain)

  # With one latent variable the trace constraint pins its variance at one, so
  # the two are the same model. They are fitted by different samplers, so they
  # agree up to Monte Carlo error.
  expect_gt(stats::cor(stats::predict(probit, type = "link"),
                       stats::predict(binary, type = "link")),
            0.98)
  expect_lt(max(abs(stats::predict(probit, type = "prob")[, "yes"] -
                      stats::predict(binary, type = "prob")[, "yes"])),
            0.1)

  # And there is no free covariance element to report.
  expect_null(probit[["aux"]])
  expect_identical(probit[["num_forest"]], 1L)
})

test_that("the fit recovers the latent means and the category probabilities", {
  d <- sim_x(n = 800, seed = 203)
  latent <- cbind(1.5 * sin(pi * d$x1), 1.4 * d$x2 - 0.7)
  sigma <- matrix(c(1, 0.5, 0.5, 1), 2L)
  d$y <- sim_mnp(latent, sigma, seed = 1203)

  fit <- bartisan(y ~ ., d, family = multinomial("probit", reference = "ref"),
                  control = quick_control(num_trees = 50L, num_burn = 400L,
                                          num_draws = 400L))

  # One forest per non-reference category, named for the contrast it carries.
  expect_identical(fit[["num_forest"]], 2L)
  expect_identical(names(fit[["eta"]]), c("a-ref", "b-ref"))

  predictor <- stats::predict(fit, type = "link")
  expect_gt(stats::cor(predictor[, 1L], latent[, 1L]), 0.9)
  expect_gt(stats::cor(predictor[, 2L], latent[, 2L]), 0.9)

  probs <- stats::predict(fit, type = "prob")
  truth <- true_probs(latent, sigma)
  expect_identical(colnames(probs), c("ref", "a", "b"))
  expect_equal(rowSums(probs), rep.int(1, nrow(d)))
  expect_lt(sqrt(mean((probs - truth)^2)), 0.06)
  expect_predictor_invariant(fit, d)
})

test_that("the latent covariance respects the trace constraint and finds the correlation", {
  d <- sim_x(n = 900, seed = 205)
  latent <- cbind(1.5 * sin(pi * d$x1), 1.4 * d$x2 - 0.7)

  chain <- quick_control(num_trees = 50L, num_burn = 400L, num_draws = 400L)
  estimated <- vapply(c(-0.8, 0.8), function(rho) {
    sigma <- matrix(c(1, rho, rho, 1), 2L)
    d$y <- sim_mnp(latent, sigma, seed = 1205)
    fit <- bartisan(y ~ ., d, family = multinomial("probit", reference = "ref"),
                    control = chain)

    covariance <- fit[["aux"]]
    expect_identical(colnames(covariance),
                     c("sigma11", "sigma21", "sigma22"))

    # trace(Sigma) = C at every draw, which is the identification.
    expect_equal(covariance[, "sigma11"] + covariance[, "sigma22"],
                 rep.int(2, nrow(covariance)))
    expect_true(all(covariance[, "sigma11"] > 0))
    expect_true(all(covariance[, "sigma22"] > 0))

    # A covariance matrix, so the implied correlation is inside the unit ball.
    correlation <- covariance[, "sigma21"] /
      sqrt(covariance[, "sigma11"] * covariance[, "sigma22"])
    expect_true(all(abs(correlation) < 1))

    mean(correlation)
  }, numeric(1L))

  # Only the deterministic constraints above and the extremes are worth
  # asserting. The correlation enters the likelihood through orthant
  # probabilities of a distribution whose location is a sum of trees, and a
  # nonparametric mean absorbs much of the dependence, so it is weakly
  # identified: measured at this sample size a true correlation of zero came
  # back anywhere from -0.57 to 0.33 depending on the draw of the data, and the
  # sweep over the truth was not monotone. What does hold is that the two
  # extremes are separated and each has the right sign.
  expect_lt(estimated[1L], 0)
  expect_gt(estimated[2L], 0)
  expect_gt(estimated[2L] - estimated[1L], 0.5)
})

test_that("the reported likelihood is the simulated multinomial probit one", {
  d <- sim_x(n = 300, seed = 209)
  latent <- cbind(1.2 * d$x1, -1.1 * d$x2 + 0.4)
  d$y <- sim_mnp(latent, matrix(c(1, 0.3, 0.3, 1), 2L), seed = 1209)

  fit <- bartisan(y ~ ., d, family = multinomial("probit", reference = "ref"),
                  control = quick_control(num_trees = 20L, num_burn = 200L,
                                          num_draws = 200L))

  # Not the augmented Gaussian the sampler works with: a log probability, so
  # every term is negative and the total is on the scale a categorical
  # likelihood is.
  expect_true(all(fit[["loglik"]] < 0))
  expect_true(all(is.finite(fit[["loglik"]])))
  expect_gt(mean(fit[["loglik"]]), nrow(d) * log(1 / 3) * 2)

  # The density of each observation's own category, which has to agree with the
  # probability matrix at the same draws.
  probs <- stats::predict(fit, type = "prob", draws = TRUE)
  density <- stats::predict(fit, type = "density", draws = TRUE)
  observed <- match(as.character(d$y), c("ref", "a", "b"))
  by_hand <- vapply(seq_along(observed),
                    function(i) probs[, i, observed[i]],
                    numeric(dim(probs)[1L]))

  # Both are simulated with fresh draws, so they agree in distribution rather
  # than exactly. A few hundred replicates put the two within a few hundredths.
  expect_lt(max(abs(colMeans(density) - colMeans(by_hand))), 0.05)

  # And it works on new data, which needs the outcome.
  held <- d[1:20, ]
  expect_length(stats::predict(fit, newdata = held, type = "density"), 20L)
  expect_error(stats::predict(fit, newdata = held[, c("x1", "x2", "x3")],
                              type = "density"),
               "must contain the outcome")
})

test_that("the probit target is quadratic, and its derivatives are the analytic ones", {
  d <- sim_x(n = 200, seed = 211)
  latent <- cbind(d$x1, -d$x2)
  d$y <- sim_mnp(latent, matrix(c(1, 0.4, 0.4, 1), 2L), seed = 1211)

  fit <- bartisan(y ~ ., d, family = multinomial("probit", reference = "ref"),
                  control = quick_control())

  # Conditional on the latent variables the target is exactly quadratic in every
  # component, which is what makes the leaf draw a closed form rather than a
  # Laplace approximation. Checked through the same route the other families
  # use: the analytic score and information against central differences of the
  # family's own log density.
  eta <- stats::predict(fit, type = "link", draws = TRUE)

  for (h in seq_len(2L)) {
    analytic <- .bartisan_derivs(y = fit[["y"]],
                                 weights = fit[["prior_weights"]],
                                eta_draws = eta, family_name = "mnp",
                                link = "probit",
                                family_opts = fit[["family_opts"]],
                                aux = fit[["aux"]], component = h - 1L,
                                by_difference = FALSE)
    numeric <- .bartisan_derivs(y = fit[["y"]],
                                weights = fit[["prior_weights"]],
                               eta_draws = eta, family_name = "mnp",
                               link = "probit",
                               family_opts = fit[["family_opts"]],
                               aux = fit[["aux"]], component = h - 1L,
                               by_difference = TRUE)

    expect_equal(analytic$d1, numeric$d1, tolerance = 1e-5)
    expect_equal(analytic$info, numeric$info, tolerance = 1e-5)
  }
})

test_that("the categorical prediction types and the interop methods all work", {
  skip_if_not_installed("rstantools")

  d <- sim_x(n = 300, seed = 213)
  latent <- cbind(1.2 * d$x1, -1.1 * d$x2 + 0.4)
  d$y <- factor(c("1", "2", "4")[
    as.integer(sim_mnp(latent, matrix(c(1, 0.3, 0.3, 1), 2L), seed = 1213))],
    levels = c("1", "2", "4"))

  fit <- bartisan(y ~ ., d, family = multinomial("probit"),
                  control = quick_control(num_trees = 20L, num_burn = 200L,
                                          num_draws = 200L))

  expect_s3_class(stats::predict(fit, type = "class"), "factor")
  expect_identical(dim(stats::predict(fit, type = "response")), c(300L, 3L))

  # The labels are numbers, so the mean of the response is available.
  means <- stats::predict(fit, type = "mean")
  expect_true(all(means > 1 & means < 4))

  # A replicate is a category index, and every category is reachable.
  replicates <- rstantools::posterior_predict(fit)
  expect_true(all(replicates %in% 1:3))
  expect_identical(dim(replicates), c(200L, 300L))

  expect_s3_class(stats::simulate(fit, nsim = 2L)$sim_1, "factor")
  expect_error(stats::residuals(fit), "no mean")
  expect_identical(dim(rstantools::log_lik(fit)), c(200L, 300L))
})

test_that("a covariance is not something augment can be asked for", {
  # The latent variables are the model rather than a rewriting of it, so `mnp`
  # is not one of the names `augment` takes.
  expect_error(bartisan_control(augment = "mnp"), "must be one of")

  d <- sim_x(n = 200, seed = 215)
  latent <- cbind(d$x1, -d$x2)
  d$y <- sim_mnp(latent, diag(2), seed = 1215)

  # And turning augmentation off changes nothing for it.
  chain <- quick_control(num_trees = 20L, num_burn = 200L, num_draws = 200L)
  set.seed(6)
  on <- bartisan(y ~ ., d, family = multinomial("probit"), control = chain)
  set.seed(6)
  off <- bartisan(y ~ ., d, family = multinomial("probit"),
                  control = quick_control(num_trees = 20L, num_burn = 200L,
                                          num_draws = 200L, augment = FALSE))
  expect_equal(on[["eta"]], off[["eta"]])
})
