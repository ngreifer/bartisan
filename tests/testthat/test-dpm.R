# The Dirichlet process mixture for the error distribution. What it adds over
# gaussian() is a *shape* for the errors, so the tests are about the mixture
# being a real posterior over shapes -- and about the reporting chart, which puts
# the mixture at mean zero so that the predictor is the conditional mean.

test_that("the mixture is stored, weighted and shaped as a mixture", {
  d <- sim_x(n = 400, seed = 301)
  set.seed(1301)
  d$y <- 3 * sin(pi * d$x1) + 2 * stats::rt(nrow(d), 3)

  fit <- bartisan(y ~ ., d, family = dpm(),
                  control = quick_control(num_trees = 20L, num_burn = 200L,
                                          num_save = 200L))

  expect_identical(colnames(fit[["aux"]]),
                   c("alpha", "clusters", "center", "error_sd"))
  expect_true(all(fit[["aux"]][, "alpha"] > 0))
  expect_true(all(fit[["aux"]][, "error_sd"] > 0))

  # More than one component, and never more than one per observation.
  clusters <- fit[["aux"]][, "clusters"]
  expect_true(all(clusters >= 1 & clusters <= nrow(d)))
  expect_gt(max(clusters), 1)

  # The components are stored flat with one offset per draw, because the count
  # changes from draw to draw. Their weights are a probability distribution and
  # their number matches what `aux` reports.
  starts <- fit[["mixture_start"]]
  expect_length(starts, nrow(fit[["aux"]]) + 1L)
  expect_identical(diff(starts) / 3, clusters)

  for (s in c(1L, 100L, 200L)) {
    components <- mixture_at(fit, s)
    expect_identical(nrow(components), as.integer(clusters[s]))
    expect_equal(sum(components[, "weight"]), 1)
    expect_true(all(components[, "sd"] > 0))
  }
})

test_that("the reporting chart puts the mixture at zero and the mean on the predictor", {
  d <- data.frame(x = stats::runif(500, -1, 1))
  truth <- 10 * d$x^3
  set.seed(1303)
  # A skewed error, so that the raw mixture is a long way off centre and the
  # shift being taken out is a real one rather than rounding.
  d$y <- truth + 3 * (stats::rgamma(nrow(d), 1.5, 1.5) - 1)

  fit <- bartisan(y ~ x, d, family = dpm(),
                  control = quick_control(num_trees = 50L, num_burn = 400L,
                                          num_save = 400L))

  # The sampler works in a chart where nothing forces the mixture to be centred;
  # reporting is done in the one where it is, so the error mean is zero exactly
  # and the whole conditional mean sits on the predictor.
  on_link <- stats::predict(fit, type = "link")
  on_response <- stats::predict(fit, type = "response")
  expect_identical(on_response, on_link)

  # `center` is the shift that was taken out, so it is not zero when the error is
  # skewed -- that is the quantity, not an estimate of an error mean.
  expect_gt(stats::sd(fit[["aux"]][, "center"]), 0)

  # The level of the predictor no longer trades off against anything, which is
  # the point: before centring its standard deviation across draws was the size
  # of the drift rather than of the posterior.
  level <- rowMeans(stats::predict(fit, type = "link", draws = TRUE))
  expect_lt(stats::sd(level), 0.5)
  expect_lt(abs(mean(on_link - truth)), 0.5)
  expect_lt(sqrt(mean((on_link - truth)^2)), 1.5)
})

test_that("`error_density()` reports a centred density", {
  d <- sim_x(n = 400, seed = 307)
  set.seed(1307)
  d$y <- 2 * d$x1 + 2 * (stats::rgamma(nrow(d), 1.5, 1.5) - 1)

  fit <- bartisan(y ~ ., d, family = dpm(),
                  control = quick_control(num_trees = 20L, num_burn = 300L,
                                          num_save = 300L))

  grid <- seq(-12, 12, length.out = 1201L)
  density <- error_density(fit, at = grid)[, "mean"]
  step <- diff(grid)[1L]

  expect_equal(sum(density) * step, 1, tolerance = 1e-3)
  expect_equal(sum(grid * density) * step, 0, tolerance = 0.05)
  # The shape is what the family is for, so the skewness should survive.
  expect_gt(sum(grid^3 * density) * step, 0.5)
})

test_that("the reported likelihood is the mixture's own predictive", {
  d <- sim_x(n = 250, seed = 305)
  set.seed(1305)
  d$y <- 2 * d$x1 + stats::rnorm(nrow(d))

  fit <- bartisan(y ~ ., d, family = dpm(),
                  control = quick_control(num_trees = 20L, num_burn = 200L,
                                          num_save = 200L))

  # The C++ sums the same expression `predict(type = "density")` builds in R from
  # the stored components, so the two have to agree exactly rather than closely.
  # An error in either would show up here.
  density <- stats::predict(fit, type = "density", draws = TRUE, log = TRUE)
  expect_equal(rowSums(density), fit[["loglik"]], tolerance = 1e-8)

  # And it works on new data, which needs the outcome.
  held <- d[1:25, ]
  expect_length(stats::predict(fit, newdata = held, type = "density"), 25L)
  expect_error(stats::predict(fit, newdata = held[, c("x1", "x2", "x3")],
                              type = "density"),
               "must contain the outcome")
})

test_that("error_density integrates to one and finds the shape of the errors", {
  d <- data.frame(x = stats::runif(600, -1, 1))
  set.seed(1307)
  # Two well-separated normals, which no single normal can look like.
  side <- stats::rbinom(nrow(d), 1L, 0.5)
  d$y <- 5 * d$x + ifelse(side == 1L, stats::rnorm(nrow(d), 3, 0.5),
                          stats::rnorm(nrow(d), -3, 0.5))

  fit <- bartisan(y ~ x, d, family = dpm(),
                  control = quick_control(num_trees = 50L, num_burn = 500L,
                                          num_save = 500L))

  estimated <- error_density(fit, at = seq(-8, 8, length.out = 401L))
  expect_identical(nrow(estimated), 401L)
  expect_named(estimated, c("at", "mean", "lower", "upper"))
  expect_true(all(estimated$mean >= 0))
  expect_true(all(estimated$lower <= estimated$mean))
  expect_true(all(estimated$upper >= estimated$mean))

  step <- diff(estimated$at)[1L]
  expect_equal(sum(estimated$mean) * step, 1, tolerance = 0.05)

  # Bimodal: the density at the centre has to be well below the density at the
  # two modes, which is exactly what a single normal could not produce.
  centre <- estimated$mean[which.min(abs(estimated$at))]
  peaks <- max(estimated$mean)
  expect_lt(centre, 0.5 * peaks)

  # Refused where there is no error distribution to estimate.
  plain <- bartisan(y ~ x, d, family = stats::gaussian(),
                    control = quick_control())
  expect_error(error_density(plain), "estimated error distribution")
})

test_that("the mixture adapts to heavy tails and stays put when the errors are normal", {
  skip_on_cran()

  d <- data.frame(x = stats::runif(800, -1, 1))
  truth <- 10 * d$x^3
  chain <- quick_control(num_trees = 50L, num_burn = 500L, num_save = 500L)

  # Normal errors: the mixture should not need many components, and should agree
  # with the Gaussian family on the scale it reports.
  set.seed(1309)
  d$y <- truth + stats::rnorm(nrow(d), 0, 2)
  normal_fit <- bartisan(y ~ x, d, family = dpm(), control = chain)
  plain <- bartisan(y ~ x, d, family = stats::gaussian(), control = chain)

  expect_equal(stats::sigma(normal_fit), stats::sigma(plain), tolerance = 0.15)
  expect_lt(sqrt(mean((stats::predict(normal_fit, type = "response") -
                         truth)^2)),
            2 * sqrt(mean((stats::predict(plain, type = "response") -
                             truth)^2)))

  # Heavy tails: the estimated density has to be *peaked* relative to the normal
  # with the same standard deviation, which is what a heavy tail looks like from
  # the middle and is the thing the Gaussian family cannot represent.
  set.seed(1310)
  d$y <- truth + 2 * stats::rt(nrow(d), 3)
  heavy <- bartisan(y ~ x, d, family = dpm(), control = chain)

  spread <- stats::sigma(heavy)
  estimated <- error_density(heavy, at = c(0))
  expect_gt(estimated$mean, 1.3 * stats::dnorm(0, 0, spread))
})

test_that("dpm refuses prior weights, and says why", {
  d <- sim_x(n = 150, seed = 311)
  set.seed(1311)
  d$y <- stats::rnorm(nrow(d))
  d$w <- stats::runif(nrow(d), 0.5, 1.5)

  expect_error(bartisan(y ~ x1 + x2, d, family = dpm(), weights = w,
                       control = quick_control()),
               "does not take prior weights")

  # Unit weights are fine, since they are what no weights means.
  d$w <- 1
  expect_no_error(bartisan(y ~ x1 + x2, d, family = dpm(), weights = w,
                          control = quick_control()))
})

test_that("a fixed concentration is honored and the drawn one moves", {
  d <- sim_x(n = 300, seed = 313)
  set.seed(1313)
  d$y <- 2 * d$x1 + 2 * stats::rt(nrow(d), 3)

  chain <- quick_control(num_trees = 20L, num_burn = 200L, num_save = 200L)

  fixed <- bartisan(y ~ ., d, family = dpm(alpha = 2), control = chain)
  expect_true(all(fixed[["aux"]][, "alpha"] == 2))

  drawn <- bartisan(y ~ ., d, family = dpm(), control = chain)
  expect_gt(length(unique(drawn[["aux"]][, "alpha"])), 1L)

  # A larger concentration means more components, which is what it is for.
  many <- bartisan(y ~ ., d, family = dpm(alpha = 50), control = chain)
  few <- bartisan(y ~ ., d, family = dpm(alpha = 0.05), control = chain)
  expect_gt(mean(many[["aux"]][, "clusters"]),
            mean(few[["aux"]][, "clusters"]))
})

test_that("the interop methods work on a mixture fit", {
  skip_if_not_installed("rstantools")

  d <- sim_x(n = 250, seed = 315)
  set.seed(1315)
  d$y <- 2 * d$x1 + 2 * stats::rt(nrow(d), 3)

  fit <- bartisan(y ~ ., d, family = dpm(),
                  control = quick_control(num_trees = 20L, num_burn = 200L,
                                          num_save = 200L))

  # A replicate picks a component -- one of the occupied ones or a fresh draw
  # from the baseline -- so its spread has to be at least the error
  # distribution's and its centre near the fit's.
  replicates <- rstantools::posterior_predict(fit)
  expect_identical(dim(replicates), c(200L, 250L))
  expect_gt(stats::sd(as.vector(replicates)), stats::sigma(fit))
  expect_lt(abs(mean(replicates) - mean(d$y)), 1)

  expect_equal(stats::residuals(fit), d$y - stats::fitted(fit))
  expect_identical(dim(rstantools::log_lik(fit)), c(200L, 250L))
  expect_predictor_invariant(fit, d)
})
