# The second shape the sampler can exploit. Hill et al. (2020, sec. 3.1.5)
# observe that the count models share the form c + a*eta + b*exp(s*eta); like the
# quadratic form, three numbers from one pass determine it everywhere, so the
# Laplace fit is iterated on scalars rather than on the data.

test_that("the families that claim the exponential form really have it", {
  # The premise, checked directly: fit the three coefficients from three
  # evaluations and see whether they predict the log density elsewhere. If the
  # form is right this is exact; if it is wrong the prediction is nonsense.
  set.seed(111)
  n <- 40
  probe <- c(-0.4, 0.1, 0.7)
  held_out <- c(-2, -1, 0, 0.5, 1.5, 2.5)

  fitted_error <- function(family, y, opts, aux, sign) {
    at <- function(e) {
      .genbart_logdens(y, rep(1, n),
                       list(matrix(e, nrow = 1L, ncol = n)), family, "log",
                       opts, aux)[1L, ]
    }

    design <- cbind(1, probe, exp(sign * probe))
    coefficients <- t(solve(design, t(vapply(probe, at, numeric(n)))))

    max(vapply(held_out, function(e) {
      predicted <- coefficients[, 1L] + coefficients[, 2L] * e +
        coefficients[, 3L] * exp(sign * e)
      max(abs(predicted - at(e)))
    }, numeric(1)))
  }

  expect_lt(fitted_error("poisson", stats::rpois(n, 3), list(),
                         matrix(0, 1L, 0L), 1), 1e-10)
  expect_lt(fitted_error("Gamma", stats::rgamma(n, 3, 1),
                         list(shape = 3, shape_prior_shape = 0.01,
                              shape_prior_rate = 0.01, update_shape = FALSE),
                         matrix(3, 1L, 1L), -1), 1e-10)

  # A family that is not of this form, so that the check above is known to be
  # capable of failing.
  expect_gt(fitted_error("binomial", stats::rbinom(n, 1, 0.5), list(),
                         matrix(0, 1L, 0L), 1), 0.1)
})

test_that("the exponential shortcut reproduces the general path", {
  d <- sim_x(n = 250, seed = 112)
  linear <- 1 + 1.5 * d$x1

  # The tolerance is not the same for all three, and the reason is the *general*
  # path rather than the closed form. Fisher scoring there stops once the score
  # falls below a hundredth of a standard error, so its mode is off by up to
  # that much, and the closed form -- which iterates to the machine's limit --
  # is the more accurate of the two. Where the general path happens to land on
  # the mode anyway the two agree to the last bits; where it stops short they
  # agree to the tolerance. The augmented negative binomial has the sharpest
  # target of the three and so shows the largest residual.
  settings <- list(
    list(label = "poisson", family = poisson(),
         y = stats::rpois(250, exp(linear)), tol = 1e-10),
    list(label = "Gamma(log)", family = Gamma("log"),
         y = stats::rgamma(250, 3, rate = 3 / exp(linear)), tol = 1e-10),
    list(label = "negbin augmented", family = negbin(),
         y = stats::rnbinom(250, mu = exp(linear), size = 2),
         augment = "negbin", tol = 1e-3))

  for (s in settings) {
    dd <- d
    dd$y <- s$y

    ctrl <- function(exact) {
      genbart_control(num_trees = 10, num_burn = 0L, num_save = 1L,
                      soft = FALSE, verbose = FALSE, exact_quadratic = exact,
                      augment = if (is.null(s$augment)) FALSE else s$augment)
    }

    set.seed(5)
    general <- genbart(y ~ ., data = dd, family = s$family, control = ctrl(FALSE))
    set.seed(5)
    closed <- genbart(y ~ ., data = dd, family = s$family, control = ctrl(TRUE))

    # One sweep, so the two are compared before any accept-or-reject decision
    # has had the chance to amplify a rounding difference.
    expect_lt(max(abs(general[["eta"]][["eta"]] - closed[["eta"]][["eta"]])),
              s$tol, label = s$label)
  }
})

test_that("the exponential shortcut stands down for soft rules", {
  # A soft rule gives each observation its own exponent, and a sum of those is
  # not a function of three numbers. The shortcut has to notice and decline, so
  # the draws must be bitwise identical either way.
  d <- sim_x(n = 200, seed = 113)
  d$y <- stats::rpois(200, exp(1 + d$x1))

  ctrl <- function(exact) {
    genbart_control(num_trees = 10, num_burn = 0L, num_save = 20L,
                    verbose = FALSE, exact_quadratic = exact)
  }

  set.seed(5)
  a <- genbart(y ~ ., data = d, family = poisson(), control = ctrl(FALSE))
  set.seed(5)
  b <- genbart(y ~ ., data = d, family = poisson(), control = ctrl(TRUE))

  expect_identical(a[["eta"]], b[["eta"]])
})

test_that("the statically dispatched accumulators match the virtual ones", {
  # Every family that can reach a closed form has its leaf sums accumulated by a
  # loop against its own concrete type, so the compiler can inline the
  # arithmetic. That is a second implementation of the same sums, and it has to
  # agree to the bit.
  d <- sim_x(n = 200, seed = 114)
  linear <- 1 + 1.5 * d$x1

  settings <- list(
    list("gaussian soft", gaussian(), linear + stats::rnorm(200, sd = 0.4),
         list()),
    list("gaussian hard", gaussian(), linear + stats::rnorm(200, sd = 0.4),
         list(soft = FALSE)),
    list("location-scale", location_scale(),
         linear + stats::rnorm(200, sd = 0.4), list()),
    list("poisson hard", poisson(), stats::rpois(200, exp(linear)),
         list(soft = FALSE)),
    list("Gamma hard", Gamma("log"),
         stats::rgamma(200, 3, rate = 3 / exp(linear)), list(soft = FALSE)),
    list("probit augmented", binomial("probit"),
         stats::rbinom(200, 1, stats::pnorm(linear - 2)), list(augment = TRUE)),
    list("logit augmented", binomial("logit"),
         stats::rbinom(200, 1, stats::plogis(linear - 2)),
         list(augment = TRUE, soft = FALSE)),
    list("negbin augmented", negbin(),
         stats::rnbinom(200, mu = exp(linear), size = 2),
         list(augment = "negbin", soft = FALSE)))

  for (s in settings) {
    dd <- d
    dd$y <- s[[3L]]

    control <- function(generic) {
      args <- c(list(num_trees = 10L, num_burn = 0L, num_save = 25L,
                     verbose = FALSE, generic_accumulate = generic), s[[4L]])
      do.call(genbart_control, args)
    }

    set.seed(5)
    virtual <- genbart(y ~ ., data = dd, family = s[[2L]],
                       control = control(TRUE))
    set.seed(5)
    static <- genbart(y ~ ., data = dd, family = s[[2L]],
                      control = control(FALSE))

    expect_identical(virtual[["eta"]], static[["eta"]], label = s[[1L]])
    expect_identical(virtual[["loglik"]], static[["loglik"]], label = s[[1L]])
  }
})

test_that("the gamma family reports its true curvature", {
  # It was reporting the expected information, which is just the shape. The
  # observed second derivative is shape * y / mu, always positive here because
  # the response is, and it is what the exponential form reads the
  # coefficients off.
  set.seed(115)
  n <- 30
  y <- stats::rgamma(n, 3, 1)
  eta <- matrix(stats::rnorm(2 * n, 0.3, 0.5), nrow = 2L)
  opts <- list(shape = 3, shape_prior_shape = 0.01, shape_prior_rate = 0.01,
               update_shape = FALSE)

  got <- .genbart_derivs(y, rep(1, n), list(eta), "Gamma", "log", opts,
                         matrix(3, 2L, 1L), 0L, FALSE, FALSE)

  y_wide <- matrix(y, nrow = 2L, ncol = n, byrow = TRUE)
  expect_equal(got[["info"]], 3 * y_wide * exp(-eta), tolerance = 1e-12)
  expect_true(all(got[["info"]] > 0))
})

test_that("the augmented negative binomial targets the same posterior", {
  skip_on_cran()

  # The rewriting here is not a Polya-Gamma one: the count is written as a
  # Poisson whose rate is drawn from a gamma, which puts the target in the
  # exponential form and costs one gamma draw per observation. That it is the
  # same posterior is the thing worth checking.
  set.seed(116)
  n <- 500
  X <- matrix(stats::runif(n * 4), n)
  colnames(X) <- paste0("x", 1:4)
  d <- as.data.frame(X)
  truth <- 1.5 + 1.5 * sin(pi * X[, 1]) + 1.5 * (X[, 2] - 0.5)
  d$y <- stats::rnbinom(n, mu = exp(truth), size = 2)

  ctrl <- function(augment) {
    genbart_control(num_trees = 50, num_burn = 1000, num_save = 1000,
                    soft = FALSE, augment = augment, verbose = FALSE)
  }

  set.seed(5)
  direct <- genbart(y ~ ., data = d, family = negbin(), control = ctrl(FALSE))
  set.seed(5)
  rewritten <- genbart(y ~ ., data = d, family = negbin(),
                       control = ctrl("negbin"))

  a <- predict(direct, type = "link")
  b <- predict(rewritten, type = "link")

  expect_gt(stats::cor(a, b), 0.98)
  expect_equal(sqrt(mean((b - truth)^2)), sqrt(mean((a - truth)^2)),
               tolerance = 0.2)

  # The dispersion is a global parameter and the one most likely to shift if the
  # rewriting were wrong.
  expect_equal(mean(rewritten[["aux"]][, "theta"]),
               mean(direct[["aux"]][, "theta"]), tolerance = 0.15)

  # And the reported log likelihood is still the negative binomial one.
  eta <- predict(rewritten, type = "link", draws = TRUE)
  theta <- rewritten[["aux"]][, "theta"]
  by_hand <- vapply(seq_len(nrow(eta)), function(s) {
    sum(stats::dnbinom(d$y, mu = exp(eta[s, ]), size = theta[s], log = TRUE))
  }, numeric(1))

  expect_equal(rewritten[["loglik"]], by_hand, tolerance = 1e-8)
})
