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
      .bartisan_logdens(y, rep(1, n),
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
    list(label = "Gamma(log)", family = stats::Gamma("log"),
         y = stats::rgamma(250, 3, rate = 3 / exp(linear)), tol = 1e-10),
    list(label = "negbin augmented", family = negbin(),
         y = stats::rnbinom(250, mu = exp(linear), size = 2),
         augment = "negbin", tol = 1e-3))

  # The two families the *rate* unlocked, which the form could not reach while it
  # was fixed at exp(+-eta). Both carry a rate that is not +-1, and the Weibull's
  # moves with the scale parameter from sweep to sweep.
  weibull <- exp(linear + 0.5 * log(-log(stats::runif(250))))
  extra <- list(
    list(label = "weibull AFT, rate -1/sigma",
         family = weibull_aft(), tol = 1e-10,
         # A two-column matrix rather than a Surv object, so that the check does
         # not rest on a suggested package.
         y = cbind(pmin(weibull, 12), as.numeric(weibull <= 12))),
    list(label = "location-scale, rate -2", family = location_scale(),
         y = linear + stats::rnorm(250, sd = exp(-1 + d$x2)), tol = 1e-4),
    # Proportional hazards: the piecewise-exponential likelihood is
    # delta * eta - Lambda_0(y) exp(eta), which is the form at rate +1.
    # The tolerance is the general path's, as for the location-scale row above:
    # Fisher scoring there stops within a hundredth of a standard error, and the
    # closed form iterates to the machine's limit.
    list(label = "proportional hazards, rate +1", family = ph(), tol = 1e-4,
         y = cbind(pmin(stats::rexp(250, exp(linear)), 3),
                   as.numeric(stats::rexp(250, exp(linear)) <= 3))))

  for (s in c(settings, extra)) {
    dd <- d
    dd$y <- s$y

    ctrl <- function(exact) {
      bartisan_control(num_trees = 10, num_burn = 0L, num_draws = 1L,
                       gate = "hard", verbose = FALSE, exact_quadratic = exact,
                      augment = if (is.null(s$augment)) FALSE else s$augment)
    }

    set.seed(5)
    general <- bartisan(y ~ ., data = dd, family = s$family, control = ctrl(FALSE))
    set.seed(5)
    closed <- bartisan(y ~ ., data = dd, family = s$family, control = ctrl(TRUE))

    # One sweep, so the two are compared before any accept-or-reject decision
    # has had the chance to amplify a rounding difference.
    expect_lt(max(abs(unlist(general[["eta"]]) - unlist(closed[["eta"]]))),
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
    bartisan_control(num_trees = 10, num_burn = 0L, num_draws = 20L,
                     verbose = FALSE, exact_quadratic = exact)
  }

  set.seed(5)
  a <- bartisan(y ~ ., data = d, family = poisson(), control = ctrl(FALSE))
  set.seed(5)
  b <- bartisan(y ~ ., data = d, family = poisson(), control = ctrl(TRUE))

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
         list(gate = "hard")),
    list("location-scale", location_scale(),
         linear + stats::rnorm(200, sd = 0.4), list()),
    list("poisson hard", poisson(), stats::rpois(200, exp(linear)),
         list(gate = "hard")),
    list("Gamma hard", stats::Gamma("log"),
         stats::rgamma(200, 3, rate = 3 / exp(linear)), list(gate = "hard")),
    list("probit augmented", binomial("probit"),
         stats::rbinom(200, 1, stats::pnorm(linear - 2)), list(augment = TRUE)),
    list("logit augmented", binomial("logit"),
         stats::rbinom(200, 1, stats::plogis(linear - 2)),
         list(augment = TRUE, gate = "hard")),
    list("negbin augmented", negbin(),
         stats::rnbinom(200, mu = exp(linear), size = 2),
         list(augment = "negbin", gate = "hard")),
    list("lognormal AFT augmented", lognormal_aft(),
         cbind(exp(pmin(linear + stats::rnorm(200, sd = 0.5), 3)),
               as.numeric(linear + stats::rnorm(200, sd = 0.5) <= 3)),
         list(augment = "aft")),
    list("loglogistic AFT augmented", loglogistic_aft(),
         cbind(exp(pmin(linear + 0.5 * stats::rlogis(200), 3)),
               as.numeric(linear + 0.5 * stats::rlogis(200) <= 3)),
         list(augment = "aft")))

  for (s in settings) {
    dd <- d
    dd$y <- s[[3L]]

    control <- function(generic) {
      args <- c(list(num_trees = 10L, num_burn = 0L, num_draws = 25L,
                     verbose = FALSE, generic_accumulate = generic), s[[4L]])
      do.call(bartisan_control, args)
    }

    set.seed(5)
    virtual <- bartisan(y ~ ., data = dd, family = s[[2L]],
                        control = control(TRUE))
    set.seed(5)
    static <- bartisan(y ~ ., data = dd, family = s[[2L]],
                       control = control(FALSE))

    expect_identical(virtual[["eta"]], static[["eta"]], label = s[[1L]])
    expect_identical(virtual[["loglik"]], static[["loglik"]], label = s[[1L]])
  }
})

test_that("the proportional hazards derivatives are the ones it claims", {
  # The score and information of delta * eta - Lambda_0(y) exp(eta) are
  # delta - mu and mu with mu = Lambda_0(y) exp(eta), so the information is the
  # observed one and is positive whenever any exposure is accrued.
  set.seed(118)
  n <- 30
  time <- stats::rexp(n, 1)
  event <- stats::rbinom(n, 1L, 0.7)
  edges <- c(0, unname(stats::quantile(time, c(0.25, 0.5, 0.75))))
  lambda <- c(0.8, 1.2, 0.6, 1.5)
  eta <- matrix(stats::rnorm(2 * n, 0, 0.6), nrow = 2L)

  opts <- list(event = event, edges = edges, lambda_shape = 1,
               lambda_rate = 1, update_lambda = FALSE)
  aux <- matrix(c(lambda, 1), nrow = 2L, ncol = 5L, byrow = TRUE)

  got <- .bartisan_derivs(time, rep(1, n), list(eta), "ph", "log", opts, aux,
                          0L, FALSE, FALSE)

  span <- c(diff(edges), Inf)
  cum <- vapply(time, function(u) sum(lambda * pmin(pmax(u - edges, 0), span)),
                numeric(1))
  cum_wide <- matrix(cum, nrow = 2L, ncol = n, byrow = TRUE)
  event_wide <- matrix(event, nrow = 2L, ncol = n, byrow = TRUE)
  mu <- cum_wide * exp(eta)

  expect_equal(got[["d1"]], event_wide - mu, tolerance = 1e-12)
  expect_equal(got[["info"]], mu, tolerance = 1e-12)
  expect_true(all(got[["info"]] > 0))
})

test_that("the beta families' tabulated derivatives match the exact ones", {
  # Both digamma combinations the beta score and information need are functions
  # of the single scalar mu once phi is fixed, which it is for the whole of a
  # sweep, so they are tabulated on a grid in the additive predictor and
  # interpolated. That is a proposal, not the target, so the standard it has to
  # meet is closeness rather than exactness -- but it does have to be close.
  set.seed(117)
  n <- 40
  y <- stats::rbeta(n, 3, 5)
  phi <- 8
  eta <- matrix(stats::rnorm(2 * n, 0, 2.5), nrow = 2L)
  opts <- list(phi = phi, phi_prior_shape = 0.01, phi_prior_rate = 0.01,
               update_phi = FALSE)

  got <- .bartisan_derivs(y, rep(1, n), list(eta), "beta", "logit", opts,
                          matrix(phi, 2L, 1L), 0L, FALSE, FALSE)

  mu <- stats::plogis(eta)
  a <- mu * phi
  b <- phi - a
  slope <- phi * mu * (1 - mu)
  y_wide <- matrix(y, nrow = 2L, ncol = n, byrow = TRUE)

  exact_score <- slope * (log(y_wide) - log1p(-y_wide) - digamma(a) + digamma(b))
  exact_info <- slope^2 * (trigamma(a) + trigamma(b))

  # Interpolation error on a grid of 2049 points over (-8, 8), scaled by the
  # slope, lands several orders of magnitude below the quantities themselves.
  expect_equal(got[["d1"]], exact_score, tolerance = 1e-5)
  expect_equal(got[["info"]], exact_info, tolerance = 1e-5)
  expect_true(all(got[["info"]] > 0))

  # Outside the grid the exact functions are used, so agreement there is to the
  # last bits rather than to the interpolation tolerance.
  far <- matrix(c(-14, 14), nrow = 2L, ncol = n)
  got_far <- .bartisan_derivs(y, rep(1, n), list(far), "beta", "logit", opts,
                              matrix(phi, 2L, 1L), 0L, FALSE, FALSE)
  mu_f <- stats::plogis(far)
  a_f <- mu_f * phi
  b_f <- phi - a_f
  slope_f <- phi * mu_f * (1 - mu_f)
  expect_equal(got_far[["info"]],
               slope_f^2 * (trigamma(a_f) + trigamma(b_f)), tolerance = 1e-12)
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

  got <- .bartisan_derivs(y, rep(1, n), list(eta), "Gamma", "log", opts,
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
    bartisan_control(num_trees = 50, num_burn = 1000, num_draws = 1000,
                     gate = "hard", augment = augment, verbose = FALSE)
  }

  set.seed(5)
  direct <- bartisan(y ~ ., data = d, family = negbin(), control = ctrl(FALSE))
  set.seed(5)
  rewritten <- bartisan(y ~ ., data = d, family = negbin(),
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
