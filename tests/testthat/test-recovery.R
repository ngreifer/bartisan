# Statistical behavior, as opposed to plumbing. These take longer than the
# other files, so they are skipped on CRAN.

test_that("a flat likelihood reproduces the tree prior", {
  skip_on_cran()

  # The sharpest available check on the reversible-jump moves. Shrinking the
  # prior weights to nothing makes the likelihood constant in the leaf values,
  # so the target collapses to the tree prior. If birth, death and change
  # satisfy detailed balance, the sampled trees must reproduce the prior's
  # branching process; if any acceptance ratio is wrong, the tree size drifts.
  set.seed(101)

  n <- 120
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- stats::rnorm(n)

  gamma <- 0.95
  beta <- 2
  num_trees <- 30

  fit <- bartisan(y ~ ., data = d, family = gaussian(),
                  weights = rep(1e-10, n),
                 control = bartisan_control(num_trees = num_trees,
                                            num_burn = 1000, num_draws = 2000,
                                           gate = "hard", sigma_mu = 0.4,
                                           update_sigma_mu = FALSE,
                                           update_s = FALSE,
                                           update_alpha = FALSE,
                                           sigma_mu_ramp = 0, verbose = FALSE))

  # Expected leaves per tree, by backward recursion over depth.
  expected <- 1
  for (depth in 60:0) {
    grow <- gamma * (1 + depth)^(-beta)
    expected <- (1 - grow) + grow * 2 * expected
  }

  # A binary tree with L leaves has L - 1 splitting rules.
  observed <- mean(rowSums(fit[["counts"]][["eta"]])) / num_trees + 1

  expect_equal(observed, expected, tolerance = 0.05)
})

test_that("a flat likelihood reproduces the leaf prior", {
  skip_on_cran()

  set.seed(102)

  n <- 120
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- stats::rnorm(n)

  num_trees <- 30
  sigma_mu <- 0.4

  fit <- bartisan(y ~ ., data = d, family = gaussian(),
                  weights = rep(1e-10, n),
                 control = bartisan_control(num_trees = num_trees,
                                            num_burn = 1000, num_draws = 2000,
                                           gate = "hard", sigma_mu = sigma_mu,
                                           update_sigma_mu = FALSE,
                                           update_s = FALSE,
                                           update_alpha = FALSE,
                                           sigma_mu_ramp = 0, verbose = FALSE))

  # Each leaf value is an independent prior draw, so the forest's contribution
  # at any point has standard deviation sqrt(num_trees) * sigma_mu.
  centered <- fit[["eta"]][["eta"]] - fit[["intercept"]][1L]

  expect_equal(stats::sd(as.vector(centered)), sqrt(num_trees) * sigma_mu,
               tolerance = 0.1)
})

test_that("nonlinear signal is recovered across families", {
  skip_on_cran()

  set.seed(103)

  n <- 300
  X <- matrix(stats::runif(n * 5), nrow = n)
  colnames(X) <- paste0("x", 1:5)
  truth <- 2 * sin(pi * X[, 1] * X[, 2]) + 2 * (X[, 3] - 0.5)^2
  d <- as.data.frame(X)

  ctrl <- bartisan_control(num_trees = 20, num_burn = 400, num_draws = 400,
                           verbose = FALSE)

  d$y <- truth + stats::rnorm(n, sd = 0.3)
  fit <- bartisan(y ~ ., data = d, control = ctrl)
  expect_gt(stats::cor(predict(fit), truth), 0.9)

  d$y <- stats::rpois(n, exp(truth))
  fit <- bartisan(y ~ ., data = d, family = poisson(), control = ctrl)
  expect_gt(stats::cor(predict(fit), truth), 0.8)

  d$y <- stats::rgamma(n, 4, rate = 4 / exp(truth))
  fit <- bartisan(y ~ ., data = d, family = stats::Gamma("log"), control = ctrl)
  expect_gt(stats::cor(predict(fit), truth), 0.8)
  # The shape is the inverse dispersion and should land near its true value.
  expect_gt(mean(fit[["aux"]][, "shape"]), 2)
  expect_lt(mean(fit[["aux"]][, "shape"]), 8)
})

test_that("the residual scale is recovered for a Gaussian response", {
  skip_on_cran()

  set.seed(104)

  n <- 400
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- 2 * d$x1 + stats::rnorm(n, sd = 0.5)

  fit <- bartisan(y ~ ., data = d, family = gaussian(),
                  control = bartisan_control(num_trees = 20, num_burn = 400,
                                             num_draws = 400, verbose = FALSE))

  expect_equal(mean(fit[["aux"]][, "sigma"]), 0.5, tolerance = 0.15)
})

test_that("ordinal cutpoint spacing is recovered", {
  skip_on_cran()

  set.seed(105)

  n <- 500
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  eta <- 2 * d$x1 - 1
  cuts <- c(-0.5, 1, 2.5)

  probs <- vapply(seq_len(n), function(i) {
    diff(c(0, stats::plogis(cuts - eta[i]), 1))
  }, numeric(4))

  d$y <- factor(apply(probs, 2L, function(p) sample.int(4, 1, prob = p)),
                levels = 1:4, ordered = TRUE)

  fit <- bartisan(y ~ ., data = d, family = ordinal(),
                  control = bartisan_control(num_trees = 20, num_burn = 500,
                                             num_draws = 500, verbose = FALSE))

  # Only the gaps between cutpoints are identified; the level is absorbed by
  # the intercept, which is why the first cutpoint is pinned at zero.
  gaps <- colMeans(fit[["aux"]])[-1L]

  expect_equal(unname(gaps), cuts[-1L] - cuts[1L], tolerance = 0.4)
})

test_that("the beta precision and mean are recovered", {
  skip_on_cran()

  set.seed(108)

  n <- 600
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  mu <- stats::plogis(-0.3 + 1.5 * sin(pi * d$x1))
  phi <- 12
  d$y <- stats::rbeta(n, mu * phi, phi - mu * phi)

  fit <- bartisan(y ~ ., data = d, family = Beta(),
                  control = bartisan_control(num_trees = 50, num_burn = 500,
                                             num_draws = 500, verbose = FALSE))

  expect_equal(mean(fit[["aux"]][, "phi"]), phi, tolerance = 0.2)

  # The mean is on the logit scale, so the response scale is what to check it on.
  expect_lt(sqrt(mean((drop(predict(fit, type = "response")) - mu)^2)), 0.05)

  # The reported log likelihood is the beta one, rebuilt from the stored pieces.
  eta <- predict(fit, type = "link", draws = TRUE)
  drawn <- fit[["aux"]][, "phi"]
  by_hand <- vapply(seq_len(nrow(eta)), function(s) {
    m <- stats::plogis(eta[s, ])
    sum(stats::dbeta(d$y, m * drawn[s], drawn[s] - m * drawn[s], log = TRUE))
  }, numeric(1))

  expect_equal(fit[["loglik"]], by_hand, tolerance = 1e-8)
})

test_that("proportional hazards recovers the log hazard ratio and the baseline", {
  skip_on_cran()

  set.seed(109)

  # A baseline hazard that rises and then falls, which is what the family is for:
  # no accelerated failure time family here can represent it, since the Weibull's
  # hazard is monotone and the other two are unimodal on the log-time scale.
  n <- 800
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  x3 = stats::runif(n))
  truth <- 1.4 * sin(pi * d$x1) + 0.7 * (d$x2 - 0.5)
  truth <- truth - mean(truth)

  base <- function(t) 0.3 + 2 * exp(-((t - 1) / 0.35)^2)
  grid <- seq(0, 15, length.out = 3001)
  cum <- c(0, cumsum((base(grid[-1]) + base(grid[-length(grid)])) / 2 *
                       diff(grid)))
  event_time <- stats::approx(cum, grid,
                              -log(stats::runif(n)) / exp(truth), rule = 2)$y
  cens <- stats::runif(n, 0, stats::quantile(event_time, 0.85) * 2)

  d$time <- pmin(event_time, cens)
  d$status <- as.numeric(event_time <= cens)

  # The predictors are named rather than taken with `.`, which would hand the
  # model the response's own columns.
  fit <- bartisan(cbind(time, status) ~ x1 + x2 + x3, data = d, family = ph(),
                  control = bartisan_control(num_burn = 500, num_draws = 500,
                                             verbose = FALSE))

  # One hazard per bin, plus the rate its own prior is given.
  lambda_cols <- grep("^lambda[0-9]+$", colnames(fit[["aux"]]))
  expect_gt(length(lambda_cols), 2L)
  expect_true("lambda_rate" %in% colnames(fit[["aux"]]))
  expect_true(all(fit[["aux"]] > 0))

  # The predictor is a log hazard ratio, identified only against the baseline, so
  # it is compared after centering.
  r_hat <- drop(predict(fit, type = "link"))
  expect_gt(stats::cor(r_hat, truth), 0.9)
  expect_lt(sqrt(mean((r_hat - mean(r_hat) - truth)^2)), 0.5 * stats::sd(truth))

  # The fitted baseline has to turn over, which is the whole point. The last bin
  # runs to infinity and holds the fewest observations, so it is left out of this:
  # its hazard is the noisiest of them.
  lambda <- colMeans(fit[["aux"]])[lambda_cols]
  interior <- lambda[-length(lambda)]
  expect_gt(which.max(interior), 1L)
  expect_gt(max(interior) / interior[1L], 1.5)

  # And the reported log likelihood is the piecewise-exponential one.
  eta <- predict(fit, type = "link", draws = TRUE)
  drawn <- fit[["aux"]][, lambda_cols, drop = FALSE]
  edges <- fit[["family_opts"]][["edges"]]
  span <- c(diff(edges), Inf)
  bin <- findInterval(d$time, edges)

  by_hand <- vapply(seq_len(nrow(eta)), function(s) {
    lam <- drawn[s, ]
    cumulative <- vapply(d$time, function(u) {
      sum(lam * pmin(pmax(u - edges, 0), span))
    }, numeric(1))
    sum(d$status * (log(lam[bin]) + eta[s, ]) - cumulative * exp(eta[s, ]))
  }, numeric(1))

  expect_equal(fit[["loglik"]], by_hand, tolerance = 1e-8)

  # `type = "response"` is the median survival time, as it is for the
  # accelerated failure time families.
  med <- drop(predict(fit, type = "response"))
  expect_true(all(med > 0))

  # A rank correlation, because the median is a monotone but strongly nonlinear
  # function of the predictor: it inverts a cumulative baseline that turns over.
  expect_gt(stats::cor(med, -truth, method = "spearman"), 0.95)
})

test_that("the proportional hazards estimates do not depend on the bin count", {
  skip_on_cran()

  # `num_bins` is exposed for checking this rather than for tuning, so the check
  # belongs in the suite: over a wide range of bin counts the estimates should
  # move less than the replicate-to-replicate noise, while the effective number
  # of parameters grows with the count -- which is why the default is where it is.
  set.seed(112)

  n <- 500
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  truth <- 1.4 * sin(pi * d$x1)
  truth <- truth - mean(truth)

  base <- function(t) 0.3 + 2 * exp(-((t - 1) / 0.35)^2)
  grid <- seq(0, 20, length.out = 3001)
  cum <- c(0, cumsum((base(grid[-1]) + base(grid[-length(grid)])) / 2 *
                       diff(grid)))
  event_time <- stats::approx(cum, grid,
                              -log(stats::runif(n)) / exp(truth), rule = 2)$y
  cens <- stats::runif(n, 0, stats::quantile(event_time, 0.85) * 2)

  d$time <- pmin(event_time, cens)
  d$status <- as.numeric(event_time <= cens)

  fit_with <- function(bins) {
    bartisan(cbind(time, status) ~ x1 + x2, data = d,
             family = ph(num_bins = bins),
            control = bartisan_control(num_burn = 300, num_draws = 300,
                                       verbose = FALSE))
  }

  errors <- vapply(c(5L, 10L, 40L), function(bins) {
    r_hat <- drop(predict(fit_with(bins), type = "link"))
    sqrt(mean((r_hat - mean(r_hat) - truth)^2))
  }, numeric(1))

  # An eightfold range of bin counts, and the error barely moves.
  expect_lt(max(errors) - min(errors), 0.4 * stats::sd(truth))
  expect_true(all(errors < 0.6 * stats::sd(truth)))
})

test_that("the survival scale is recovered under censoring", {
  skip_on_cran()

  set.seed(106)

  n <- 500
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  log_t <- 2 * d$x1 + 0.7 * stats::rlogis(n)
  cens <- stats::quantile(log_t, 0.75)

  d$time <- exp(pmin(log_t, cens))
  d$status <- as.numeric(log_t <= cens)

  fit <- bartisan(survival::Surv(time, status) ~ x1 + x2, data = d,
                  family = loglogistic_aft(),
                 control = bartisan_control(num_trees = 20, num_burn = 500,
                                            num_draws = 500, verbose = FALSE))

  expect_equal(mean(fit[["aux"]][, "sigma"]), 0.7, tolerance = 0.2)
})

test_that("the augmented survival models target the same posterior", {
  skip_on_cran()
  skip_if_not_installed("survival")

  # Right-censoring is what makes the direct likelihood awkward: a failure
  # contributes a density and a censored observation a survival function, and the
  # two have different shapes in the predictor. Imputing the failure time above
  # its censoring time makes every contribution the same shape, and quadratic.
  # Heavy censoring is where the imputation carries the most weight, so that is
  # where the two samplers are compared.
  set.seed(107)

  n <- 400
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  x3 = stats::runif(n))
  truth <- 1 + 1.2 * sin(pi * d$x1)

  for (family in list(lognormal_aft(), loglogistic_aft())) {
    error <- if (identical(family[["link"]], "lognormal")) {
      stats::rnorm(n)
    } else {
      stats::rlogis(n)
    }

    time <- exp(truth + 0.5 * error)
    # Half the sample censored.
    cut <- stats::quantile(time, 0.5)
    dd <- d
    dd$time <- pmin(time, cut)
    dd$status <- as.numeric(time <= cut)

    ctrl <- function(augment) {
      bartisan_control(num_trees = 50, num_burn = 1000, num_draws = 1500,
                       augment = augment, verbose = FALSE)
    }

    set.seed(5)
    direct <- bartisan(survival::Surv(time, status) ~ x1 + x2 + x3, dd,
                       family = family, control = ctrl(FALSE))
    set.seed(5)
    imputed <- bartisan(survival::Surv(time, status) ~ x1 + x2 + x3, dd,
                        family = family, control = ctrl("aft"))

    label <- family[["link"]]

    # The scale is the global parameter and the one most likely to shift if the
    # imputation were wrong. The tolerance is a few times its posterior standard
    # deviation, which is about 0.025 here.
    expect_equal(mean(imputed[["aux"]][, "sigma"]),
                 mean(direct[["aux"]][, "sigma"]),
                 tolerance = 0.05, label = label)
    expect_equal(stats::sd(imputed[["aux"]][, "sigma"]),
                 stats::sd(direct[["aux"]][, "sigma"]),
                 tolerance = 0.3, label = label)

    # And the predictor, read against its own posterior spread rather than an
    # absolute scale: the two chains are different samplers, so they agree only
    # to Monte Carlo error, and measured in posterior standard deviations that
    # error is the same size whichever pair of chains is compared.
    a <- predict(imputed, type = "link", draws = TRUE)
    b <- predict(direct, type = "link", draws = TRUE)
    spread <- apply(b, 2L, stats::sd)
    expect_lt(mean(abs(colMeans(a) - colMeans(b)) / spread), 0.4)

    # The reported log likelihood is still the observed-data one, censoring and
    # all, not the complete-data likelihood the sampler works with.
    sigma <- imputed[["aux"]][, "sigma"]
    by_hand <- vapply(seq_len(nrow(a)), function(s) {
      r <- (log(dd$time) - a[s, ]) / sigma[s]
      log_dens <- if (identical(label, "lognormal")) {
        stats::dnorm(r, log = TRUE) - log(sigma[s])
      } else {
        stats::dlogis(r, log = TRUE) - log(sigma[s])
      }
      log_surv <- if (identical(label, "lognormal")) {
        stats::pnorm(r, lower.tail = FALSE, log.p = TRUE)
      } else {
        stats::plogis(r, lower.tail = FALSE, log.p = TRUE)
      }
      sum(ifelse(dd$status > 0, log_dens, log_surv))
    }, numeric(1))

    expect_equal(imputed[["loglik"]], by_hand, tolerance = 1e-8, label = label)
  }
})

test_that("irrelevant predictors are used less than relevant ones", {
  skip_on_cran()

  set.seed(107)

  n <- 300
  X <- matrix(stats::runif(n * 8), nrow = n)
  colnames(X) <- paste0("x", 1:8)
  d <- as.data.frame(X)
  d$y <- 3 * X[, 1] + 3 * X[, 2] + stats::rnorm(n, sd = 0.3)

  fit <- bartisan(y ~ ., data = d,
                  control = bartisan_control(num_trees = 20, num_burn = 500,
                                             num_draws = 500, verbose = FALSE))

  used <- colMeans(fit[["counts"]][["eta"]])

  expect_gt(min(used[c("x1", "x2")]), max(used[paste0("x", 3:8)]))
})

test_that("a flat likelihood reproduces the whole prior distribution of tree sizes", {
  skip_on_cran()

  # Stronger than matching the mean tree size: the sampler must reproduce the
  # entire distribution of the number of leaves. A wrong acceptance ratio can
  # leave the mean intact while distorting the shape.
  #
  # The reference is exact, not simulated. Writing P_d for the distribution of
  # the number of leaves in a subtree rooted at depth d, and rho_d for the
  # branching probability there,
  #
  #   P_d(1) = 1 - rho_d,   P_d(L) = rho_d * sum_j P_{d+1}(j) P_{d+1}(L - j)
  #
  # since a branch's two subtrees are independent and their leaves add.
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

  set.seed(11)

  n <- 120
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  x3 = stats::runif(n))
  d$y <- stats::rnorm(n)

  # A single tree, so the recorded split count is that tree's alone: a binary
  # tree with L leaves has L - 1 splits.
  fit <- bartisan(y ~ ., data = d, family = gaussian(),
                  weights = rep(1e-10, n),
                 control = bartisan_control(num_trees = 1, num_burn = 2000,
                                            num_draws = 20000, gate = "hard",
                                           sigma_mu = 0.4,
                                           update_sigma_mu = FALSE,
                                           update_s = FALSE,
                                           update_alpha = FALSE,
                                           sigma_mu_ramp = 0, verbose = FALSE))
  sampled <- rowSums(fit[["counts"]][["eta"]]) + 1L

  cells <- 1:6
  reference <- leaf_pmf()
  expected <- c(reference[cells], 1 - sum(reference[cells]))
  observed <- c(vapply(cells, function(j) sum(sampled == j), numeric(1)),
                sum(sampled > max(cells)))

  stat <- suppressWarnings(
    stats::chisq.test(observed, p = expected / sum(expected)))

  # Successive draws are strongly dependent, so the nominal chi-square is
  # anti-conservative; scale the statistic to the effective sample size.
  rho <- stats::acf(sampled, plot = FALSE, lag.max = 200)$acf[-1]
  cut <- which(rho < 0.02)[1]

  if (is.na(cut)) {
    cut <- 200L
  }

  ess <- length(sampled) / (1 + 2 * sum(rho[1:cut]))
  scaled <- as.numeric(stat[["statistic"]]) * ess / length(sampled)

  expect_gt(stats::pchisq(scaled, stat[["parameter"]], lower.tail = FALSE), 0.01)

  # And the mean, which is the sharpest single number.
  expect_equal(mean(sampled), sum(reference * seq_along(reference)),
               tolerance = 0.02)
})

test_that("the two multinomial codings put the same prior on the log odds", {
  skip_on_cran()

  # Shrinking the weights to nothing makes the likelihood constant, so what the
  # sampler draws is the prior. The symmetric coding writes each log-odds
  # contrast as a difference of two forests, doubling its variance, which is why
  # the leaf scale is divided by sqrt(2); the two should then agree.
  set.seed(41)

  n <- 150
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- factor(sample(c("a", "b", "c"), n, replace = TRUE))

  flat <- function(family) {
    bartisan(y ~ ., data = d, family = family, weights = rep(1e-10, n),
             control = bartisan_control(num_trees = 20, num_burn = 300,
                                        num_draws = 1500, gate = "hard",
                                      update_sigma_mu = FALSE,
                                      update_s = FALSE, update_alpha = FALSE,
                                      sigma_mu_ramp = 0, verbose = FALSE))
  }

  # Contrasts of b and c against a, under each coding.
  sym <- flat(multinomial())
  ref <- flat(multinomial(reference = "a"))

  sd_sym <- c(stats::sd(sym[["eta"]][["b"]] - sym[["eta"]][["a"]]),
              stats::sd(sym[["eta"]][["c"]] - sym[["eta"]][["a"]]))
  sd_ref <- c(stats::sd(ref[["eta"]][["b"]]), stats::sd(ref[["eta"]][["c"]]))

  expect_equal(mean(sd_sym), mean(sd_ref), tolerance = 0.1)

  # The unidentified direction -- adding the same function to every category --
  # is pinned only by the prior. It must stay put rather than drift, which is
  # what makes the posterior proper.
  level <- (sym[["eta"]][["a"]] + sym[["eta"]][["b"]] +
              sym[["eta"]][["c"]]) / 3
  by_half <- rowMeans(level)
  half <- seq_len(length(by_half) / 2)

  expect_lt(abs(mean(by_half[half]) - mean(by_half[-half])),
            3 * stats::sd(by_half))
  expect_lt(max(abs(by_half)), 10 * stats::sd(by_half) + 1)
})

test_that("dpm_aft recovers a bimodal error distribution under censoring", {
  skip_on_cran()

  # Henderson, Louis, Rosner and Varadhan (2020): log T = m(x) + W with W a
  # mean-constrained Dirichlet process mixture, and censored log-times imputed.
  # A two-component error is the case no other survival family here can
  # represent, so it is the one to check.
  set.seed(112)

  n <- 900
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  x3 = stats::runif(n))
  truth <- 1 + 1.2 * sin(pi * d$x1) + 0.6 * (d$x2 - 0.5)

  heavy <- stats::rbinom(n, 1L, 0.35)
  error <- ifelse(heavy == 1, stats::rnorm(n, 1.1, 0.25),
                  stats::rnorm(n, -0.6, 0.35))
  error <- error - (0.35 * 1.1 + 0.65 * -0.6)

  event_time <- exp(truth + error)
  cens <- stats::runif(n, 0, stats::quantile(event_time, 0.8) * 2.5)
  d$time <- pmin(event_time, cens)
  d$status <- as.numeric(event_time <= cens)

  fit <- bartisan(cbind(time, status) ~ x1 + x2 + x3, d, family = dpm_aft(),
                  control = bartisan_control(num_burn = 500, num_draws = 500,
                                             verbose = FALSE))

  # The predictor is the conditional mean of log T, as it is for `dpm()`, so it
  # is comparable against the truth on its own level rather than only centered.
  fitted <- drop(predict(fit, type = "link"))
  expect_gt(stats::cor(fitted, truth), 0.95)
  expect_lt(sqrt(mean((fitted - truth)^2)), 0.3 * stats::sd(truth))

  # The mixture found more than one component, which is the whole point.
  expect_gt(mean(fit[["aux"]][, "clusters"]), 1.5)
  expect_identical(colnames(fit[["aux"]]),
                   c("alpha", "clusters", "center", "error_sd"))

  # And the fitted error density has both modes.
  at <- seq(-2, 2, by = 0.1)
  density <- error_density(fit, at = at)[["mean"]]
  peaks <- at[which(diff(sign(diff(density))) < 0) + 1L]
  expect_gte(length(peaks), 2L)

  # The reported log likelihood is the *observed-data* one -- a density for an
  # event, a survival probability for a censoring -- and the density route
  # reproduces it draw by draw. That is the check that the two independent
  # implementations of it agree, censoring included.
  by_row <- rowSums(predict(fit, type = "density", draws = TRUE, log = TRUE))
  expect_equal(by_row, fit[["loglik"]], tolerance = 1e-8)

  # Prior weights are refused, as they are for `dpm()`.
  d$w <- stats::runif(n, 0.5, 1.5)
  expect_error(bartisan(cbind(time, status) ~ x1 + x2 + x3, d, family = dpm_aft(),
                       weights = w, control = quick_control()),
               "does not take prior weights")
})

test_that("dpm_aft costs nothing when a single normal error is right", {
  skip_on_cran()

  # The property that makes it a sensible default rather than a specialist tool:
  # on log-normal errors, where `lognormal_aft()` is correctly specified, the
  # mixture should be level with it rather than paying for its flexibility.
  set.seed(113)

  n <- 700
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  truth <- 1 + 1.2 * sin(pi * d$x1)
  event_time <- exp(truth + stats::rnorm(n, 0, 0.5))
  cens <- stats::runif(n, 0, stats::quantile(event_time, 0.8) * 2.5)
  d$time <- pmin(event_time, cens)
  d$status <- as.numeric(event_time <= cens)

  held <- d[1:200, ]
  train <- d[201:n, ]
  ctrl <- bartisan_control(num_burn = 400, num_draws = 400, verbose = FALSE)

  score <- function(family) {
    fit <- bartisan(cbind(time, status) ~ x1 + x2, train, family = family,
                    control = ctrl)
    sum(log(predict(fit, newdata = held, type = "density")))
  }

  mixture <- score(dpm_aft())
  single <- score(lognormal_aft())

  # Both densities are with respect to log T, so the two scores are comparable.
  # Within a few log points of each other on 200 held-out observations.
  expect_lt(abs(mixture - single), 15)
})
