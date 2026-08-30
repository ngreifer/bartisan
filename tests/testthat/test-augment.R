# Rewriting the probit likelihood as the margin of a Gaussian one changes the
# sampler, not the target. These check that the target really is unchanged, that
# the reported log likelihood is still the probit one rather than the augmented
# density, and that the rewriting declines where it does not apply.

test_that("the Polya-Gamma sampler has the moments it should", {
  # E[PG(b, c)] = (b / 2c) tanh(c / 2) and Var[PG(b, c)] = (b / 4c^3)
  # (sinh c - c) sech^2(c / 2), with the c = 0 limits b/4 and b/24. Integer b
  # goes through Devroye's exact method and non-integer b through the series, so
  # both routes are checked.
  set.seed(91)

  pg_mean <- function(b, c) {
    if (abs(c) < 1e-8) b / 4 else b * tanh(c / 2) / (2 * c)
  }

  pg_var <- function(b, c) {
    if (abs(c) < 1e-8) return(b / 24)
    b * (sinh(c) - c) / (4 * c^3) / cosh(c / 2)^2
  }

  for (b in c(1, 3, 0.5, 7.3)) {
    for (cc in c(0, 1.5, 8)) {
      x <- .bartisan_rpg(20000L, b, cc)
      label <- paste("b =", b, "c =", cc)

      expect_true(all(x > 0), info = label)
      # Four standard errors of the mean, which is a two-sided level of 6e-5 per
      # cell and so about 7e-4 over the twelve of them.
      expect_lt(abs(mean(x) - pg_mean(b, cc)) /
                  (stats::sd(x) / sqrt(length(x))), 4)
      expect_equal(stats::var(x), pg_var(b, cc), tolerance = 0.06,
                   info = label)
    }
  }
})

test_that("augment resolves to the families it says it does", {
  every <- c("binomial", "ordinal", "multinomial", "zip", "zinb", "aft")

  # The negative binomial is the one whose rewriting depends on the rules: its
  # gain is the exponential form, which only hard rules get.
  expect_setequal(bartisan_control(augment = TRUE)[["augment"]], every)
  expect_setequal(bartisan_control(augment = TRUE, gate = "hard")[["augment"]],
                  c(every, "negbin"))
  expect_identical(bartisan_control(augment = FALSE)[["augment"]],
                   character())
  expect_identical(bartisan_control(augment = "multinomial")[["augment"]],
                   "multinomial")
  expect_setequal(bartisan_control(augment = c("binomial", "negbin"))[["augment"]],
                  c("binomial", "negbin"))
  expect_setequal(bartisan_control(augment = c("zip", "zinb"))[["augment"]],
                  c("zip", "zinb"))
  expect_error(bartisan_control(augment = "poisson"), "must be one of")
})

test_that("the reported log likelihood is the probit one, not the augmented one", {
  d <- sim_x(n = 120, seed = 71)
  d$y <- stats::rbinom(120, 1, stats::plogis(2 * d$x1 - 1))

  fit <- bartisan(y ~ ., data = d, family = binomial("probit"),
                  control = quick_control(augment = TRUE))

  # Sum over observations of log Phi(eta) or log Phi(-eta), which is the probit
  # log likelihood and has nothing to do with the latent normal the sampler
  # actually works with.
  eta <- predict(fit, type = "link", draws = TRUE)
  by_hand <- apply(eta, 1L, function(e) {
    sum(ifelse(d$y == 1, stats::pnorm(e, log.p = TRUE),
               stats::pnorm(-e, log.p = TRUE)))
  })

  expect_equal(fit[["loglik"]], by_hand, tolerance = 1e-10)
})

test_that("augmentation declines where it does not apply", {
  d <- sim_x(n = 80, seed = 72)

  # A binomial count response needs one latent per trial, so there is no single
  # latent to draw and the direct family is used. The check is that the reported
  # log likelihood still carries the binomial coefficients, which the augmented
  # family knows nothing about.
  trials <- rep(4, 80)
  d$y <- stats::rbinom(80, 4, 0.4) / 4

  plain <- bartisan(y ~ ., data = d, family = binomial("probit"),
                    weights = trials, control = quick_control())
  asked <- bartisan(y ~ ., data = d, family = binomial("probit"),
                    weights = trials, control = quick_control(augment = TRUE))

  # Both fits used the same family, so the reported likelihoods live on the same
  # scale; the augmented one would be missing the binomial coefficients and sit
  # about 80 * log(6) higher.
  expect_lt(abs(mean(plain[["loglik"]]) - mean(asked[["loglik"]])),
            0.2 * abs(mean(plain[["loglik"]])))

  # A logit link is not a probit one, so nothing is rewritten.
  d$y <- stats::rbinom(80, 1, 0.4)
  a <- bartisan(y ~ ., data = d, family = binomial("logit"),
                control = quick_control(augment = TRUE))
  expect_true(all(is.finite(a[["loglik"]])))
})

test_that("augmentation targets the same posterior as the direct sampler", {
  skip_on_cran()

  set.seed(73)
  n <- 500
  X <- matrix(stats::runif(n * 4), n)
  colnames(X) <- paste0("x", 1:4)
  d <- as.data.frame(X)
  truth <- 1.5 * sin(pi * X[, 1]) + 1.5 * (X[, 2] - 0.5) - 0.3
  d$y <- stats::rbinom(n, 1, stats::pnorm(truth))

  ctrl <- function(augment) {
    bartisan_control(num_trees = 50, num_burn = 800, num_save = 800,
                     augment = augment, verbose = FALSE)
  }

  set.seed(5)
  direct <- bartisan(y ~ ., data = d, family = binomial("probit"),
                     control = ctrl(FALSE))
  set.seed(5)
  augmented <- bartisan(y ~ ., data = d, family = binomial("probit"),
                        control = ctrl(TRUE))

  p_direct <- predict(direct, type = "response")
  p_augmented <- predict(augmented, type = "response")
  p_true <- stats::pnorm(truth)

  # Two chains at the same target, so they agree to Monte Carlo error rather
  # than exactly, and both sit the same distance from the truth.
  expect_equal(sqrt(mean((p_augmented - p_true)^2)),
               sqrt(mean((p_direct - p_true)^2)), tolerance = 0.15)
  expect_gt(stats::cor(p_direct, p_augmented), 0.97)

  # And the posterior spread is the same, which is the part a broken
  # augmentation would get wrong even while the mean looked fine.
  sd_direct <- apply(predict(direct, type = "response", draws = TRUE), 2L,
                     stats::sd)
  sd_augmented <- apply(predict(augmented, type = "response", draws = TRUE),
                        2L, stats::sd)
  expect_equal(mean(sd_augmented), mean(sd_direct), tolerance = 0.15)
})

test_that("the logit rewriting reports the binomial log likelihood", {
  d <- sim_x(n = 120, seed = 92)
  d$y <- stats::rbinom(120, 1, stats::plogis(2 * d$x1 - 1))

  fit <- bartisan(y ~ ., data = d, family = binomial("logit"),
                  control = quick_control(augment = TRUE))

  eta <- predict(fit, type = "link", draws = TRUE)
  by_hand <- apply(eta, 1L, function(e) {
    sum(d$y * stats::plogis(e, log.p = TRUE) +
          (1 - d$y) * stats::plogis(-e, log.p = TRUE))
  })

  expect_equal(fit[["loglik"]], by_hand, tolerance = 1e-9)
  expect_predictor_invariant(fit, d)
})

test_that("the logit rewriting handles binomial counts", {
  d <- sim_x(n = 100, seed = 93)
  trials <- sample(1:5, 100, replace = TRUE)
  successes <- stats::rbinom(100, trials, stats::plogis(2 * d$x1 - 1))
  d$y <- successes / trials

  fit <- bartisan(y ~ ., data = d, family = binomial("logit"), weights = trials,
                  control = quick_control(augment = TRUE))

  eta <- predict(fit, type = "link", draws = TRUE)
  by_hand <- apply(eta, 1L, function(e) {
    sum(stats::dbinom(successes, trials, stats::plogis(e), log = TRUE))
  })

  expect_equal(fit[["loglik"]], by_hand, tolerance = 1e-9)
  expect_predictor_invariant(fit, d)
})

test_that("the negative binomial rewriting reports its own log likelihood", {
  d <- sim_x(n = 120, seed = 94)
  d$y <- stats::rnbinom(120, mu = exp(1 + d$x1), size = 3)

  fit <- bartisan(y ~ ., data = d, family = negbin(theta = 3),
                  control = quick_control(augment = "negbin"))

  eta <- predict(fit, type = "link", draws = TRUE)
  by_hand <- apply(eta, 1L, function(e) {
    sum(stats::dnbinom(d$y, mu = exp(e), size = 3, log = TRUE))
  })

  expect_equal(fit[["loglik"]], by_hand, tolerance = 1e-8)
  expect_predictor_invariant(fit, d)

  # The dispersion is still drawn when it is not fixed.
  drawn <- bartisan(y ~ ., data = d, family = negbin(),
                    control = quick_control(augment = "negbin"))
  expect_identical(colnames(drawn[["aux"]]), "theta")
  expect_true(all(drawn[["aux"]][, "theta"] > 0))
})

test_that("the multinomial rewriting reports its own log likelihood", {
  d <- sim_x(n = 120, seed = 95)
  d$y <- factor(sample(c("a", "b", "c"), 120, replace = TRUE))

  for (reference in list(NULL, "a")) {
    fit <- bartisan(y ~ ., data = d,
                    family = multinomial(reference = reference),
                   control = quick_control(augment = "multinomial"))

    probs <- predict(fit, type = "prob", draws = TRUE)
    observed <- match(d$y, dimnames(probs)[[3L]])
    by_hand <- vapply(seq_len(dim(probs)[1L]), function(s) {
      one <- probs[s, , , drop = FALSE][1L, , ]
      sum(log(one[cbind(seq_along(observed), observed)]))
    }, numeric(1))

    expect_equal(fit[["loglik"]], by_hand, tolerance = 1e-8)
    expect_predictor_invariant(fit, d)
  }
})

test_that("a rewriting is declined where it cannot apply", {
  d <- sim_x(n = 80, seed = 96)

  # Both the negative binomial and the multinomial augmentations need one trial
  # per observation, so a fractional weight sends them back to the direct
  # family. The check is that the fit still runs and reports a sane likelihood.
  w <- stats::runif(80, 0.5, 1.5)
  d$y <- stats::rpois(80, 2)
  fit <- bartisan(y ~ ., data = d, family = negbin(), weights = w,
                  control = quick_control(augment = "negbin"))
  expect_true(all(is.finite(fit[["loglik"]])))

  # A cloglog link has no Polya-Gamma form, so asking for the binomial
  # rewriting leaves it alone.
  d$y <- stats::rbinom(80, 1, 0.4)
  fit <- bartisan(y ~ ., data = d, family = binomial("cloglog"),
                  control = quick_control(augment = TRUE))
  expect_true(all(is.finite(fit[["loglik"]])))
  expect_predictor_invariant(fit, d)
})

test_that("the logit rewriting targets the same posterior", {
  skip_on_cran()

  set.seed(97)
  n <- 500
  X <- matrix(stats::runif(n * 4), n)
  colnames(X) <- paste0("x", 1:4)
  d <- as.data.frame(X)
  truth <- 1.5 * sin(pi * X[, 1]) + 1.5 * (X[, 2] - 0.5) - 0.3
  d$y <- stats::rbinom(n, 1, stats::plogis(truth))

  ctrl <- function(a) {
    bartisan_control(num_trees = 50, num_burn = 1000, num_save = 1000,
                     augment = a, verbose = FALSE)
  }

  set.seed(5)
  direct <- bartisan(y ~ ., data = d, family = binomial("logit"),
                     control = ctrl(FALSE))
  set.seed(5)
  rewritten <- bartisan(y ~ ., data = d, family = binomial("logit"),
                        control = ctrl(TRUE))

  p_true <- stats::plogis(truth)
  p_direct <- predict(direct, type = "response")
  p_rewritten <- predict(rewritten, type = "response")

  expect_gt(stats::cor(p_direct, p_rewritten), 0.97)
  expect_equal(sqrt(mean((p_rewritten - p_true)^2)),
               sqrt(mean((p_direct - p_true)^2)), tolerance = 0.2)

  sd_direct <- mean(apply(predict(direct, type = "response", draws = TRUE), 2L,
                          stats::sd))
  sd_rewritten <- mean(apply(predict(rewritten, type = "response",
                                     draws = TRUE), 2L, stats::sd))
  expect_equal(sd_rewritten, sd_direct, tolerance = 0.2)
})

test_that("the truncated normal draw matches its exact moments, including in the tails", {
  skip_on_cran()
  set.seed(4)

  # Exact mean and sd of a standard normal restricted to (lo, hi), computed in
  # whichever tail is representable -- the naive formula divides by
  # pnorm(hi) - pnorm(lo), which is exactly zero once both round to one.
  moments <- function(lo, hi) {
    if (lo + hi > 0) {
      m <- moments(-hi, -lo)
      return(c(-m[1L], m[2L]))
    }

    z <- exp(stats::pnorm(hi, log.p = TRUE)) - exp(stats::pnorm(lo, log.p = TRUE))
    lower <- if (is.infinite(lo)) 0 else lo * stats::dnorm(lo)
    upper <- if (is.infinite(hi)) 0 else hi * stats::dnorm(hi)
    m <- (stats::dnorm(lo) - stats::dnorm(hi)) / z
    c(m, sqrt(1 + (lower - upper) / z - m^2))
  }

  bounds <- list(c(-1, 1), c(0.5, 2), c(3, 4), c(-4, -3), c(5, 6),
                 c(-8, -7.5), c(8, 9), c(10, 12), c(20, 21))

  for (b in bounds) {
    z <- .bartisan_rtruncnorm(40000L, b[1L], b[2L])
    expect_true(all(z >= b[1L] & z <= b[2L]))
    exact <- moments(b[1L], b[2L])
    # Three standard errors of the mean, and a loose relative check on the sd.
    expect_lt(abs(mean(z) - exact[1L]), 3 * exact[2L] / sqrt(40000))
    expect_lt(abs(stats::sd(z) / exact[2L] - 1), 0.05)
  }

  # One-sided and unbounded reduce correctly.
  expect_true(all(.bartisan_rtruncnorm(1000L, 0, Inf) >= 0))
  expect_true(all(.bartisan_rtruncnorm(1000L, -Inf, 0) <= 0))
  expect_lt(abs(mean(.bartisan_rtruncnorm(40000L, -Inf, Inf))), 3 / sqrt(40000))
})

test_that("the ordinal probit augmentation targets the same posterior as the direct fit", {
  skip_on_cran()

  set.seed(11)
  n <- 400
  x <- sim_x(n = n, p = 3, seed = 11)
  lin <- 2 * sin(pi * x$x1) - x$x2
  cuts <- stats::quantile(lin, c(1, 2) / 3)
  d <- x
  d$y <- ordered(rowSums(outer(lin + stats::rnorm(n), cuts, ">")) + 1L)

  fit <- function(augment) {
    set.seed(5)
    bartisan(y ~ ., d, family = ordinal("probit"), gate = "hard",
             num_trees = 20, num_burn = 400, num_save = 400, augment = augment)
  }

  direct <- fit(FALSE)
  aug <- fit("ordinal")

  # The two samplers explore the same target, so their posterior means agree to
  # Monte Carlo error. The predictor is only identified up to what the cutpoints
  # absorb, so it is compared through its correlation with the direct fit's.
  expect_gt(stats::cor(colMeans(direct[["eta"]][[1L]]),
                       colMeans(aug[["eta"]][[1L]])), 0.95)

  # The cutpoints are identified, and are the sharper comparison.
  expect_equal(colMeans(direct[["aux"]]), colMeans(aug[["aux"]]),
               tolerance = 0.15)

  # The reported log likelihood is the ordinal one either way, not the augmented
  # target, so the two are on the same scale.
  expect_lt(abs(mean(direct[["loglik"]]) - mean(aug[["loglik"]])),
            0.05 * abs(mean(direct[["loglik"]])))

  expect_predictor_invariant(aug, d)
})

test_that("the ordinal augmentation applies only where it is exact", {
  d <- sim_x(n = 150, seed = 12)
  set.seed(1012)
  # With noise, so the response is not separable by the predictors and the fits
  # below have no reason to warn about a runaway leaf scale.
  z <- d$x1 - d$x2 + stats::rnorm(nrow(d))
  d$y <- ordered(rowSums(outer(z, stats::quantile(z, c(1, 2) / 3), ">")) + 1L)

  # A logit link has no Gaussian margin, so it falls back on the direct family
  # rather than silently fitting a probit.
  logit <- bartisan(y ~ ., d, family = ordinal("logit"), control = quick_control(),
                    augment = "ordinal")
  probit <- bartisan(y ~ ., d, family = ordinal("probit"),
                     control = quick_control(), augment = "ordinal")

  expect_s3_class(logit, "bartisan_fit")
  expect_s3_class(probit, "bartisan_fit")
  expect_predictor_invariant(logit, d)
  expect_predictor_invariant(probit, d)

  # Prior weights other than one break the latent-normal representation, so they
  # fall back too.
  w <- rep(c(1, 2), length.out = nrow(d))
  weighted <- bartisan(y ~ ., d, family = ordinal("probit"), weights = w,
                       control = quick_control(), augment = "ordinal")
  expect_s3_class(weighted, "bartisan_fit")
  expect_predictor_invariant(weighted, d)

  # "ordinal" is in the default set, and is accepted by name.
  expect_true("ordinal" %in% bartisan_control(gate = "smoothstep")$augment)
  expect_true("ordinal" %in% bartisan_control(gate = "hard")$augment)
  expect_identical(bartisan_control(augment = "ordinal")$augment, "ordinal")
})

test_that("the ordinal augmentation recovers known cutpoints", {
  skip_on_cran()

  set.seed(21)
  n <- 800
  x <- sim_x(n = n, p = 3, seed = 21)
  lin <- 2 * x$x1 - x$x2
  truth <- c(0, 1.2)
  z <- lin + stats::rnorm(n)
  d <- x
  d$y <- ordered(rowSums(outer(z, truth, ">")) + 1L)

  fit <- bartisan(y ~ ., d, family = ordinal("probit"), gate = "hard",
                  num_trees = 20, num_burn = 400, num_save = 400,
                 augment = "ordinal")

  # The cutpoints are reported in the chart where the predictor has mean zero
  # over the sample, so the truth has to be moved into the same chart before it
  # can be compared -- the two differ by the mean of the generating predictor.
  expect_equal(colMeans(fit[["aux"]]), truth - mean(lin), tolerance = 0.25,
               ignore_attr = TRUE)

  # The gaps between cutpoints are what the chart does not touch.
  expect_equal(diff(colMeans(fit[["aux"]])), diff(truth), tolerance = 0.25,
               ignore_attr = TRUE)
})

test_that("the logistic density is the Polya-Gamma normal scale mixture the ordinal logit augmentation rests on", {
  skip_on_cran()
  set.seed(5)

  # Polson, Scott and Windle (2013) Theorem 1 at a = 1, b = 2 says
  #   e^x / (1 + e^x)^2  =  (1/4) E[exp(-w x^2 / 2)],  w ~ PG(2, 0),
  # and the left side is the standard logistic density. Everything the ordinal
  # logit augmentation does follows from this one identity, so it is checked
  # against the density directly rather than assumed.
  w <- .bartisan_rpg(2e6L, 2, 0)

  expect_equal(mean(w), 0.5, tolerance = 0.01)

  for (x in c(0, 0.25, 0.5, 1, 2, 3)) {
    expect_equal(mean(exp(-w * x^2 / 2)) / 4, stats::dlogis(x),
                 tolerance = 0.005, info = paste("x =", x))
  }

  # And the tilting property the conditional draw uses: PG(b, c) is PG(b, 0)
  # reweighted by exp(-c^2 w / 2), which is why the conditional of the precision
  # given a residual r is exactly PG(2, |r|).
  for (cc in c(0.5, 1.5, 3)) {
    tilted <- sum(w * exp(-cc^2 * w / 2)) / sum(exp(-cc^2 * w / 2))
    expect_equal(mean(.bartisan_rpg(2e5L, 2, cc)), tilted, tolerance = 0.02,
                 info = paste("c =", cc))
  }
})

test_that("the ordinal logit augmentation targets the same posterior as the direct fit", {
  skip_on_cran()

  set.seed(31)
  n <- 400
  x <- sim_x(n = n, p = 3, seed = 31)
  lin <- 2 * sin(pi * x$x1) - x$x2
  cuts <- stats::quantile(lin, c(1, 2) / 3)
  d <- x
  d$y <- ordered(rowSums(outer(lin + stats::rlogis(n), cuts, ">")) + 1L)

  fit <- function(augment) {
    set.seed(6)
    bartisan(y ~ ., d, family = ordinal("logit"), gate = "hard", num_trees = 20,
             num_burn = 400, num_save = 400, augment = augment)
  }

  direct <- fit(FALSE)
  aug <- fit("ordinal")

  expect_gt(stats::cor(colMeans(direct[["eta"]][[1L]]),
                       colMeans(aug[["eta"]][[1L]])), 0.95)

  # The cutpoints are the identified part and so the sharper comparison.
  expect_equal(colMeans(direct[["aux"]]), colMeans(aug[["aux"]]),
               tolerance = 0.2)

  # Reported on the ordinal scale either way, not on the augmented one.
  expect_lt(abs(mean(direct[["loglik"]]) - mean(aug[["loglik"]])),
            0.05 * abs(mean(direct[["loglik"]])))

  expect_predictor_invariant(aug, d)
})

test_that("both ordinal links are augmented, and only where it is exact", {
  d <- sim_x(n = 250, seed = 33)
  set.seed(1033)
  # Enough noise that the categories are not close to separable, so the fits
  # below have no reason to warn about a runaway leaf scale.
  z <- d$x1 - d$x2 + 2 * stats::rlogis(nrow(d))
  d$y <- ordered(rowSums(outer(z, stats::quantile(z, c(1, 2) / 3), ">")) + 1L)

  # Both links now have a Gaussian margin, so `augment = TRUE` covers both. The
  # leaf scale is held fixed because this is a plumbing check on a 30-draw run,
  # which is far too short for it to settle and would warn on any family.
  for (link in c("logit", "probit")) {
    fit <- bartisan(y ~ ., d, family = ordinal(link),
                    control = quick_control(update_sigma_mu = FALSE))
    expect_s3_class(fit, "bartisan_fit")
    expect_predictor_invariant(fit, d)
  }

  # Frequency weights break the one-latent-per-observation representation, so
  # they fall back on the direct family.
  w <- rep(c(1, 2), length.out = nrow(d))
  weighted <- bartisan(y ~ ., d, family = ordinal("logit"), weights = w,
                       control = quick_control(update_sigma_mu = FALSE),
                      augment = "ordinal")
  expect_s3_class(weighted, "bartisan_fit")
  expect_predictor_invariant(weighted, d)
})

test_that("the ordinal logit augmentation recovers known cutpoints", {
  skip_on_cran()

  set.seed(35)
  n <- 800
  x <- sim_x(n = n, p = 3, seed = 35)
  lin <- 2 * x$x1 - x$x2
  truth <- c(0, 1.5)
  d <- x
  d$y <- ordered(rowSums(outer(lin + stats::rlogis(n), truth, ">")) + 1L)

  fit <- bartisan(y ~ ., d, family = ordinal("logit"), gate = "hard",
                  num_trees = 20, num_burn = 400, num_save = 400,
                 augment = "ordinal")

  expect_equal(colMeans(fit[["aux"]]), truth - mean(lin), tolerance = 0.3,
               ignore_attr = TRUE)
  expect_equal(diff(colMeans(fit[["aux"]])), diff(truth), tolerance = 0.3,
               ignore_attr = TRUE)
})

test_that("the zero-inflated rewriting reports its own log likelihood", {
  d <- sim_x(n = 150, seed = 151)
  set.seed(1151)
  structural <- stats::runif(nrow(d)) < 0.3
  d$y <- ifelse(structural, 0L, stats::rpois(nrow(d), exp(0.7 + d$x1)))

  fit <- bartisan(y ~ ., data = d, family = zi_poisson(),
                  control = quick_control(augment = "zip"))

  # The sampler works with the two latent variables; the likelihood it reports
  # has to be the mixture the caller asked for, with both integrated out.
  eta <- predict(fit, type = "link", draws = TRUE)
  by_hand <- vapply(seq_len(nrow(eta[[1L]])), function(s) {
    pi <- stats::plogis(eta[[2L]][s, ])
    mu <- exp(eta[[1L]][s, ])
    sum(log(ifelse(d$y == 0, pi + (1 - pi) * stats::dpois(0, mu),
                   (1 - pi) * stats::dpois(d$y, mu))))
  }, numeric(1L))

  expect_equal(fit[["loglik"]], by_hand, tolerance = 1e-8)
  expect_predictor_invariant(fit, d)

  # And the negative binomial version, whose dispersion is still drawn from the
  # true zero-inflated likelihood rather than from the augmented one.
  set.seed(1152)
  d$ynb <- ifelse(structural, 0L,
                  stats::rnbinom(nrow(d), mu = exp(0.7 + d$x1), size = 3))

  fixed <- bartisan(ynb ~ x1 + x2 + x3, data = d, family = zi_negbin(theta = 3),
                    control = quick_control(augment = "zinb"))

  eta <- predict(fixed, type = "link", draws = TRUE)
  by_hand <- vapply(seq_len(nrow(eta[[1L]])), function(s) {
    pi <- stats::plogis(eta[[2L]][s, ])
    mu <- exp(eta[[1L]][s, ])
    zero <- stats::dnbinom(0, mu = mu, size = 3)
    sum(log(ifelse(d$ynb == 0, pi + (1 - pi) * zero,
                   (1 - pi) * stats::dnbinom(d$ynb, mu = mu, size = 3))))
  }, numeric(1L))

  expect_equal(fixed[["loglik"]], by_hand, tolerance = 1e-8)

  drawn <- bartisan(ynb ~ x1 + x2 + x3, data = d, family = zi_negbin(),
                    control = quick_control(augment = "zinb"))
  expect_identical(colnames(drawn[["aux"]]), "theta")
  expect_true(all(drawn[["aux"]][, "theta"] > 0))
})

test_that("the zero-inflated rewriting targets the same posterior", {
  skip_on_cran()

  set.seed(153)
  n <- 500
  d <- as.data.frame(matrix(stats::runif(n * 4), n))
  names(d) <- paste0("x", 1:4)
  count_mean <- 0.4 + 2 * sin(pi * d$x1) - d$x2
  zero_prob <- stats::plogis(-0.8 + d$x2)
  set.seed(1153)
  d$y <- ifelse(stats::runif(n) < zero_prob, 0L,
                stats::rpois(n, exp(count_mean)))
  truth <- exp(count_mean) * (1 - zero_prob)

  ctrl <- function(a) {
    bartisan_control(num_trees = 50, num_burn = 800, num_save = 800,
                     augment = a, verbose = FALSE)
  }

  set.seed(5)
  direct <- bartisan(y ~ ., data = d, family = zi_poisson(),
                     control = ctrl(FALSE))
  set.seed(5)
  rewritten <- bartisan(y ~ ., data = d, family = zi_poisson(),
                        control = ctrl("zip"))

  mean_direct <- predict(direct, type = "response")
  mean_rewritten <- predict(rewritten, type = "response")

  expect_gt(stats::cor(mean_direct, mean_rewritten), 0.95)
  expect_equal(sqrt(mean((mean_rewritten - truth)^2)),
               sqrt(mean((mean_direct - truth)^2)), tolerance = 0.25)

  # The posterior spread has to agree too: a rewriting that sped things up by
  # understating the uncertainty would pass a comparison of point estimates.
  spread <- function(fit) {
    mean(apply(predict(fit, type = "response", draws = TRUE), 2L, stats::sd))
  }
  expect_equal(spread(rewritten), spread(direct), tolerance = 0.3)
})

test_that("the multinomial rewriting works under either kind of rule", {
  d <- sim_x(n = 150, seed = 155)
  set.seed(1155)
  d$y <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))

  for (gate in c("smoothstep", "hard")) {
    fit <- bartisan(y ~ ., data = d, family = multinomial(),
                    control = quick_control(augment = "multinomial",
                                            gate = gate))

    probs <- predict(fit, type = "prob", draws = TRUE)
    observed <- match(d$y, dimnames(probs)[[3L]])
    by_hand <- vapply(seq_len(dim(probs)[1L]), function(s) {
      one <- probs[s, , , drop = FALSE][1L, , ]
      sum(log(one[cbind(seq_along(observed), observed)]))
    }, numeric(1L))

    expect_equal(fit[["loglik"]], by_hand, tolerance = 1e-8)
    expect_predictor_invariant(fit, d)
  }
})
