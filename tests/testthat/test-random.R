# Random intercepts. The thing to check is not that a fit runs but that the
# intercepts are the parameters they claim to be: recovered from data that has
# them, on the right scale, in the right forest, and consistent with the stored
# forest when replayed.

test_that("a random intercept is recovered, and its scale with it", {
  skip_on_cran()

  set.seed(101)
  ng <- 60L
  n <- ng * 10L
  g <- factor(rep(seq_len(ng), each = 10))
  b <- stats::rnorm(ng, sd = 1.2)
  x <- matrix(stats::runif(n * 4), n, 4)
  colnames(x) <- paste0("x", 1:4)
  d <- as.data.frame(x)
  d$g <- g
  f <- 2 * sin(pi * x[, 1]) - x[, 2]
  d$y <- f + b[as.integer(g)] + stats::rnorm(n)

  fit <- bartisan(y ~ x1 + x2 + x3 + x4 + (1 | g), d, num_trees = 50,
                  num_burn = 400, num_draws = 400)

  expect_identical(dim(fit[["ranef"]][[1L]]), c(400L, ng))
  expect_identical(dim(fit[["tau"]][[1L]]), c(400L, 1L))
  expect_identical(colnames(fit[["tau"]][[1L]]), "g")

  est <- colMeans(fit[["ranef"]][[1L]])
  expect_gt(stats::cor(est, b), 0.85)
  # Shrunk towards zero, so the fitted spread is at or below the true one.
  expect_lt(stats::sd(est), 1.15 * stats::sd(b))
  expect_gt(stats::sd(est), 0.6 * stats::sd(b))
  expect_equal(mean(fit[["tau"]][[1L]]), 1.2, tolerance = 0.3)

  # The stored forest plus the stored intercepts reproduce the predictor the
  # sampler recorded, which is what says the random part is in the right place.
  expect_predictor_invariant(fit, d)

  # And it is worth having here: without it the group variation is unexplained.
  plain <- bartisan(y ~ x1 + x2 + x3 + x4, d, num_trees = 50, num_burn = 400,
                    num_draws = 400)
  truth <- f + b[as.integer(g)]
  expect_lt(sqrt(mean((fit[["fitted"]] - truth)^2)),
            0.6 * sqrt(mean((plain[["fitted"]] - truth)^2)))
})

test_that("several grouping factors each get their own scale", {
  skip_on_cran()

  set.seed(102)
  n <- 800
  ng <- 40L
  nh <- 25L
  g <- factor(sample(ng, n, TRUE))
  h <- factor(sample(nh, n, TRUE))
  bg <- stats::rnorm(ng, sd = 1)
  bh <- stats::rnorm(nh, sd = 0.6)
  x <- matrix(stats::runif(n * 3), n, 3)
  colnames(x) <- paste0("x", 1:3)
  d <- as.data.frame(x)
  d$g <- g
  d$h <- h
  d$y <- 2 * sin(pi * x[, 1]) + bg[as.integer(g)] + bh[as.integer(h)] +
    stats::rnorm(n)

  fit <- bartisan(y ~ x1 + x2 + x3 + (1 | g) + (1 | h), d, num_trees = 50,
                  num_burn = 400, num_draws = 400)

  expect_identical(names(fit[["random"]]), c("g", "h"))
  expect_identical(ncol(fit[["ranef"]][[1L]]), ng + nh)
  expect_identical(colnames(fit[["tau"]][[1L]]), c("g", "h"))

  est <- colMeans(fit[["ranef"]][[1L]])
  expect_gt(stats::cor(est[seq_len(ng)], bg), 0.85)
  expect_gt(stats::cor(est[ng + seq_len(nh)], bh), 0.8)

  # The two scales are separately estimated and ordered as the truth is.
  taus <- colMeans(fit[["tau"]][[1L]])
  expect_gt(taus[["g"]], taus[["h"]])

  expect_predictor_invariant(fit, d)
})

test_that("nesting expands as it does in lme4", {
  d <- sim_x(n = 200, seed = 103)
  set.seed(1103)
  d$a <- factor(sample(5, nrow(d), TRUE))
  d$b <- factor(sample(3, nrow(d), TRUE))
  d$y <- stats::rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + (1 | a / b), d, control = quick_control())

  # `(1 | a/b)` is two intercept terms, the inner one being the interaction.
  expect_identical(length(fit[["random"]]), 2L)
  expect_true("a" %in% names(fit[["random"]]))
  expect_predictor_invariant(fit, d)
})

test_that("a random slope is refused, and says what to do instead", {
  d <- sim_x(n = 100, seed = 104)
  d$g <- factor(sample(5, nrow(d), TRUE))
  d$y <- stats::rnorm(nrow(d))

  expect_error(bartisan(y ~ x1 + (1 + x1 | g), d, control = quick_control()),
               "[Oo]nly random intercepts")
  expect_error(bartisan(y ~ x1 + (x1 | g), d, control = quick_control()),
               "random slope")

  # A grouping factor with one level has nothing to vary over.
  d$one <- factor("a")
  expect_error(bartisan(y ~ x1 + (1 | one), d, control = quick_control()),
               "nothing for a group effect")

  # And a missing group cannot be given a group effect.
  d$h <- factor(sample(4, nrow(d), TRUE))
  d$h[3] <- NA
  expect_error(bartisan(y ~ x1 + (1 | h), d, control = quick_control()),
               "missing values")
})

test_that("the grouping factor is not also a predictor", {
  d <- sim_x(n = 200, seed = 105)
  set.seed(1105)
  d$g <- factor(sample(6, nrow(d), TRUE))
  d$y <- stats::rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + x2 + (1 | g), d, control = quick_control())

  # `g` is in the random part, so no splitting rule may use it.
  expect_false("g" %in% fit[["group_names"]])
  expect_setequal(fit[["term_labels"]], c("x1", "x2"))

  # And `y ~ . + (1 | g)` puts g in both, which is the caller's business but has
  # to work: the dot expands over the data, g included.
  both <- bartisan(y ~ . - g + (1 | g), d, control = quick_control())
  expect_false("g" %in% both[["group_names"]])
  expect_predictor_invariant(both, d)
})

test_that("each additive predictor gets its own set of intercepts", {
  skip_on_cran()

  set.seed(106)
  n <- 1000
  ng <- 40
  g <- factor(sample(ng, n, TRUE))
  b1 <- stats::rnorm(ng, sd = 0.8)
  b2 <- stats::rnorm(ng, sd = 0.5)
  x <- matrix(stats::runif(n * 3), n, 3)
  colnames(x) <- paste0("x", 1:3)
  d <- as.data.frame(x)
  d$g <- g

  # A group effect on the mean and a different one on the log standard
  # deviation, which is the case that says the two sets are separate.
  mu <- 2 * sin(pi * x[, 1]) + b1[as.integer(g)]
  lsd <- -0.5 + 0.5 * x[, 2] + b2[as.integer(g)]
  d$y <- stats::rnorm(n, mu, exp(lsd))

  fit <- bartisan(y ~ x1 + x2 + x3 + (1 | g), d, family = location_scale(),
                  num_trees = 30, num_burn = 400, num_draws = 400)

  expect_identical(length(fit[["ranef"]]), 2L)
  expect_identical(length(fit[["tau"]]), 2L)
  expect_gt(stats::cor(colMeans(fit[["ranef"]][[1L]]), b1), 0.8)
  expect_gt(stats::cor(colMeans(fit[["ranef"]][[2L]]), b2), 0.8)

  # Each set is matched to its own predictor, not swapped.
  expect_gt(stats::cor(colMeans(fit[["ranef"]][[1L]]), b1),
            stats::cor(colMeans(fit[["ranef"]][[1L]]), b2))
  expect_gt(stats::cor(colMeans(fit[["ranef"]][[2L]]), b2),
            stats::cor(colMeans(fit[["ranef"]][[2L]]), b1))

  expect_predictor_invariant(fit, d)
})

test_that("the intercepts work through the general and exponential target paths", {
  skip_on_cran()

  # The leaf update has three paths -- a closed form for a quadratic target, a
  # scalar iteration for an exponential one, a Laplace-and-Metropolis step
  # otherwise -- and a random intercept goes through whichever the family calls
  # for. One family per path.
  set.seed(107)
  n <- 800
  ng <- 40
  g <- factor(sample(ng, n, TRUE))
  b <- stats::rnorm(ng, sd = 0.8)
  x <- matrix(stats::runif(n * 2), n, 2)
  colnames(x) <- paste0("x", 1:2)
  d <- as.data.frame(x)
  d$g <- g
  lin <- 2 * sin(pi * x[, 1]) + b[as.integer(g)]

  d$yb <- stats::rbinom(n, 1, stats::pnorm(lin))          # quadratic, augmented
  d$yc <- stats::rpois(n, exp(1 + 0.5 * lin))             # exponential form
  d$yo <- ordered(rowSums(outer(lin + stats::rlogis(n),
                                c(-0.5, 1), ">")) + 1L)   # general

  fits <- list(
    probit = bartisan(yb ~ x1 + x2 + (1 | g), d, family = binomial("probit"),
                      gate = "hard", num_trees = 30, num_burn = 300,
                     num_draws = 300),
    poisson = bartisan(yc ~ x1 + x2 + (1 | g), d, family = poisson(),
                       gate = "hard", num_trees = 30, num_burn = 300,
                      num_draws = 300),
    ordinal = bartisan(yo ~ x1 + x2 + (1 | g), d, family = ordinal("logit"),
                       gate = "hard", augment = FALSE, num_trees = 30,
                      num_burn = 300, num_draws = 300))

  for (nm in names(fits)) {
    est <- colMeans(fits[[nm]][["ranef"]][[1L]])
    expect_gt(stats::cor(est, b), 0.6, label = nm)
    expect_predictor_invariant(fits[[nm]], d)
  }
})

test_that("a level not seen at fitting time gets the prior mean", {
  d <- sim_x(n = 200, seed = 108)
  set.seed(1108)
  d$g <- factor(sample(letters[1:6], nrow(d), TRUE))
  d$y <- stats::rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + x2 + (1 | g), d, control = quick_control())

  nd <- d[1:5, ]
  nd$g <- factor("z", levels = c(levels(d$g), "z"))

  expect_warning(p <- stats::predict(fit, newdata = nd, type = "link"),
                 "not present when the model was fit")
  expect_true(all(is.finite(p)))

  # Two copies of the *same* row, differing only in the group, differ by exactly
  # that group's intercept -- so the unseen level really did get zero. The
  # covariates have to be identical or the forest contributes to the difference
  # as well.
  mixed <- d[c(1L, 1L), ]
  mixed$g <- factor(c(levels(d$g)[1L], "z"), levels = c(levels(d$g), "z"))
  suppressWarnings(pm <- stats::predict(fit, newdata = mixed, type = "link",
                                        draws = TRUE))
  b_first <- fit[["ranef"]][[1L]][, 1L]
  expect_equal(as.vector(pm[, 1L] - pm[, 2L]), as.vector(b_first),
               tolerance = 1e-8)
})

test_that("chains pool the intercepts and diagnose them", {
  skip_on_cran()

  d <- sim_x(n = 300, seed = 109)
  set.seed(1109)
  d$g <- factor(sample(12, nrow(d), TRUE))
  d$y <- d$x1 + stats::rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + x2 + (1 | g), d, chains = 3L, num_trees = 20L,
                  num_burn = 150L, num_draws = 150L)

  expect_identical(nrow(fit[["ranef"]][[1L]]), 450L)
  expect_identical(nrow(fit[["tau"]][[1L]]), 450L)

  # The scale and the intercepts both appear in the diagnostics table.
  quantities <- fit[["rhat"]]$quantity
  expect_true(any(grepl("^tau\\.", quantities)))
  expect_true(any(grepl("^ranef\\.", quantities)))
  expect_predictor_invariant(fit, d)
})

test_that("update_tau = FALSE holds the scale where it started", {
  d <- sim_x(n = 200, seed = 110)
  set.seed(1110)
  d$g <- factor(sample(8, nrow(d), TRUE))
  d$y <- stats::rnorm(nrow(d))

  fit <- bartisan(y ~ x1 + (1 | g), d, control = quick_control(update_tau = FALSE))
  expect_identical(stats::sd(as.vector(fit[["tau"]][[1L]])), 0)

  drawn <- bartisan(y ~ x1 + (1 | g), d, control = quick_control())
  expect_gt(stats::sd(as.vector(drawn[["tau"]][[1L]])), 0)
})

test_that("a formula with no bars is unaffected", {
  d <- sim_x(n = 100, seed = 111)
  d$y <- stats::rnorm(nrow(d))

  set.seed(5)
  a <- bartisan(y ~ x1 + x2, d, control = quick_control())

  expect_null(a[["random"]])
  expect_null(a[["ranef"]])
  expect_null(a[["tau"]])
  expect_null(summary(a)[["tau"]])
})
