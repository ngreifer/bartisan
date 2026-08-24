# The sampler's tree bookkeeping was reworked for speed: nodes are recycled
# rather than allocated per proposal, the child weights come out of the split
# that produces them, the supports are snapshotted so a rejected bandwidth move
# can be undone by copying, and a tree with no splits draws its bandwidth from
# its prior directly. None of that is supposed to change what the sampler
# computes, so this file checks the things that would break if it did.

test_that("recycling nodes leaves the fit unchanged", {
  # Every field a recycled node starts with comes from the same function a fresh
  # one uses, so a stale field is the failure mode to look for. It would show up
  # as the stored forest disagreeing with the predictor the sampler recorded --
  # a wrong cutpoint, variable or leaf value replays differently.
  d <- sim_x(n = 200, seed = 51)
  set.seed(1051)
  d$y <- 2 * sin(pi * d$x1) - d$x2 + stats::rnorm(nrow(d))

  for (soft in c(TRUE, FALSE)) {
    fit <- genbart(y ~ ., d, control = quick_control(
      soft = soft, num_trees = 25L, num_burn = 200L, num_save = 200L))
    expect_predictor_invariant(fit, d)
  }

  # A long run with many births and deaths, so the pool is exercised heavily,
  # and with a deep tree prior so nodes are taken and given at several depths.
  deep <- genbart(y ~ ., d, control = quick_control(
    soft = FALSE, gamma = 0.99, beta = 0.5, num_trees = 10L,
    num_burn = 300L, num_save = 300L))
  expect_predictor_invariant(deep, d)
})

test_that("the bandwidth of a tree with no splits is drawn from its prior", {
  skip_on_cran()

  # With the branching probability at zero the tree-shape prior forbids splits,
  # so every tree stays a single leaf, the bandwidth never enters the likelihood,
  # and its conditional is exactly the Exp(bandwidth) prior. That makes this the
  # one part of the sampler with a closed-form answer to check against.
  set.seed(61)
  d <- data.frame(x1 = stats::runif(300), x2 = stats::runif(300))
  d$y <- stats::rnorm(300)

  fit <- genbart(y ~ ., d, gamma = 1e-10, bandwidth = 0.1, num_trees = 20L,
                 num_burn = 200L, num_save = 2000L)

  b <- as.vector(fit[["bandwidth"]])
  expect_gt(length(b), 1000)

  # Mean and standard deviation of Exp(mean = 0.1) are both 0.1.
  expect_equal(mean(b), 0.1, tolerance = 0.03)
  expect_equal(stats::sd(b), 0.1, tolerance = 0.05)

  # And the whole distribution, not just its first two moments.
  expect_gt(suppressWarnings(stats::ks.test(b, "pexp", 10))$p.value, 0.01)

  # A different prior scale moves it, so the draw is using the prior rather than
  # anything hard-coded.
  wide <- genbart(y ~ ., d, gamma = 1e-10, bandwidth = 0.4, num_trees = 20L,
                  num_burn = 200L, num_save = 2000L)
  expect_equal(mean(as.vector(wide[["bandwidth"]])), 0.4, tolerance = 0.03)
})

test_that("a rejected bandwidth move leaves the tree exactly as it was", {
  # The snapshot has to be a faithful stand-in for recomputing every gate. If it
  # were not, the memberships and the predictor would drift apart, which is what
  # the predictor invariant detects.
  d <- sim_x(n = 250, seed = 71)
  set.seed(1071)
  d$y <- 2 * d$x1 * d$x2 + stats::rnorm(nrow(d))

  for (gate in c("logistic", "smoothstep", "smootherstep")) {
    fit <- genbart(y ~ ., d, control = quick_control(
      gate = gate, num_trees = 20L, num_burn = 300L, num_save = 300L))
    expect_predictor_invariant(fit, d)
    # The bandwidth actually moved, so rejections really happened.
    expect_gt(stats::sd(as.vector(fit[["bandwidth"]])), 0)
  }
})

test_that("the child weights the split records are the ones the target needs", {
  # Dividing a node's support and recording what each side got used to be two
  # passes over the same gates. Fusing them is only safe if the recorded weights
  # still sum to the parent's, which is what the whole membership scheme rests
  # on -- if they did not, the predictor and the stored forest would disagree.
  d <- sim_x(n = 200, seed = 81)
  set.seed(1081)
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(2 * d$x1 - d$x2))

  for (soft in c(TRUE, FALSE)) {
    fit <- genbart(y ~ ., d, family = binomial("probit"),
                   control = quick_control(soft = soft, num_trees = 20L,
                                           num_burn = 200L, num_save = 200L))
    expect_predictor_invariant(fit, d)
  }

  # Including with missing values, where a rule decides the side outright and
  # the weight is hard even when the rules are soft. `na.pass` is what keeps the
  # rows, rather than the default `na.omit` dropping them.
  d2 <- d
  d2$x3[c(4, 19, 60, 130)] <- NA
  fit <- genbart(y ~ ., d2, family = binomial("probit"),
                 na.action = stats::na.pass,
                 control = quick_control(num_trees = 20L, num_burn = 200L,
                                         num_save = 200L))
  expect_predictor_invariant(fit, d2)
})
