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

  for (gate in c("smoothstep", "hard")) {
    fit <- bartisan(y ~ ., d, control = quick_control(
      gate = gate, num_trees = 25L, num_burn = 200L, num_draws = 200L))
    expect_predictor_invariant(fit, d)
  }

  # A long run with many births and deaths, so the pool is exercised heavily,
  # and with a deep tree prior so nodes are taken and given at several depths.
  deep <- bartisan(y ~ ., d, control = quick_control(
    gate = "hard", gamma = 0.99, beta = 0.5, num_trees = 10L,
    num_burn = 300L, num_draws = 300L))
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

  fit <- bartisan(y ~ ., d, gamma = 1e-10, bandwidth = 0.1, num_trees = 20L,
                  num_burn = 200L, num_draws = 2000L)

  b <- as.vector(fit[["bandwidth"]])
  expect_gt(length(b), 1000)

  # Mean and standard deviation of Exp(mean = 0.1) are both 0.1.
  expect_equal(mean(b), 0.1, tolerance = 0.03)
  expect_equal(stats::sd(b), 0.1, tolerance = 0.05)

  # And the whole distribution, not just its first two moments.
  expect_gt(suppressWarnings(stats::ks.test(b, "pexp", 10))$p.value, 0.01)

  # A different prior scale moves it, so the draw is using the prior rather than
  # anything hard-coded.
  wide <- bartisan(y ~ ., d, gamma = 1e-10, bandwidth = 0.4, num_trees = 20L,
                   num_burn = 200L, num_draws = 2000L)
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
    fit <- bartisan(y ~ ., d, control = quick_control(
      gate = gate, num_trees = 20L, num_burn = 300L, num_draws = 300L))
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

  for (gate in c("smoothstep", "hard")) {
    fit <- bartisan(y ~ ., d, family = binomial("probit"),
                    control = quick_control(gate = gate, num_trees = 20L,
                                            num_burn = 200L, num_draws = 200L))
    expect_predictor_invariant(fit, d)
  }

  # Including with missing values, where a rule decides the side outright and
  # the weight is hard even when the rules are soft. `na.pass` is what keeps the
  # rows, rather than the default `na.omit` dropping them.
  d2 <- d
  d2$x3[c(4, 19, 60, 130)] <- NA
  fit <- bartisan(y ~ ., d2, family = binomial("probit"),
                  na.action = stats::na.pass,
                 control = quick_control(num_trees = 20L, num_burn = 200L,
                                         num_draws = 200L))
  expect_predictor_invariant(fit, d2)
})

# `get_varnames()` is `all.vars()` except that an extraction with `$`, `[[` or
# `[` is one variable rather than the pieces it is written from, which is how
# `model.frame()` names it.
test_that("get_varnames keeps an extraction whole where all.vars splits it", {
  expect_identical(get_varnames(quote(x)), "x")
  expect_identical(get_varnames(quote(log(x))), "x")
  expect_identical(get_varnames(quote(poly(x, 2))), "x")
  expect_identical(get_varnames(quote(cbind(a, b))), c("a", "b"))

  # The three extraction operators, each of which `all.vars()` mishandles: `$`
  # splits into the object and the column, and `[[` and `[` lose the column.
  expect_identical(get_varnames(quote(d$g)), "d$g")
  expect_identical(get_varnames(quote(d[["g"]])), 'd[["g"]]')
  expect_identical(get_varnames(quote(log(d$g + 1))), "d$g")
  expect_identical(get_varnames(quote(x + d$g)), c("x", "d$g"))

  expect_identical(all.vars(quote(d$g)), c("d", "g"))
  expect_identical(all.vars(quote(d[["g"]])), "d")
})

test_that("get_varnames walks a terms object, which as.list() would derail", {
  d <- data.frame(y = 1:3, x1 = 1:3, x2 = 1:3)
  tt <- stats::terms(y ~ x1 + d$x2, data = d)

  expect_identical(get_varnames(tt), c("y", "x1", "d$x2"))
  expect_identical(get_varnames(stats::delete.response(tt)), c("x1", "d$x2"))

  # A terms object is a call carrying attributes, and `as.list()` on one does
  # not drop the operator the way it does on a plain call, so a traversal built
  # on `as.list(e)[-1L]` walks the wrong tree. Indexing positions avoids it.
  expect_length(as.list(tt)[-1L], 3L)

  # The same formula without the attributes must give the same answer.
  expect_identical(get_varnames(y ~ x1 + d$x2), c("y", "x1", "d$x2"))
})
