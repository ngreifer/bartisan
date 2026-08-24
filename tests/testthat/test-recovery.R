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

  fit <- genbart(y ~ ., data = d, weights = rep(1e-10, n),
                 control = genbart_control(num_trees = num_trees,
                                           num_burn = 1000, num_save = 2000,
                                           soft = FALSE, sigma_mu = 0.4,
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

  fit <- genbart(y ~ ., data = d, weights = rep(1e-10, n),
                 control = genbart_control(num_trees = num_trees,
                                           num_burn = 1000, num_save = 2000,
                                           soft = FALSE, sigma_mu = sigma_mu,
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

  ctrl <- genbart_control(num_trees = 20, num_burn = 400, num_save = 400,
                          verbose = FALSE)

  d$y <- truth + stats::rnorm(n, sd = 0.3)
  fit <- genbart(y ~ ., data = d, control = ctrl)
  expect_gt(stats::cor(predict(fit), truth), 0.9)

  d$y <- stats::rpois(n, exp(truth))
  fit <- genbart(y ~ ., data = d, family = poisson(), control = ctrl)
  expect_gt(stats::cor(predict(fit), truth), 0.8)

  d$y <- stats::rgamma(n, 4, rate = 4 / exp(truth))
  fit <- genbart(y ~ ., data = d, family = Gamma("log"), control = ctrl)
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

  fit <- genbart(y ~ ., data = d,
                 control = genbart_control(num_trees = 20, num_burn = 400,
                                           num_save = 400, verbose = FALSE))

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

  fit <- genbart(y ~ ., data = d, family = ordinal(),
                 control = genbart_control(num_trees = 20, num_burn = 500,
                                           num_save = 500, verbose = FALSE))

  # Only the gaps between cutpoints are identified; the level is absorbed by
  # the intercept, which is why the first cutpoint is pinned at zero.
  gaps <- colMeans(fit[["aux"]])[-1L]

  expect_equal(unname(gaps), cuts[-1L] - cuts[1L], tolerance = 0.4)
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

  fit <- genbart(survival::Surv(time, status) ~ x1 + x2, data = d,
                 family = loglogistic_aft(),
                 control = genbart_control(num_trees = 20, num_burn = 500,
                                           num_save = 500, verbose = FALSE))

  expect_equal(mean(fit[["aux"]][, "sigma"]), 0.7, tolerance = 0.2)
})

test_that("irrelevant predictors are used less than relevant ones", {
  skip_on_cran()

  set.seed(107)

  n <- 300
  X <- matrix(stats::runif(n * 8), nrow = n)
  colnames(X) <- paste0("x", 1:8)
  d <- as.data.frame(X)
  d$y <- 3 * X[, 1] + 3 * X[, 2] + stats::rnorm(n, sd = 0.3)

  fit <- genbart(y ~ ., data = d,
                 control = genbart_control(num_trees = 20, num_burn = 500,
                                           num_save = 500, verbose = FALSE))

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
  fit <- genbart(y ~ ., data = d, weights = rep(1e-10, n),
                 control = genbart_control(num_trees = 1, num_burn = 2000,
                                           num_save = 20000, soft = FALSE,
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
    genbart(y ~ ., data = d, family = family, weights = rep(1e-10, n),
            control = genbart_control(num_trees = 20, num_burn = 300,
                                      num_save = 1500, soft = FALSE,
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
