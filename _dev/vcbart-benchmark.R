# bartisan's varying coefficients against VCBART's own.
#
# The design carries four coefficient functions of deliberately different kinds:
# one strongly varying, one mildly varying, one constant and nonzero, and one
# identically zero. The last two are the interesting ones. A constant
# coefficient is what a varying-coefficient model has to *not* invent structure
# in, and a null one is where Ghosh, Bhogale and Deshpande (2025) claim a
# global-local prior gives "substantially narrower and better-calibrated"
# intervals than the plain VCBART prior -- so the null column is also the
# measurement of whether bartisan's own sparsity settings do enough.
#
# Paired within replicate; the standard error reported is of the paired
# difference, not of either method's own spread.

library(bartisan)
suppressMessages(library(VCBART))

REPS <- 8L
N <- 500L
BURN <- 1000L
DRAWS <- 1000L
TREES <- 50L

# beta_0 is the control function; beta_1 varies strongly, beta_2 mildly,
# beta_3 is constant, beta_4 is zero everywhere.
betas <- function(z) {
  cbind(3 + 2 * z[, 1] - 3 * z[, 2]^2,
        2 * sin(pi * z[, 1]),
        z[, 3],
        rep.int(1, nrow(z)),
        rep.int(0, nrow(z)))
}

LABELS <- c("beta0 (control)", "beta1 (strong)", "beta2 (mild)",
            "beta3 (constant 1)", "beta4 (null 0)")

simulate <- function(rep) {
  set.seed(1000 + rep)
  z <- matrix(runif(N * 5, -1, 1), N, 5, dimnames = list(NULL, paste0("z", 1:5)))
  x <- matrix(rnorm(N * 4), N, 4, dimnames = list(NULL, paste0("x", 1:4)))
  b <- betas(z)
  mu <- b[, 1] + rowSums(x * b[, -1, drop = FALSE])
  list(z = z, x = x, b = b, mu = mu, y = mu + rnorm(N))
}

# Recovery of each coefficient function, the coverage of its 95% interval, and
# how wide that interval is. Width only means something next to coverage, which
# is why the two are reported together.
score <- function(draws, truth) {
  est <- colMeans(draws)
  lo <- apply(draws, 2L, stats::quantile, 0.025)
  hi <- apply(draws, 2L, stats::quantile, 0.975)

  c(rmse = sqrt(mean((est - truth)^2)),
    coverage = mean(truth >= lo & truth <= hi),
    width = mean(hi - lo))
}

run_rep <- function(rep) {
  d <- simulate(rep)
  frame <- data.frame(y = d[["y"]], d[["x"]], d[["z"]])
  modifiers <- ~ z1 + z2 + z3 + z4 + z5

  form <- stats::reformulate(
    c(paste0("z", 1:5),
      sprintf('vc(x%d, ~ z1 + z2 + z3 + z4 + z5, center = "zero")', 1:4)),
    response = quote(y))

  t_bartisan <- system.time({
    fit <- bartisan(form, data = frame, family = gaussian(),
                    num_trees = TREES, num_burn = BURN, num_draws = DRAWS,
                    verbose = FALSE)
  })[["elapsed"]]

  # The control function is the prediction with every covariate at zero, which
  # is what `center = "zero"` makes it.
  at_zero <- frame
  at_zero[paste0("x", 1:4)] <- 0
  b_bartisan <- c(list(predict(fit, newdata = at_zero, type = "link",
                              draws = TRUE)),
                  unname(coef(fit, draws = TRUE)))

  t_vcbart <- system.time({
    vfit <- VCBART::VCBART_ind(Y_train = d[["y"]],
                               subj_id_train = seq_len(N),
                               ni_train = rep.int(1L, N),
                               X_train = d[["x"]],
                               Z_cont_train = d[["z"]],
                               M = TREES, nd = DRAWS, burn = BURN,
                               save_samples = TRUE, save_trees = FALSE,
                               verbose = FALSE)
  })[["elapsed"]]

  b_vcbart <- asplit(vfit[["betahat.train"]], 3L)

  out <- vapply(seq_len(5L), function(j) {
    c(bartisan = score(b_bartisan[[j]], d[["b"]][, j]),
      VCBART = score(b_vcbart[[j]], d[["b"]][, j]))
  }, numeric(6L))

  colnames(out) <- LABELS
  list(scores = out, time = c(bartisan = t_bartisan, VCBART = t_vcbart))
}

results <- lapply(seq_len(REPS), run_rep)
scores <- simplify2array(lapply(results, `[[`, "scores"))
times <- vapply(results, `[[`, numeric(2L), "time")

cat(sprintf("\nVCBART comparison: %d replicates, n = %d, %d trees, %d draws\n",
            REPS, N, TREES, DRAWS))

for (metric in c("rmse", "coverage", "width")) {
  cat(sprintf("\n%s\n", switch(metric,
    rmse = "Root mean squared error of the coefficient function",
    coverage = "Coverage of the 95% interval (nominal 0.95)",
    width = "Mean width of the 95% interval")))
  cat(sprintf("  %-20s %10s %10s %12s %9s\n",
              "coefficient", "bartisan", "VCBART", "difference", "se"))

  for (j in seq_along(LABELS)) {
    a <- scores[paste0("bartisan.", metric), j, ]
    b <- scores[paste0("VCBART.", metric), j, ]
    cat(sprintf("  %-20s %10.4f %10.4f %+12.4f %9.4f\n", LABELS[j],
                mean(a), mean(b), mean(a - b),
                stats::sd(a - b) / sqrt(REPS)))
  }
}

cat(sprintf("\nElapsed seconds per fit: bartisan %.1f, VCBART %.1f\n",
            mean(times["bartisan", ]), mean(times["VCBART", ])))

saveRDS(list(scores = scores, times = times), "_dev/vcbart-benchmark.rds")
