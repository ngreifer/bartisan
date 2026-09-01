# Do bartisan's existing sparsity settings do what a global-local prior would?
#
# Ghosh, Bhogale and Deshpande (2025) add two things to VCBART: a DART sparsity
# prior over the modifiers, which bartisan already has per forest, and a
# regularized horseshoe on the leaf jumps, which bartisan has only half of -- a
# per-forest half-Cauchy scale, with no shared global parameter tying the
# ensembles together and no slab. Their headline claim is narrower and
# better-calibrated intervals, "especially for null covariate effects".
#
# So the measurement that matters is the null coefficient: how wide its interval
# is, whether it covers, and whether the sparsity prior is what closes the gap.
# Width alone proves nothing -- an interval of zero width covers nothing -- so
# the two are always read together.

library(bartisan)

REPS <- 10L
N <- 500L
BURN <- 1000L
DRAWS <- 1000L

betas <- function(z) {
  cbind(3 + 2 * z[, 1] - 3 * z[, 2]^2,
        2 * sin(pi * z[, 1]),
        z[, 3],
        rep.int(1, nrow(z)),
        rep.int(0, nrow(z)))
}

LABELS <- c("beta1 (strong)", "beta2 (mild)", "beta3 (constant 1)",
            "beta4 (null 0)")

ARMS <- list(
  `sparsity = TRUE (default)` = list(sparsity = TRUE),
  `sparsity = FALSE` = list(sparsity = FALSE),
  `sparsity = "strong"` = list(sparsity = "strong"),
  `sparsity = "weak"` = list(sparsity = "weak")
)

score <- function(draws, truth) {
  est <- colMeans(draws)
  lo <- apply(draws, 2L, stats::quantile, 0.025)
  hi <- apply(draws, 2L, stats::quantile, 0.975)

  c(rmse = sqrt(mean((est - truth)^2)),
    coverage = mean(truth >= lo & truth <= hi),
    width = mean(hi - lo))
}

run_rep <- function(rep) {
  set.seed(1000 + rep)
  z <- matrix(runif(N * 5, -1, 1), N, 5, dimnames = list(NULL, paste0("z", 1:5)))
  x <- matrix(rnorm(N * 4), N, 4, dimnames = list(NULL, paste0("x", 1:4)))
  b <- betas(z)
  y <- b[, 1] + rowSums(x * b[, -1, drop = FALSE]) + rnorm(N)
  frame <- data.frame(y = y, x, z)

  form <- stats::reformulate(
    c(paste0("z", 1:5),
      sprintf('vc(x%d, ~ z1 + z2 + z3 + z4 + z5, center = "zero")', 1:4)),
    response = quote(y))

  vapply(ARMS, function(arm) {
    fit <- do.call(bartisan, c(
      list(formula = form, data = frame, family = gaussian(),
           num_trees = 50L, num_burn = BURN, num_draws = DRAWS,
           verbose = FALSE),
      arm))

    cf <- coef(fit, draws = TRUE)
    as.vector(vapply(seq_len(4L), function(j) score(cf[[j]], b[, j + 1L]),
                     numeric(3L)))
  }, numeric(12L))
}

each <- vapply(seq_len(REPS), run_rep, matrix(0, 12L, length(ARMS)))

cat(sprintf("\nSparsity and the null coefficient: %d replicates, n = %d\n",
            REPS, N))

for (m in seq_len(3L)) {
  metric <- c("rmse", "coverage", "width")[m]
  cat(sprintf("\n%s\n", switch(metric,
    rmse = "Root mean squared error",
    coverage = "Coverage of the 95% interval (nominal 0.95)",
    width = "Mean width of the 95% interval")))
  cat(sprintf("  %-26s", "arm"))
  cat(sprintf("%20s", LABELS), "\n")

  for (a in seq_along(ARMS)) {
    cat(sprintf("  %-26s", names(ARMS)[a]))
    for (j in seq_len(4L)) {
      v <- each[(j - 1L) * 3L + m, a, ]
      base <- each[(j - 1L) * 3L + m, 1L, ]
      cat(sprintf("%12.3f%s", mean(v),
                  if (a == 1L) "        "
                  else sprintf(" (%+.3f)", mean(v - base))))
    }
    cat("\n")
  }
}

cat("\nParentheses are the paired difference from the default arm.\n")
cat("Standard errors of those differences, for the null coefficient:\n")
for (m in seq_len(3L)) {
  for (a in 2:length(ARMS)) {
    d <- each[9L + m, a, ] - each[9L + m, 1L, ]
    cat(sprintf("  %-12s %-26s %+.4f  (se %.4f)\n",
                c("rmse", "coverage", "width")[m], names(ARMS)[a],
                mean(d), stats::sd(d) / sqrt(REPS)))
  }
}

saveRDS(each, "_dev/sparsity-null.rds")
