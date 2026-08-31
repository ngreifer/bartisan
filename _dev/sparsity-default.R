# Should `sparsity = TRUE` stay the default?
#
# The atom-at-zero finding says the sparsity prior costs something for a
# contrast estimand. It does not say what the prior buys elsewhere, and the
# default should turn on the trade rather than on either half of it. Two
# questions, measured separately because they have different answers:
#
#   (1) Prediction. Does concentrating splits improve accuracy, and does the
#       answer depend on how many irrelevant predictors there are?
#   (2) A contrast. What does each setting do to the ATE: its bias, its interval,
#       and the mass the prior puts at exactly zero?
#
# Everything below uses default `num_burn` and `num_draws` so that this script
# does not depend on their names.

library(bartisan)
library(marginaleffects)

LEVELS <- c("none", "weak", "moderate", "strong")
REPS <- 3

friedman <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] + 5 * x[, 5]
}

# ---- (1) prediction ------------------------------------------------------

predict_rep <- function(p, rep) {
  set.seed(1000 * p + rep)
  n <- 500
  xtr <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  xte <- matrix(runif(1000 * p), 1000, p, dimnames = list(NULL, colnames(xtr)))
  ftr <- friedman(xtr)
  fte <- friedman(xte)
  dtr <- data.frame(y = ftr + rnorm(n), xtr)
  dte <- as.data.frame(xte)

  vapply(LEVELS, function(lv) {
    set.seed(7)
    fit <- bartisan(y ~ ., data = dtr, family = gaussian(), chains = 2,
                    control = bartisan_control(sparsity = lv, verbose = FALSE))
    # Scored against the true regression function, not the noisy outcome, so
    # this is error in the fit rather than irreducible noise.
    sqrt(mean((predict(fit, newdata = dte) - fte)^2))
  }, numeric(1))
}

cat("== Prediction: RMSE against the true function, mean of", REPS, "reps ==\n")
cat(sprintf("%-8s %s\n", "p", paste(sprintf("%8s", LEVELS), collapse = "")))
for (p in c(10L, 50L)) {
  m <- rowMeans(vapply(seq_len(REPS), function(r) predict_rep(p, r),
                       numeric(length(LEVELS))))
  cat(sprintf("%-8d %s\n", p, paste(sprintf("%8.3f", m), collapse = "")))
}

# ---- (2) a contrast ------------------------------------------------------

TAU <- 0.5

contrast_rep <- function(p, rep) {
  set.seed(2000 * p + rep)
  n <- 800
  x <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  # Confounded assignment, and a constant treatment effect so the truth is one
  # number the interval either covers or does not.
  z <- rbinom(n, 1, plogis(2 * (x[, 1] - 0.5)))
  d <- data.frame(y = friedman(x) / 5 + TAU * z + rnorm(n), z = z, x)

  t(vapply(LEVELS, function(lv) {
    set.seed(7)
    fit <- bartisan(y ~ ., data = d, family = gaussian(), chains = 2,
                    control = bartisan_control(sparsity = lv, verbose = FALSE))
    a <- avg_comparisons(fit, variables = "z")
    post <- attr(a, "posterior_draws")
    vi <- variable_importance(fit)

    c(est = a$estimate, lo = a$conf.low, hi = a$conf.high,
      width = a$conf.high - a$conf.low,
      covers = as.numeric(a$conf.low <= TAU && a$conf.high >= TAU),
      atom = if (is.null(post)) NA_real_ else mean(post == 0),
      prop_used = vi[["prop_used"]][vi[["variable"]] == "z"])
  }, numeric(7)))
}

for (p in c(10L, 50L)) {
  acc <- Reduce(`+`, lapply(seq_len(REPS), function(r) contrast_rep(p, r))) / REPS
  cat(sprintf("\n== A contrast, p = %d, truth %.2f, mean of %d reps ==\n",
              p, TAU, REPS))
  cat(sprintf("%-10s %7s %7s %7s %7s %7s %9s\n",
              "sparsity", "est", "lo", "hi", "width", "atom", "prop_used"))
  for (i in seq_along(LEVELS)) {
    cat(sprintf("%-10s %7.3f %7.3f %7.3f %7.3f %7.3f %9.3f\n",
                LEVELS[i], acc[i, "est"], acc[i, "lo"], acc[i, "hi"],
                acc[i, "width"], acc[i, "atom"], acc[i, "prop_used"]))
  }
  cat(sprintf("coverage of the truth: %s\n",
              paste(sprintf("%s %.2f", LEVELS, acc[, "covers"]), collapse = "  ")))
}
