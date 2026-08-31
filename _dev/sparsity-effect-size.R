# The first pass found no atom at any sparsity level with a treatment effect of
# 0.5 against residual sd 1: `prop_used` was 1.000 everywhere. The RHC fit has
# `prop_used` 0.898 and an atom of 0.102. So the atom is not a property of the
# sparsity prior alone, it is what the prior does to a predictor whose signal is
# weak enough to be droppable.
#
# That makes effect size the variable to sweep, and it turns the default question
# into a measurable one: when the signal is weak, does the sparsity prior only
# put a bound at exactly zero, or does it also bias the estimate and break the
# interval's coverage? The first is a reporting quirk. The second is a reason to
# change the default.

library(bartisan)
library(marginaleffects)

LEVELS <- c("none", "moderate", "strong")
TAUS <- c(0.05, 0.10, 0.20, 0.50)
REPS <- 5
P <- 20

friedman <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + 10 * x[, 4]
}

one <- function(tau, lv, rep) {
  set.seed(rep * 977 + round(tau * 1000))
  n <- 800
  x <- matrix(runif(n * P), n, P, dimnames = list(NULL, paste0("x", seq_len(P))))
  z <- rbinom(n, 1, plogis(2 * (x[, 1] - 0.5)))
  d <- data.frame(y = friedman(x) / 5 + tau * z + rnorm(n), z = z, x)

  set.seed(7)
  fit <- bartisan(y ~ ., data = d, family = gaussian(),
                  control = bartisan_control(sparsity = lv, chains = 2,
                                             verbose = FALSE))
  a <- avg_comparisons(fit, variables = "z")
  post <- attr(a, "posterior_draws")
  vi <- variable_importance(fit)

  c(est = a$estimate,
    atom = if (is.null(post)) NA_real_ else mean(post == 0),
    lo_is_zero = as.numeric(isTRUE(all.equal(a$conf.low, 0))),
    covers = as.numeric(a$conf.low <= tau && a$conf.high >= tau),
    width = a$conf.high - a$conf.low,
    prop_used = vi[["prop_used"]][vi[["variable"]] == "z"])
}

cat(sprintf("Treatment effect on a continuous outcome, residual sd 1, n = 800, p = %d.\n", P))
cat(sprintf("Mean of %d reps. `lo=0` is the share of reps whose lower bound is exactly zero.\n\n", REPS))
cat(sprintf("%6s %-10s %7s %7s %7s %7s %7s %9s\n",
            "tau", "sparsity", "est", "bias", "atom", "lo=0", "covers", "prop_used"))

for (tau in TAUS) {
  for (lv in LEVELS) {
    m <- rowMeans(vapply(seq_len(REPS), function(r) one(tau, lv, r), numeric(6)))
    cat(sprintf("%6.2f %-10s %7.3f %+7.3f %7.3f %7.2f %7.2f %9.3f\n",
                tau, lv, m[["est"]], m[["est"]] - tau, m[["atom"]],
                m[["lo_is_zero"]], m[["covers"]], m[["prop_used"]]))
  }
  cat("\n")
}
