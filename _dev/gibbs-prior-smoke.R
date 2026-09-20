# Does the `sparsity = "gibbs"` prototype do what the Gibbs prior is supposed to?
#
# WHAT IS BEING TESTED
#
# Linero and Du (2023) claim the Gibbs prior makes fewer false variable
# selections than DART because it penalizes the number of predictors directly,
# where the Dirichlet prior reaches sparsity only indirectly and "is more
# tolerant of variables which have miniscule impact on the outcome". Their
# Table 1 at P = 7, sigma = 3 puts DART's precision at 0.79 against the Gibbs
# prior's 0.98, with recall 1.00 for both, so the whole difference is false
# positives.
#
# `_dev/varsel-check.R` measured the same failure here from the other side: the
# median probability model fires on every replicate of a pure null, 7.1 of 25
# predictors selected. If the prototype is wired up correctly, it should move
# that number and should not move recall.
#
# This is a smoke test, not the measurement that would decide anything. It is
# four replicates at a small size, and its job is to establish that the urn is
# reached, that it changes the answer in the direction the paper says, and that
# the guards fire. A real comparison against DART belongs in its own script with
# enough replicates to separate 0.79 from 0.98.
#
# WHAT THE OUTCOME MEASUREMENT IS
#
# Three designs, because the first one tried could not tell the two priors
# apart: Friedman at n = 400 with p = 20 gave both a clean 5 of 5 with no false
# positives, so there was nothing for the Gibbs prior to improve. That is worth
# recording rather than quietly swapping out -- DART fails at *small* p, which
# is the paper's own finding, and on a null, which is where `_dev/varsel-check.R`
# measured it here.
#
#   sparse-20   Friedman, n = 400, p = 20, 5 real. The easy case, kept as the
#               check that neither prior loses a real predictor.
#   sparse-7    Friedman, n = 400, p = 7, 5 real. The paper's own worst case for
#               DART, where its precision reads 0.79 against the Gibbs 0.98.
#   null-15     n = 400, p = 15, nothing matters. Everything selected is false.
#
# Per fit, read off `variable_importance()` at the median probability model
# (`prop_used >= 0.5`):
#
#   noise      how many predictors with no signal are selected. The decisive
#              column, and the one the paper's precision gap lives in.
#   recovered  how many real ones it finds. The guard against a prototype that
#              is "sparse" only because it is broken.
#
# WHAT EACH OUTCOME WOULD MEAN
#
# Fewer noise predictors selected at equal recovery: the urn is doing the thing
# the paper describes and the prototype is worth measuring properly.
#
# Identical to DART: the urn is not being reached, or the counts it reads are
# stale, and the wiring is wrong rather than the prior being ineffective.
#
# Lower recovery: the prior is over-penalizing, which at zeta = 1 would be a bug
# rather than a property, since that is Linero and Du's own default.

library(bartisan)

REPS <- 6L
N <- 400L

# name -> (p, number of real predictors)
DESIGNS <- list("sparse-20" = c(p = 20L, real = 5L),
                "sparse-7"  = c(p =  7L, real = 5L),
                "null-15"   = c(p = 15L, real = 0L))

sim <- function(rep, p, real) {
  set.seed(4100 + 17 * rep + p)
  X <- matrix(stats::runif(N * p), N, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))

  truth <- if (real == 0L) {
    rep(0, N)
  } else {
    10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
      10 * X[, 4] + 5 * X[, 5]
  }

  data.frame(y = truth + stats::rnorm(N, 0, 3), X)
}

score <- function(fit, real) {
  vi <- variable_importance(fit)
  used <- vi[["prop_used"]]
  names(used) <- vi[["variable"]]
  picked <- names(used)[used >= 0.5]
  truth <- if (real == 0L) character(0) else paste0("x", seq_len(real))
  c(noise = sum(!picked %in% truth), recovered = sum(truth %in% picked))
}

ctrl <- function(sparsity) {
  bartisan_control(num_trees = 20L, num_burn = 200L, num_draws = 400L,
                   sparsity = sparsity)
}

rows <- list()

for (nm in names(DESIGNS)) {
  spec <- DESIGNS[[nm]]

  for (rep in seq_len(REPS)) {
    d <- sim(rep, spec[["p"]], spec[["real"]])

    for (sp in c("moderate", "gibbs")) {
      set.seed(10L * rep)
      fit <- bartisan(y ~ ., d, family = gaussian(), control = ctrl(sp))
      rows[[length(rows) + 1L]] <- data.frame(design = nm, rep = rep,
                                              sparsity = sp,
                                              t(score(fit, spec[["real"]])))
    }
  }
}

res <- do.call(rbind, rows)

cat("\n--- means over", REPS, "replicates ---\n")
print(aggregate(cbind(noise, recovered) ~ design + sparsity, res, mean),
      row.names = FALSE)

cat("\n--- paired difference in false positives (gibbs - dart) ---\n")
for (nm in names(DESIGNS)) {
  w <- res[res[["design"]] == nm, ]
  g <- w[w[["sparsity"]] == "gibbs", "noise"]
  m <- w[w[["sparsity"]] == "moderate", "noise"]
  cat(sprintf("%-10s %+6.2f  (dart %.2f, gibbs %.2f)\n",
              nm, mean(g - m), mean(m), mean(g)))
}

# The guards, which are cheap and are the other half of the wiring.
cat("\n--- guards ---\n")

expect_error <- function(label, expr) {
  ok <- inherits(try(force(expr), silent = TRUE), "try-error")
  cat(sprintf("%-28s %s\n", label, if (ok) "errors (good)" else "SILENT (bad)"))
}

expect_error("gibbs + split_prior",
             bartisan_control(sparsity = "gibbs",
                              split_prior = c(x1 = 3)))
expect_error("gibbs + share_sparsity",
             bartisan_control(sparsity = "gibbs", share_sparsity = TRUE))

cat("\n--- reproducibility ---\n")
d <- sim(1L, 20L, 5L)
set.seed(99L)
a <- bartisan(y ~ ., d, family = gaussian(), control = ctrl("gibbs"))
set.seed(99L)
b <- bartisan(y ~ ., d, family = gaussian(), control = ctrl("gibbs"))
cat("same seed, same fitted values:",
    isTRUE(all.equal(fitted(a), fitted(b))), "\n")
