# The coding of a moderated predictor, four ways.
#
# `center` decides what the control function sits at, and for a fixed choice it
# is a reparameterization of the control function alone -- the coefficient is
# untouched. `"estimate"` is different in kind: it draws the coding rather than
# fixing it (Hahn, Murray and Carvalho 2020, sec. 5.3), which makes the model
# invariant to how the levels are labelled, and for more than two levels it also
# imposes a restriction, since one shared forest is scaled by a coefficient per
# level rather than each level having a forest of its own.
#
# Paired within replicate, so the standard error is of the paired difference.

library(bartisan)

REPS <- 12L
BURN <- 500L
DRAWS <- 500L

fit_one <- function(form, d, n_trees = 50L, ...) {
  set.seed(7)
  bartisan(form, data = d, family = gaussian(), num_trees = n_trees,
           num_burn = BURN, num_draws = DRAWS, verbose = FALSE, ...)
}

report <- function(title, each) {
  means <- rowMeans(each)
  best <- which.min(means)
  cat(sprintf("\n%s\n", title))
  cat(sprintf("  %-24s %8s %9s %9s\n", "coding", "RMSE", "vs best", "se"))
  for (i in seq_along(means)) {
    d <- each[i, ] - each[best, ]
    cat(sprintf("  %-24s %8.4f %+9.4f %9s\n", rownames(each)[i], means[i],
                means[i] - means[best],
                if (i == best) "" else sprintf("%.4f", stats::sd(d) / sqrt(ncol(each)))))
  }
}

# ---- binary ---------------------------------------------------------------

binary_rep <- function(rep) {
  set.seed(200 + rep)
  n <- 800
  x <- matrix(rnorm(n * 4), n, 4, dimnames = list(NULL, paste0("x", 1:4)))
  tau <- 1 + 0.8 * x[, 3]
  z <- rbinom(n, 1, plogis(0.6 * x[, 1]))
  d <- data.frame(y = 2 * x[, 1] - x[, 2] + z * tau + rnorm(n), z = z, x)

  err <- function(f, effect) sqrt(mean((effect - tau)^2))

  c(mean = err(NULL, coef(fit_one(y ~ x1 + x2 + x3 + x4 + vc(z, center = "mean"), d))[, 1L]),
    zero = err(NULL, coef(fit_one(y ~ x1 + x2 + x3 + x4 + vc(z, center = "zero"), d))[, 1L]),
    mid = err(NULL, coef(fit_one(y ~ x1 + x2 + x3 + x4 + vc(z, center = "mid"), d))[, 1L]),
    estimate = err(NULL, coef(fit_one(y ~ x1 + x2 + x3 + x4 + vc(z, center = "estimate"), d))[, 1L]))
}

report("Binary moderated predictor, RMSE of the effect function",
       vapply(seq_len(REPS), binary_rep, numeric(4L)))

# ---- binary: the invariance the estimated coding exists to buy -------------

cat("\nInvariance to how the treatment is labelled.\n")
cat("The same data with the treatment coded 0/1 and 1/0 is the same model, so\n")
cat("the effect should come back the same size with the sign flipped. Run at two\n")
cat("signal strengths, because the codings differ only in the prior and a clear\n")
cat("signal swamps that difference.\n")

invariance <- function(n, effect_size, noise, rep) {
  set.seed(900 + rep)
  x <- matrix(rnorm(n * 4), n, 4, dimnames = list(NULL, paste0("x", 1:4)))
  tau <- effect_size * (1 + 0.8 * x[, 3])
  z <- rbinom(n, 1, plogis(0.6 * x[, 1]))
  d <- data.frame(y = 2 * x[, 1] - x[, 2] + z * tau + rnorm(n, 0, noise),
                  z = z, x)
  flipped <- transform(d, z = 1 - z)

  vapply(c("zero", "mean", "estimate"), function(center) {
    form <- stats::reformulate(
      c("x1", "x2", "x3", "x4", sprintf('vc(z, center = "%s")', center)),
      response = quote(y))
    mean(coef(fit_one(form, d))[, 1L]) + mean(coef(fit_one(form, flipped))[, 1L])
  }, numeric(1L))
}

for (design in list(list(lab = "n = 800, effect 1.0, noise 1", n = 800L, e = 1, s = 1),
                    list(lab = "n = 150, effect 0.2, noise 2", n = 150L, e = 0.2, s = 2))) {
  each <- vapply(seq_len(REPS), function(r) {
    invariance(design[["n"]], design[["e"]], design[["s"]], r)
  }, numeric(3L))
  cat(sprintf("\n  %s -- mean |ATE(0/1) + ATE(1/0)|, which is 0 if the coding does not matter\n",
              design[["lab"]]))
  for (i in seq_len(nrow(each))) {
    cat(sprintf("    %-9s %7.4f  (se %.4f)\n", rownames(each)[i],
                mean(abs(each[i, ])), stats::sd(abs(each[i, ])) / sqrt(ncol(each))))
  }
}

# ---- categorical ----------------------------------------------------------

categorical_rep <- function(rep, shared_shape) {
  set.seed(400 + rep)
  n <- 900
  x <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, paste0("x", 1:3)))
  g <- factor(sample(c("a", "b", "c"), n, TRUE))

  # Rank one: every level's effect is the same shape times a scalar, which is
  # exactly what a drawn coding assumes. Otherwise each level has its own shape,
  # which it cannot represent.
  shape <- 1 + 0.8 * x[, 3]
  effect <- if (shared_shape) {
    c(a = 0, b = 1, c = 2)[as.character(g)] * shape
  } else {
    cbind(0, 1 + 0.8 * x[, 3], 2 - 0.9 * x[, 1])[
      cbind(seq_len(n), as.integer(g))]
  }

  d <- data.frame(y = 2 * x[, 1] + effect + rnorm(n), g = g, x)

  err <- function(cf) {
    # Compare the fitted level contrasts against the truth's.
    truth <- if (shared_shape) {
      cbind(1 * shape, 2 * shape)
    } else {
      cbind(1 + 0.8 * x[, 3], 2 - 0.9 * x[, 1])
    }
    sqrt(mean((cf - truth)^2))
  }

  symmetric <- coef(fit_one(y ~ x1 + x2 + x3 + vc(g), d))
  symmetric <- symmetric[, c("gb", "gc")] - symmetric[, "ga"]

  c(`symmetric, one forest per level` = err(symmetric),
    `estimate, one shared forest` =
      err(coef(fit_one(y ~ x1 + x2 + x3 + vc(g, center = "estimate"), d))))
}

report("Categorical, truth is rank one (every level the same shape)",
       vapply(seq_len(REPS), function(r) categorical_rep(r, TRUE), numeric(2L)))
report("Categorical, truth is not rank one (each level its own shape)",
       vapply(seq_len(REPS), function(r) categorical_rep(r, FALSE), numeric(2L)))
