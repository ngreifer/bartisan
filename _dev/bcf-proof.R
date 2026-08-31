# Does the sampler already support a varying-coefficient (Bayesian causal
# forest) model? Nothing in the package exposes one, but `custom_family()` can
# express the likelihood, so this asks the question without adding a feature.
#
# The model is
#
#   g(mu_i) = m(x_i) + z_i * tau(x_i)
#
# with `m` and `tau` each a forest. That is Hahn, Murray and Carvalho (2020):
# the point is that `tau` gets its own forest and its own prior, rather than
# being whatever difference a single forest with `z` among its predictors
# happens to produce.
#
# THE DEVICE. A custom likelihood is called once per leaf with the rows
# reaching that leaf, and it is handed `y` for those rows and nothing else. It
# is not told which rows they are, so a treatment vector held outside cannot be
# lined up with them. The way around it here is to carry `z` inside the
# response and unpack it on arrival. For a binary outcome that is exact:
# `y + 2 * z` takes four integer values and both parts come back by division.
# For a continuous outcome it is a shift large enough to separate the two,
# which is exact to rounding.
#
# This is a device for answering the question, not a proposed interface. What
# it establishes is whether the sampler moves the two forests correctly when
# one of them is scaled by a covariate. See `_dev/bcf-interfaces.md`.

library(bartisan)

# ---- data ---------------------------------------------------------------

sim <- function(n, p = 5, binary_outcome = FALSE) {
  x <- matrix(rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))

  # Prognostic and treatment-effect surfaces. `tau` varies with x3 alone, so
  # recovery is checked against something no part of the fit was told.
  m <- 1 + 2 * x[, 1] - x[, 2] + 0.5 * x[, 1] * x[, 2]
  tau <- 1 + 0.8 * x[, 3]

  # Confounded assignment: the same covariates drive treatment and outcome.
  ps <- plogis(0.6 * x[, 1] - 0.4 * x[, 2])
  z <- rbinom(n, 1, ps)

  eta <- m + z * tau

  y <- if (binary_outcome) {
    rbinom(n, 1, plogis(eta - mean(eta)))
  } else {
    eta + rnorm(n)
  }

  data.frame(y = y, z = z, x, m = m, tau = tau, ps = ps)
}

# ---- the two families ---------------------------------------------------

# Gaussian. The scale is drawn as a nuisance parameter, which forces the
# derivatives to be differenced: `custom_family()` calls the derivative
# function as `f(y, eta, h)` with no `aux`, so an analytic score cannot see
# sigma. Correct either way; three R calls per leaf visit rather than one.
SHIFT <- 1000

bcf_gaussian <- custom_family(
  logdens = function(y, eta, aux) {
    z <- round(y / SHIFT)
    dnorm(y - SHIFT * z, eta[, 1] + z * eta[, 2], exp(aux[1]), log = TRUE)
  },
  num_predictors = 2,
  aux_names = "log_sigma",
  aux_start = 0,
  name = "bcf_gaussian"
)

# Binomial. No nuisance parameter, so the analytic derivatives apply, and the
# chain rule is the whole content of them: d mu / d eta_1 = 1, d mu / d eta_2 = z.
bcf_binomial <- custom_family(
  logdens = function(y, eta) {
    z <- y %/% 2
    stats::dbinom(y %% 2, 1L, plogis(eta[, 1] + z * eta[, 2]), log = TRUE)
  },
  derivatives = function(y, eta, h) {
    z <- y %/% 2
    prob <- plogis(eta[, 1] + z * eta[, 2])
    slope <- if (h == 1L) 1 else z
    list(score = (y %% 2 - prob) * slope,
         info  = prob * (1 - prob) * slope^2)
  },
  num_predictors = 2,
  name = "bcf_binomial"
)

# ---- fitting ------------------------------------------------------------

# Fewer trees for the treatment-effect forest than the prognostic one, which is
# the usual BCF choice: the effect surface is the simpler of the two, and the
# smaller forest shrinks it harder.
ctrl <- bartisan_control(num_trees = c(50, 25), num_burn = 750, num_draws = 750)

x_form <- ~ x1 + x2 + x3 + x4 + x5

fit_bcf <- function(d, family, packed) {
  d$packed <- packed
  bartisan(stats::update(x_form, packed ~ .), data = d,
           family = family, control = ctrl, chains = 2)
}

# The comparison: one forest with `z` among its predictors. This is what a user
# would do today, and the difference between the two forests it implies is the
# thing BCF replaces with an explicit one.
fit_s <- function(d, family) {
  bartisan(stats::update(x_form, y ~ z + .), data = d, family = family,
           control = bartisan_control(num_burn = 750, num_draws = 750),
           chains = 2)
}

rhat_detail <- function(fit) {
  r <- fit$rhat
  r <- r[order(-r$rhat), ]
  utils::head(r, 4)
}

report <- function(label, tau_draws, d) {
  tau_hat <- colMeans(tau_draws)
  ate <- rowMeans(tau_draws)

  cat(sprintf("%-22s cor %5.3f  rmse %5.3f  ATE %6.3f [%6.3f, %6.3f]  truth %5.3f\n",
              label,
              cor(tau_hat, d$tau),
              sqrt(mean((tau_hat - d$tau)^2)),
              mean(ate), quantile(ate, 0.025), quantile(ate, 0.975),
              mean(d$tau)))
  invisible(tau_hat)
}

# ---- run ----------------------------------------------------------------

set.seed(20260831)
dg <- sim(1000)

set.seed(1)
g_bcf <- fit_bcf(dg, bcf_gaussian, dg$y + SHIFT * dg$z)

cat("\n== Gaussian ==\n")
print(rhat_detail(g_bcf))
cat("sigma:", round(exp(mean(g_bcf$aux[, "log_sigma"])), 3), "(truth 1)\n\n")
report("BCF", g_bcf$eta[[2]], dg)

# The s-learner has no tau forest to read, so its effect comes from the
# contrast of its own predictions under z = 1 and z = 0.
set.seed(1)
g_s <- fit_s(dg, gaussian())
d1 <- transform(dg, z = 1)
d0 <- transform(dg, z = 0)
s_draws <- predict(g_s, newdata = d1, type = "link", summary = FALSE, draws = TRUE) -
  predict(g_s, newdata = d0, type = "link", summary = FALSE, draws = TRUE)
report("one forest with z", s_draws, dg)

# Hahn, Murray and Carvalho's own remedy for the two surfaces being weakly
# separated is to put an estimate of the propensity score into the prognostic
# forest. Nothing here can give a covariate to one forest and not the other, so
# this gives it to both, which is the closest a user could come today.
set.seed(1)
dg$ps_hat <- predict(bartisan(z ~ x1 + x2 + x3 + x4 + x5, data = dg,
                              family = binomial(), chains = 2,
                              control = bartisan_control(num_burn = 500,
                                                         num_draws = 500)))

set.seed(1)
dg$packed <- dg$y + SHIFT * dg$z
g_ps <- bartisan(packed ~ x1 + x2 + x3 + x4 + x5 + ps_hat, data = dg,
                 family = bcf_gaussian, control = ctrl, chains = 2)
cat("\nworst rhat, BCF + propensity score:",
    round(max(g_ps$rhat$rhat, na.rm = TRUE), 3), "\n")
report("BCF + ps", g_ps$eta[[2]], dg)

set.seed(20260831)
db <- sim(1500, binary_outcome = TRUE)

set.seed(1)
b_bcf <- fit_bcf(db, bcf_binomial, db$y + 2L * db$z)

cat("\n== Binomial (effects on the logit scale) ==\n")
print(rhat_detail(b_bcf))
cat("\n")
report("BCF", b_bcf$eta[[2]], db)

set.seed(1)
b_s <- fit_s(db, binomial())
d1 <- transform(db, z = 1)
d0 <- transform(db, z = 0)
s_draws <- predict(b_s, newdata = d1, type = "link", summary = FALSE, draws = TRUE) -
  predict(b_s, newdata = d0, type = "link", summary = FALSE, draws = TRUE)
report("one forest with z", s_draws, db)
