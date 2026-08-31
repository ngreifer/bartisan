# Does the leaf-scale mixing problem belong to bartisan, or to the parameter?
#
# `sigma_mu` shows rhat between 1.2 and 2.0 with ess_bulk under 15 out of 3000
# draws on every design tried, while the fitted surface mixes fine. Before
# spending more on it, check the other BART packages. Each names the same
# quantity differently and only some of them draw it at all:
#
#   bartisan   sigma_mu, half-Cauchy prior, independence Metropolis step
#   dbarts     k, the leaf sd is 0.5 / (k sqrt(m)), optional chi hyperprior
#   stochtree  sigma2_leaf, inverse-gamma prior, conjugate Gibbs step
#   BART       k = 2, fixed, never drawn
#
# The samplers differ, so if all of the ones that draw it mix badly then the
# problem is the parameter and not the step. stochtree is the decisive case: its
# prior is conjugate, so its draw is exact.

library(bartisan)

CHAINS <- 4
BURN <- 750
SAVE <- 750

set.seed(1)
n <- 600
p <- 5
x <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
fx <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + 10 * x[, 4]
y <- fx + rnorm(n)
d <- data.frame(y = y, x)

# rhat and ess on a draws-by-chains matrix, using this package's own functions so
# every row of the table is computed the same way.
report <- function(label, per_chain, note = "") {
  if (is.null(per_chain)) {
    cat(sprintf("%-34s %s\n", label, note))
    return(invisible(NULL))
  }
  cat(sprintf("%-34s rhat %5.3f  ess %6.1f  %s\n", label,
              bartisan:::split_rhat(per_chain),
              bartisan:::ess_bulk(per_chain), note))
}

# ---- bartisan -----------------------------------------------------------

set.seed(2)
fit <- bartisan(y ~ ., data = d, family = gaussian(), chains = CHAINS,
                control = bartisan_control(num_burn = BURN, num_draws = SAVE))
report("bartisan sigma_mu", matrix(fit$sigma_mu[, 1], ncol = CHAINS),
       "(half-Cauchy, independence MH)")
report("bartisan sd(eta), for contrast",
       matrix(apply(fit$eta[[1]], 1, stats::sd), ncol = CHAINS))

# ---- dbarts -------------------------------------------------------------

# `k` drawn under its chi hyperprior. dbarts runs one chain per call here so the
# chains are directly comparable to the others.
dk <- vapply(seq_len(CHAINS), function(i) {
  set.seed(100 + i)
  # The docs call this `chi(degreesOfFreedom, scale)`; the generator is not
  # exported, so the class is constructed directly.
  kprior <- methods::new(
    methods::getClass("dbartsChiHyperprior", where = asNamespace("dbarts")),
    degreesOfFreedom = 1.25, scale = Inf)
  f <- dbarts::bart2(y ~ ., data = d, k = kprior,
                     n.burn = BURN, n.samples = SAVE, n.chains = 1L,
                     n.trees = 50L, verbose = FALSE, keepTrees = FALSE)
  as.numeric(f$k)
}, numeric(SAVE))
report("dbarts k", dk, "(chi hyperprior, slice)")

# ---- stochtree ----------------------------------------------------------

# The decisive one: an inverse-gamma prior, so the leaf scale is drawn exactly
# by a Gibbs step rather than by a Metropolis step.
sk <- vapply(seq_len(CHAINS), function(i) {
  set.seed(200 + i)
  f <- stochtree::bart(
    X_train = as.data.frame(x), y_train = y,
    num_gfr = 0L, num_burnin = BURN, num_mcmc = SAVE,
    general_params = list(num_chains = 1L, verbose = FALSE),
    mean_forest_params = list(num_trees = 50L, sample_sigma2_leaf = TRUE))
  as.numeric(f$sigma2_leaf_samples)
}, numeric(SAVE))
report("stochtree sigma2_leaf", sk, "(inverse-gamma, Gibbs)")

# ---- BART ---------------------------------------------------------------

report("BART k", NULL, "not drawn: k = 2 is fixed, so no diagnostic exists")
