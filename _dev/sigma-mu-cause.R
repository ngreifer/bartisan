# The ramp bug was real and fixing it changed nothing, so `sigma_mu` mixing has
# a different cause. Two candidates:
#
#   (a) the sampler. It is an independence Metropolis step whose proposal is the
#       conditional posterior of the precision under a flat prior, corrected by
#       the half-Cauchy term. If the correction is often decisive the chain
#       sticks.
#   (b) the parameterization. That proposal is Gamma(n/2 + 1, 2/sse) over the n
#       leaf parameters of the whole forest. With n in the hundreds it is very
#       tight, so `sigma_mu` is close to a deterministic readout of the root mean
#       square leaf value, and it can mix no faster than that aggregate does.
#
# (a) predicts a low acceptance rate. (b) predicts a high one, `sigma_mu` ess
# tracking the ess of a global aggregate of the fitted surface, and ess falling
# as trees are added, because more leaves make the proposal tighter and tie
# `sigma_mu` more closely to the forest.

library(bartisan)

set.seed(1)
n <- 600
x <- matrix(runif(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
d <- data.frame(y = 10 * sin(pi * x[, 1] * x[, 2]) + rnorm(n), x)

# `ess_bulk()` wants draws by chains, which is the shape both probes build.
ess <- function(m) bartisan:::ess_bulk(m)

probe <- function(trees, burn = 750, save = 750) {
  set.seed(2)
  fit <- bartisan(y ~ ., data = d, family = gaussian(), chains = 4,
                  control = bartisan_control(num_trees = trees, num_burn = burn,
                                             num_draws = save))
  sm <- fit$sigma_mu[, 1]

  # Acceptance rate: an independence proposal that is rejected leaves the value
  # untouched, so consecutive equal draws are rejections. Per chain, because the
  # draws of different chains are unrelated.
  per_chain <- matrix(sm, ncol = 4)
  accept <- mean(apply(per_chain, 2, function(v) mean(diff(v) != 0)))

  # A global aggregate of the fitted surface: its spread across observations,
  # one number per draw. Same kind of quantity as sigma_mu and no hyperparameter
  # sampler involved.
  spread <- matrix(apply(fit$eta[[1]], 1, stats::sd), ncol = 4)

  r <- fit$rhat
  cat(sprintf("trees %3d burn/save %4d  accept %.2f | sigma_mu rhat %5.3f ess %5.1f | sd(eta) ess %5.1f | eta rhat %5.3f\n",
              trees, burn, accept,
              max(r$rhat[startsWith(r$quantity, "sigma_mu")]),
              ess(per_chain), ess(spread),
              max(r$rhat[startsWith(r$quantity, "eta")])))
}

for (tr in c(10L, 50L, 200L)) probe(tr)
cat("\nLonger run, to separate slow mixing from not having converged:\n")
probe(50L, burn = 4000, save = 4000)

# Both candidates are refuted by the numbers above. Acceptance is 0.89 to 0.97,
# which rules out (a). And `sd(eta)` -- a global aggregate of the fitted surface
# -- has an ess in the thousands while `sigma_mu` has single digits, so
# `sigma_mu` is not inheriting the forest's mixing either, which rules out (b)
# as stated.
#
# The two informative facts are that adding trees makes `sigma_mu` worse while
# making everything else better, and that its ess does not grow with the length
# of the run (10.1 at 750 saved draws, 8.1 at 4000). A quantity whose ess does
# not grow with the run is not mixing slowly, it is not moving: each chain
# settles somewhere and stays.
#
# Candidate (c): an additive ensemble can trade leaf magnitude against tree
# structure without changing its sum. Many small leaves and few large ones give
# the same fitted surface, and `sigma_mu` is the root mean square leaf value, so
# it distinguishes configurations the likelihood cannot. Adding trees makes the
# proposal Gamma(n/2 + 1, 2/sse) tighter, which pins `sigma_mu` harder to
# whichever configuration its chain is in, which is why more trees is worse.
#
# The test: chains that grow more splits should show a smaller `sigma_mu`.

cat("\n-- (c): leaf magnitude against tree count, per chain --\n")

for (tr in c(50L, 200L)) {
  set.seed(2)
  fit <- bartisan(y ~ ., data = d, family = gaussian(), chains = 4,
                  control = bartisan_control(num_trees = tr, num_burn = 750,
                                             num_draws = 750))
  sm <- matrix(fit$sigma_mu[, 1], ncol = 4)
  splits <- matrix(rowSums(fit$counts$eta), ncol = 4)

  cat(sprintf("trees %3d  per-chain mean sigma_mu %s\n", tr,
              paste(sprintf("%.3f", colMeans(sm)), collapse = " ")))
  cat(sprintf("           per-chain mean splits   %s   cor %.2f\n",
              paste(sprintf("%.1f", colMeans(splits)), collapse = " "),
              stats::cor(colMeans(sm), colMeans(splits))))
}
