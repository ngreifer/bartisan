# The BCF proof showed rhat around 1.4 on the prognostic forest, with the
# treatment-effect forest near 1.10. Is that the varying coefficient, or is it
# the custom-family machinery it had to be written in? The comparison in
# `bcf-proof.R` was not fair: it put a hand-written likelihood whose scale is a
# nuisance parameter drawn by the generic Laplace-plus-Metropolis step against
# `gaussian()`, whose scale has a Gibbs step.
#
# So: same data, same custom-family route, with and without the varying
# coefficient.

library(bartisan)

set.seed(20260831)
n <- 1000
p <- 5
x <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
m <- 1 + 2 * x[, 1] - x[, 2] + 0.5 * x[, 1] * x[, 2]
tau <- 1 + 0.8 * x[, 3]
z <- rbinom(n, 1, plogis(0.6 * x[, 1] - 0.4 * x[, 2]))
y <- m + z * tau + rnorm(n)
d <- data.frame(y = y, z = z, x)

SHIFT <- 1000
d$packed <- y + SHIFT * z

ctrl <- bartisan_control(num_burn = 750, num_draws = 750)
f_x <- packed ~ x1 + x2 + x3 + x4 + x5

# 1. Hand-written Gaussian, one forest, no varying coefficient. Same nuisance
#    parameter, same generic sampler for it. This is the control.
plain <- custom_family(
  function(y, eta, aux) {
    dnorm(y - SHIFT * round(y / SHIFT), eta[, 1], exp(aux[1]), log = TRUE)
  },
  num_predictors = 1, aux_names = "log_sigma", name = "plain")

# 2. The same thing with the varying coefficient added.
vc <- custom_family(
  function(y, eta, aux) {
    zz <- round(y / SHIFT)
    dnorm(y - SHIFT * zz, eta[, 1] + zz * eta[, 2], exp(aux[1]), log = TRUE)
  },
  num_predictors = 2, aux_names = "log_sigma", name = "vc")

# 3. The compiled Gaussian, for reference on how this data mixes at all.
show <- function(label, fit) {
  r <- fit$rhat[order(-fit$rhat$rhat), ]
  cat("\n--", label, "--\n")
  print(utils::head(r, 3), row.names = FALSE)
}

set.seed(1)
show("custom, 1 forest, no varying coefficient",
     bartisan(f_x, data = d, family = plain, control = ctrl, chains = 4))

set.seed(1)
show("custom, 2 forests, varying coefficient on z",
     bartisan(f_x, data = d, family = vc,
              control = bartisan_control(num_burn = 750, num_draws = 750,
                                         num_trees = c(50, 25)),
              chains = 4))

set.seed(1)
show("compiled gaussian(), z as an ordinary predictor",
     bartisan(y ~ z + x1 + x2 + x3 + x4 + x5, data = d, family = gaussian(),
              control = ctrl, chains = 4))
