# A factor with many levels: as a predictor, or as a random intercept?
#
# The two do different things, and subset rules change the comparison. A random
# intercept shrinks every level towards one common mean under one variance, which
# is the right prior when the level effects really are exchangeable draws from a
# distribution. A subset rule pools levels into *groups* that share a leaf, which
# is the right prior when the level effects take a few distinct values. Under
# one-hot the fixed route could not express the second at all, so the comparison
# was between a random intercept and a partition into singletons.
#
# Two truths, therefore, and the interesting question is whether each route wins
# on the truth that matches its prior:
#
#   iid       level effects are independent normal draws
#   clustered level effects take one of four values, five levels each

library(bartisan)

K <- 20L
REPS <- 10L
BURN <- 500L
DRAWS <- 500L
TREES <- 50L

# The level effects belong to the replicate, not to the dataset: train and test
# have to share them or no method can recover anything. Drawing them inside
# `sim()` from a seed that included the dataset gave the training and test sets
# different level effects, and every method then scored an RMSE larger than the
# standard deviation of the truth -- which is what gave that away.
level_effects <- function(rep, truth) {
  if (truth == "clustered") {
    return(rep(c(-3, -1, 1, 3), each = K / 4L))
  }

  set.seed(9000 + rep)
  stats::rnorm(K, sd = 2)
}

sim <- function(n, rep, truth = c("iid", "clustered"), seed) {
  truth <- match.arg(truth)
  effects <- level_effects(rep, truth)
  set.seed(seed)

  level <- sample(K, n, TRUE)
  x1 <- stats::runif(n)
  mu <- effects[level] + 2 * x1

  data.frame(y = mu + stats::rnorm(n), truth = mu,
             g = factor(sprintf("L%02d", level),
                        levels = sprintf("L%02d", seq_len(K))),
             x1 = x1)
}

rmse <- function(pred, test) sqrt(mean((pred - test$truth)^2))

fit_one <- function(form, tr, te, ...) {
  set.seed(7)
  f <- bartisan(form, data = tr, family = gaussian(),
                control = bartisan_control(num_trees = TREES, num_burn = BURN,
                                           num_draws = DRAWS, verbose = FALSE,
                                           ...))
  rmse(as.numeric(predict(f, newdata = te)), te)
}

run <- function(n, truth, rep) {
  tr <- sim(n, rep, truth, seed = 4000 + rep)
  te <- sim(2000L, rep, truth, seed = 6000 + rep)

  c(`fixed, subset` = fit_one(y ~ g + x1, tr, te),
    `fixed, onehot` = fit_one(y ~ g + x1, tr, te, categorical = "onehot"),
    `random intercept` = fit_one(y ~ x1 + (1 | g), tr, te),
    `both` = fit_one(y ~ g + x1 + (1 | g), tr, te))
}

cat(sprintf("%d levels, 50 trees, %d warmup, %d draws, %d replicates.\n",
            K, BURN, DRAWS, REPS))
cat("RMSE against the true mean function. Every method sees the same data within\n")
cat("a replicate, so the comparison is paired and the standard error is of the\n")
cat("paired difference from the best method, which is what the gaps have to beat.\n")

for (truth in c("iid", "clustered")) {
  for (n in c(200L, 1000L)) {
    each <- vapply(seq_len(REPS), function(r) run(n, truth, r),
                   numeric(4L))
    means <- rowMeans(each)
    best <- which.min(means)

    cat(sprintf("\n%-9s truth, n = %d (%.0f per level)\n", truth, n, n / K))
    for (i in seq_along(means)) {
      diff <- each[i, ] - each[best, ]
      se <- stats::sd(diff) / sqrt(REPS)
      cat(sprintf("  %-18s %.4f  %+.4f %s\n", rownames(each)[i], means[i],
                  means[i] - means[best],
                  if (i == best) "" else sprintf("(se %.4f)", se)))
    }
  }
}
