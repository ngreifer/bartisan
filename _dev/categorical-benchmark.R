# Subset rules against one-hot, and against flexBART, on the case they are for:
# levels that share a mean in groups, with the data thin enough per level that
# pooling is what buys accuracy.
#
# One-hot can only express a partition with at most one cell holding more than
# one level, so to pool a group of five levels a single tree has to isolate each
# of them. Subset rules put the group in one cell. The axis that should decide
# whether that matters is observations per level: with plenty, each level's mean
# is estimated well on its own and there is nothing to pool.

library(bartisan)

CLUSTERS <- 4L
PER_CLUSTER <- 5L
K <- CLUSTERS * PER_CLUSTER
REPS <- 5L
BURN <- 500L
DRAWS <- 500L
TREES <- 50L

sim <- function(n, rep) {
  set.seed(1000 + rep)
  level <- sample(K, n, TRUE)
  cluster <- (level - 1L) %/% PER_CLUSTER + 1L
  x1 <- stats::runif(n)
  truth <- c(-3, -1, 1, 3)[cluster] + 2 * x1

  data.frame(y = truth + stats::rnorm(n), truth = truth,
             g = factor(sprintf("L%02d", level), levels = sprintf("L%02d", 1:K)),
             x1 = x1)
}

rmse <- function(pred, test) sqrt(mean((pred - test$truth)^2))

run <- function(n_train, rep) {
  tr <- sim(n_train, rep)
  te <- sim(2000L, rep + 500L)
  out <- list()

  # Both rule kinds under hard decision rules, which is what flexBART has, so
  # the comparison is of the categorical rule and not of soft against hard. The
  # soft default is reported alongside, since it is what a user actually gets.
  for (gate in c("hard", "smoothstep")) {
    for (mode in c("subset", "onehot")) {
      set.seed(7)
      t0 <- proc.time()[["elapsed"]]
      f <- bartisan(y ~ g + x1, data = tr, family = gaussian(),
                    control = bartisan_control(num_trees = TREES,
                                               num_burn = BURN,
                                               num_draws = DRAWS,
                                               categorical = mode,
                                               gate = gate,
                                               verbose = FALSE))
      secs <- proc.time()[["elapsed"]] - t0
      label <- sprintf("%s, %s", mode,
                       if (gate == "hard") "hard" else "soft")
      out[[label]] <-
        c(rmse = rmse(as.numeric(predict(f, newdata = te)), te), secs = secs)
    }
  }

  # One chain and the same iteration counts, so the comparison is of the rules
  # rather than of how much sampling each package does by default. `yhat.test`
  # is draws by observations; `yhat.test.mean` is the posterior mean already.
  set.seed(7)
  t0 <- proc.time()[["elapsed"]]
  fb <- flexBART::flexBART(
    y ~ bart(g + x1),
    train_data = tr[, c("y", "g", "x1")],
    test_data = te[, c("g", "x1")],
    M = TREES, n.chains = 1L, nd = DRAWS, burn = BURN, verbose = FALSE)
  secs <- proc.time()[["elapsed"]] - t0
  out[["flexBART"]] <- c(rmse = rmse(as.numeric(fb$yhat.test.mean), te),
                         secs = secs)

  do.call(rbind, out)
}

cat(sprintf("%d levels in %d clusters of %d. RMSE against the true mean function,\n",
            K, CLUSTERS, PER_CLUSTER))
cat(sprintf("%d trees, %d warmup, %d draws, one chain, mean of %d replicates.\n",
            TREES, BURN, DRAWS, REPS))

for (n_train in c(200L, 500L, 2000L)) {
  acc <- Reduce(`+`, lapply(seq_len(REPS), function(r) run(n_train, r))) / REPS
  cat(sprintf("\nn = %d  (%.0f observations per level)\n",
              n_train, n_train / K))
  cat(sprintf("%-18s %8s %9s\n", "rule", "RMSE", "seconds"))
  for (nm in rownames(acc)) {
    cat(sprintf("%-18s %8.4f %9.2f\n", nm, acc[nm, "rmse"], acc[nm, "secs"]))
  }
}
