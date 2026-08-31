# Does the engine's categorical rule reach the partitions the prior is meant to?
# Read straight off fitted trees: a partition of the levels is what the leaves
# induce, and with no signal in the data the fit is drawing close to the prior.

library(bartisan)

set.seed(1)
n <- 800
K <- 5
g <- factor(sample(letters[seq_len(K)], n, TRUE))
d <- data.frame(y = stats::rnorm(n), g = g)

# The partition each draw's forest induces over the levels, read off predictions:
# two levels are in the same cell of a single tree's partition exactly when that
# tree predicts the same value for them. Using a one-tree forest makes the
# forest's partition the tree's.
fit <- bartisan(y ~ g, data = d, family = gaussian(),
                control = bartisan_control(num_trees = 1, num_burn = 500,
                                           num_draws = 4000, gate = "hard",
                                           sparsity = FALSE))

one <- d[match(levels(g), d$g), , drop = FALSE]
draws <- predict(fit, newdata = one, type = "link", summary = FALSE,
                 draws = TRUE)

# Canonical label of the partition each draw induces.
labels <- apply(draws, 1L, function(v) paste(match(v, unique(v)), collapse = ""))

seen <- length(unique(labels))
co <- matrix(0, K, K)
for (b in seq_len(nrow(draws))) {
  v <- draws[b, ]
  co <- co + outer(v, v, "==")
}
co <- co / nrow(draws)

bell <- function(K) {
  b <- numeric(K + 1); b[1] <- 1
  for (m in seq_len(K)) b[m + 1] <- sum(choose(m - 1, seq_len(m) - 1) * b[seq_len(m)])
  b[K + 1]
}

cat(sprintf("K = %d, %d posterior draws of a one-tree forest on pure noise\n",
            K, nrow(draws)))
cat(sprintf("distinct partitions visited: %d  (one-hot can reach %d, all = %d)\n",
            seen, 2^K - K, bell(K)))
cat(sprintf("mean off-diagonal co-clustering: %.3f\n",
            mean(co[upper.tri(co)])))

# The partitions one-hot cannot express are those with two or more cells of size
# two or more. Visiting any of them is the thing that was impossible before.
non_onehot <- vapply(unique(labels), function(l) {
  sizes <- table(strsplit(l, "")[[1L]])
  sum(sizes >= 2) >= 2
}, logical(1L))

cat(sprintf("of those, %d are partitions one-hot encoding cannot form at all\n",
            sum(non_onehot)))
cat(sprintf("share of draws in such a partition: %.3f\n",
            mean(labels %in% unique(labels)[non_onehot])))
