# What does bartisan's decision-rule prior do to the levels of a factor?
#
# Deshpande (2024) measures this with the prior co-clustering probability: for
# each pair of levels, how often a tree drawn from the prior puts them in the
# same leaf. In the infinite-tree limit that matrix is proportional to the prior
# covariance of the level means, so it is what determines how much the fit pools
# one level towards another.
#
# This simulates the tree prior directly rather than fitting anything, so it is
# a statement about the prior alone. Four schemes:
#
#   onehot   what bartisan does now. A factor is one group in the sparsity
#            prior, then one of its K indicator columns is drawn uniformly and
#            the rule is a threshold on that column, which peels that one level
#            off. Columns already used up the path are still drawn, and the
#            split is then degenerate, so that is modelled too.
#   subset   flexBART. Each available level goes left with probability 1/2,
#            conditioned on both children being non-empty.
#   order1   shuffle the level order once, then treat the factor as numeric.
#            Rules are thresholds, so a rule cuts a contiguous block of that
#            one order.
#   ordertree the same, reshuffled for every tree, so the ensemble mixes orders
#            even though each tree sees only one.

GAMMA <- 0.95
BETA <- 2

# The branching process bartisan uses: a node at depth d has children with
# probability gamma * (1 + d)^-beta.
grow <- function(d) stats::runif(1) < GAMMA * (1 + d)^(-BETA)

# One tree, returning the partition of `levels` its leaves induce. `available`
# is the set of levels still reachable at this node; `rule` draws a split of it.
partition_tree <- function(levels, rule, depth = 0) {
  if (length(levels) <= 1L || !grow(depth)) {
    return(list(levels))
  }

  split <- rule(levels)

  # A degenerate rule sends everything one way. The node still has two children,
  # but one of them is empty, so the levels carry on down a single branch and the
  # split has bought nothing.
  if (length(split$left) == 0L || length(split$right) == 0L) {
    return(partition_tree(levels, rule, depth + 1L))
  }

  c(partition_tree(split$left, rule, depth + 1L),
    partition_tree(split$right, rule, depth + 1L))
}

# --- the four rules -------------------------------------------------------

# `all_levels` is what bartisan draws a column from: every indicator of the
# factor, whether or not this path has used it already.
make_onehot <- function(all_levels) {
  function(levels) {
    k <- sample(all_levels, 1L)
    list(left = setdiff(levels, k), right = intersect(levels, k))
  }
}

subset_rule <- function(levels) {
  repeat {
    go_left <- stats::runif(length(levels)) < 0.5
    if (any(go_left) && !all(go_left)) {
      return(list(left = levels[go_left], right = levels[!go_left]))
    }
  }
}

# A threshold on a fixed order cuts a contiguous block of it.
make_order_rule <- function(order) {
  function(levels) {
    ranked <- levels[order(match(levels, order))]
    at <- sample(seq_len(length(ranked) - 1L), 1L)
    list(left = ranked[seq_len(at)], right = ranked[-seq_len(at)])
  }
}

# --- measurement ----------------------------------------------------------

co_cluster <- function(K, rule_factory, draws = 20000) {
  out <- matrix(0, K, K)
  leaves <- numeric(draws)
  singles <- numeric(draws)

  for (b in seq_len(draws)) {
    parts <- partition_tree(seq_len(K), rule_factory())
    leaves[b] <- length(parts)
    singles[b] <- sum(lengths(parts) == 1L)
    for (p in parts) out[p, p] <- out[p, p] + 1
  }

  list(kernel = out / draws, leaves = mean(leaves), singles = mean(singles))
}

report <- function(K, draws = 20000) {
  set.seed(1)
  schemes <- list(
    onehot    = function() make_onehot(seq_len(K)),
    subset    = function() subset_rule,
    order1    = local({ ord <- sample(K); function() make_order_rule(ord) }),
    ordertree = function() make_order_rule(sample(K))
  )

  cat(sprintf("\n== K = %d levels, %d prior draws ==\n", K, draws))
  cat(sprintf("%-10s %8s %8s %10s %10s\n",
              "scheme", "leaves", "singles", "co-clust", "spread"))

  for (nm in names(schemes)) {
    r <- co_cluster(K, schemes[[nm]], draws)
    off <- r$kernel[upper.tri(r$kernel)]
    cat(sprintf("%-10s %8.2f %8.2f %10.3f %10.3f\n",
                nm, r$leaves, r$singles, mean(off), stats::sd(off)))
  }
}

# How many partitions of K levels each scheme can reach at all, against the
# Bell number. `onehot` reaches 2^K - K; a single fixed order reaches the
# interval partitions of it, 2^(K-1); subsets reach all of them.
bell <- function(K) {
  b <- numeric(K + 1); b[1] <- 1
  for (n in seq_len(K)) {
    b[n + 1] <- sum(choose(n - 1, seq_len(n) - 1) * b[seq_len(n)])
  }
  b[K + 1]
}

cat("Partitions of K levels each scheme can reach:\n")
cat(sprintf("%4s %12s %12s %12s %12s\n", "K", "onehot", "order", "subset/all", "Bell"))
for (K in c(3, 5, 10, 20)) {
  cat(sprintf("%4d %12.0f %12.0f %12s %12.0f\n",
              K, 2^K - K, 2^(K - 1), "all", bell(K)))
}

for (K in c(5L, 10L)) report(K)
