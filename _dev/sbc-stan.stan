// A sum of soft trees with known structure, Poisson response: the model of
// `_dev/sbc.R`'s soft-rule Poisson arm with the tree structures given as data,
// so that a general-purpose sampler can fit the leaf values and the per-tree
// bandwidths. See `_dev/sbc-stan.R`.
data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;
  array[N] int<lower=0> y;
  vector[N] log_offset;
  int<lower=1> T;                       // trees
  int<lower=1> K;                       // leaves, over all trees
  int<lower=1> D;                       // longest root-to-leaf path
  array[K] int<lower=1, upper=T> leaf_tree;
  array[K] int<lower=0, upper=D> depth; // path length of each leaf
  array[K, D] int<lower=0, upper=P> split_var;
  array[K, D] real split_val;
  array[K, D] int<lower=0, upper=1> goes_left;  // 1 if the path goes left there
  real<lower=0> sigma_mu;
  real<lower=0> bandwidth_mean;
  real<lower=0> half_width;
  int<lower=1, upper=N> A;
  int<lower=1, upper=N> B;
}
transformed data {
  // Each step's distance from the cutpoint, which does not depend on any
  // parameter, so the gates below cost one scaling and a clamp.
  array[K, D] vector[N] dist;
  for (k in 1:K) {
    for (s in 1:D) {
      dist[k, s] = s <= depth[k] ? split_val[k, s] - X[, split_var[k, s]]
                                 : rep_vector(0, N);
    }
  }
}
parameters {
  vector[K] mu;
  vector<lower=0>[T] bandwidth;
}
transformed parameters {
  vector[N] eta = log_offset;
  for (k in 1:K) {
    vector[N] w = rep_vector(1, N);
    real half = bandwidth[leaf_tree[k]] * half_width;
    for (s in 1:depth[k]) {
      vector[N] t = fmin(fmax(0.5 + (0.5 / half) * dist[k, s], 0), 1);
      vector[N] p = t .* t .* (3 - 2 * t);
      w = goes_left[k, s] == 1 ? w .* p : w .* (1 - p);
    }
    eta += w * mu[k];
  }
}
model {
  mu ~ normal(0, sigma_mu);
  bandwidth ~ exponential(1 / bandwidth_mean);
  y ~ poisson_log(eta);
}
generated quantities {
  real contrast = eta[A] - eta[B];
}
