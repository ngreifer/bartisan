// The leaves-only model of `_dev/sbc.R` with `SBC_STAN`: each tree's structure
// and bandwidth given, so the predictor is linear in the leaf values,
// eta = base + W mu (base the offset), with W the gate weight of each observation in each leaf.
// A Poisson or binomial (logit) response. See `_dev/TASKS.md`, 2026-11-05.
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] W;
  array[N] int<lower=0> y;
  array[N] int<lower=0> trials;
  vector[N] base;
  real<lower=0> sigma_mu;
  int<lower=1, upper=2> family;   // 1 Poisson, 2 binomial
  int<lower=1, upper=N> A;
  int<lower=1, upper=N> B;
}
parameters {
  vector[K] mu;
}
model {
  mu ~ normal(0, sigma_mu);
  if (family == 1) {
    y ~ poisson_log_glm(W, base, mu);
  } else {
    y ~ binomial_logit(trials, base + W * mu);
  }
}
generated quantities {
  real contrast = dot_product(W[A] - W[B], mu);
  real level = base[A] + dot_product(W[A], mu);
}
