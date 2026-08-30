#include "mcmc.h"

namespace bartisan {


namespace {

// Fisher scoring stops once the score falls below this many standard errors of
// the current fit, so the proposal is centered within SCORE_TOL standard errors
// of the mode and the effect on the acceptance rate is of order SCORE_TOL^2.
//
// What matters for correctness is not how tight this is but that it is *fixed*:
// together with a fixed starting value it makes the fit a deterministic function
// of the current state, so the birth and death moves build the same proposal and
// the chain stays reversible. Linero's code uses 0.1 but restarts from the
// node's previously stored mode, which makes the stopping point depend on the
// path taken to the current state.
const double SCORE_TOL = 1e-2;
const int MAX_SCORE_STEPS = 50;

// A trust region on the Fisher-scoring step, in standard errors of the current
// fit. Newton's method is only reliable near the mode: on a sharply asymmetric
// target it can be thrown far past it, and then crawl back so slowly that the
// step cap is reached with the fit nowhere near the mode -- which produces a
// proposal that is always rejected, so the leaf never moves and the same failure
// repeats every sweep. It is a permanent trap rather than slow mixing, and it
// showed up on a nuisance parameter carried as a forest pinned at depth zero,
// where one leaf holds a whole parameter instead of a small increment.
//
// Capping the step converges in a handful of iterations instead. Like SCORE_TOL,
// what matters for correctness is that the cap is *fixed*: the fit stays a
// deterministic function of the current state, so the forward and reverse moves
// still build the same proposal. And it applies only to the general path -- a
// quadratic target is reached exactly in one step from anywhere, and capping
// that would break an exactness the sampler relies on.
const double STEP_CAP = 4.0;

// A trust region on the *location* of the proposal, in its own standard errors.
// The Laplace fit is an independence proposal: it is centred on the mode and does
// not depend on where the leaf currently sits, which is what lets one fit serve
// both directions of the Metropolis ratio. That breaks down when the leaf is far
// from the mode -- the reverse density at the current value is then
// astronomically small, so the ratio is hugely negative however much the target
// improves, and the move can never be accepted. Measured on a nuisance parameter
// started twenty standard errors out: the target improved by 222 log points and
// the proposal cost 449, so the chain sat there forever.
//
// Damping the proposal toward the current value makes the step local, and the
// chain walks to the mode over a few sweeps. Reversibility is kept by building
// the reverse proposal the same way from the proposed value, which costs nothing
// because both come from the one shared fit; and when the cap does not bind both
// reduce to the fit itself, so the ratio is exactly the undamped one and every
// family that sits near its mode is unaffected.
const double PROPOSAL_CAP = 3.0;

inline double damped_mean(const Laplace1& fit, double from) {
  double delta = fit.mean - from;
  double cap = PROPOSAL_CAP * fit.sd;

  if (delta > cap) {
    return from + cap;
  }
  if (delta < -cap) {
    return from - cap;
  }

  return fit.mean;
}

// The Newton step, held inside the trust region.
inline double capped_step(double u, double j, bool quadratic) {
  double delta = u / j;

  if (quadratic) {
    return delta;
  }

  double cap = STEP_CAP / std::sqrt(j);

  if (delta > cap) {
    return cap;
  }
  if (delta < -cap) {
    return -cap;
  }

  return delta;
}

// Every Laplace fit starts from zero rather than from a neighboring node's
// value, for the same reason: the proposal must be a function of the current
// state alone, not of the path taken to it. The leaf prior shrinks towards
// zero, so this costs almost nothing in iterations.
const double SCORE_INIT = 0.0;

// Probability of attempting a birth or a death rather than a change, and the
// split between the two. A tree that is a single leaf has no branch to kill or
// rule to change, so a birth is forced; the acceptance ratios below use these
// same numbers so the boundary case stays reversible.
const double P_BIRTH_DEATH = 0.7;

double p_birth_move(const Node* root) {
  return root->is_leaf ? 1.0 : 0.5 * P_BIRTH_DEATH;
}

double p_death_move(const Node* root) {
  return root->is_leaf ? 0.0 : 0.5 * P_BIRTH_DEATH;
}

// The birth probability of the tree that results from collapsing this branch,
// which is a single leaf exactly when the branch is the root. Linero's code
// evaluates this on the pre-collapse tree and so uses 0.5 even when the root is
// collapsed, which unbalances the death move there.
double p_birth_after_death(const Node* branch) {
  return branch->is_root ? 1.0 : 0.5 * P_BIRTH_DEATH;
}

} // namespace

double Laplace2::log_dens(double x, double y) const {
  double dx = x - mean[0];
  double dy = y - mean[1];
  double quad = prec[0] * dx * dx + 2.0 * prec[1] * dx * dy + prec[2] * dy * dy;
  return -LN_2PI + 0.5 * log_det_prec - 0.5 * quad;
}

void Laplace2::draw(double* out) const {
  // If prec = L L' then mean + L^{-T} z has covariance prec^{-1}.
  double l11 = std::sqrt(prec[0]);
  double l21 = prec[1] / l11;
  double l22 = std::sqrt(std::max(prec[2] - l21 * l21, 1e-300));

  double z1 = norm_rand();
  double z2 = norm_rand();

  double d2 = z2 / l22;
  double d1 = (z1 - l21 * d2) / l11;

  out[0] = mean[0] + d1;
  out[1] = mean[1] + d2;
}

void Context::make_base(const Node* node, std::vector<double>& base) const {
  std::size_t n = node->idx.size();
  base.resize(n);
  const arma::mat& e = *eta;
  const double* wt = node->weights();

  if (wt == nullptr) {
    for (std::size_t k = 0; k < n; k++) {
      base[k] = e(h, node->idx[k]) - node->mu;
    }
    return;
  }

  for (std::size_t k = 0; k < n; k++) {
    base[k] = e(h, node->idx[k]) - wt[k] * node->mu;
  }
}

void Context::make_base_children(const Node* parent, std::vector<double>& base,
                                 std::vector<double>& w_left,
                                 std::vector<double>& w_right) const {
  std::size_t n = parent->idx.size();
  base.resize(n);
  w_left.resize(n);
  w_right.resize(n);
  const arma::mat& e = *eta;
  bool soft = parent->tree->hypers->soft;
  int gate = parent->tree->hypers->gate;
  const arma::mat& X = *parent->tree->X;
  double mu_left = parent->left->mu;
  double mu_right = parent->right->mu;

  const double* wt = parent->weights();

  for (std::size_t k = 0; k < n; k++) {
    int i = parent->idx[k];
    double g = left_prob(X(i, parent->var), parent->val,
                         parent->tree->bandwidth, soft, parent->na_rule, gate);
    double parent_wt = wt == nullptr ? 1.0 : wt[k];
    double wl = parent_wt * g;
    double wr = parent_wt - wl;
    w_left[k] = wl;
    w_right[k] = wr;
    base[k] = e(h, i) - wl * mu_left - wr * mu_right;
  }
}


// Every evaluation below exists in two forms. The per-observation loop is what
// a compiled family wants: the log density is a handful of operations, so
// materializing a vector of them costs more than it saves -- measured at 30% of
// the runtime for a Gaussian response. The blocked form builds the whole leaf's
// predictors and hands them over in one call, which is what a family whose log
// density is an R function needs, since there the fixed cost per call is
// everything. `Context::blocked` picks, and the default block implementations
// in Family sum in the same order as the loops here, so for any family that
// does not override them the two paths agree to the last bit. The test suite
// runs a fit both ways and requires exactly that.

void Context::reserve(std::size_t n, std::size_t copies) const {
  block.resize(copies * n * H);
  values.resize(copies * n);
  d1.resize(copies * n);
  d2.resize(copies * n);
}

// The predictors of the node's support, copied into a contiguous block with
// component h replaced by the leaf contribution the given parameter implies.
void Context::fill(const Node* node, const std::vector<double>& base, double mu,
                   std::size_t at) const {
  std::size_t n = node->idx.size();

  // With a single additive predictor there is nothing to copy: the one
  // component is overwritten.
  const double* wt = node->weights();

  if (H == 1) {
    double* dst = block.data() + at;

    if (wt == nullptr) {
      for (std::size_t k = 0; k < n; k++) {
        dst[k] = base[k] + mu;
      }

      return;
    }

    for (std::size_t k = 0; k < n; k++) {
      dst[k] = base[k] + wt[k] * mu;
    }

    return;
  }

  const arma::mat& e = *eta;

  if (wt == nullptr) {
    for (std::size_t k = 0; k < n; k++) {
      double* dst = block.data() + (at + k) * H;
      const double* src = e.colptr(node->idx[k]);
      std::copy(src, src + H, dst);
      dst[h] = base[k] + mu;
    }

    return;
  }

  for (std::size_t k = 0; k < n; k++) {
    double* dst = block.data() + (at + k) * H;
    const double* src = e.colptr(node->idx[k]);
    std::copy(src, src + H, dst);
    dst[h] = base[k] + wt[k] * mu;
  }
}

// The same, taking the leaf's current contribution from the predictor itself
// and shifting it, which is what the hot path wants: no base array is needed.
void Context::fill_at(const Node* node, double shift, std::size_t at) const {
  std::size_t n = node->idx.size();
  const arma::mat& e = *eta;

  const double* wt = node->weights();

  if (H == 1) {
    double* dst = block.data() + at;
    const double* src = e.memptr();

    if (wt == nullptr) {
      for (std::size_t k = 0; k < n; k++) {
        dst[k] = src[node->idx[k]] + shift;
      }

      return;
    }

    for (std::size_t k = 0; k < n; k++) {
      dst[k] = src[node->idx[k]] + wt[k] * shift;
    }

    return;
  }

  if (wt == nullptr) {
    for (std::size_t k = 0; k < n; k++) {
      double* dst = block.data() + (at + k) * H;
      const double* src = e.colptr(node->idx[k]);
      std::copy(src, src + H, dst);
      dst[h] += shift;
    }

    return;
  }

  for (std::size_t k = 0; k < n; k++) {
    double* dst = block.data() + (at + k) * H;
    const double* src = e.colptr(node->idx[k]);
    std::copy(src, src + H, dst);
    dst[h] += wt[k] * shift;
  }
}

// The two-child version, where the parent's support is divided between the
// children by the gate.
void Context::fill2(const Node* parent, const std::vector<double>& base,
                    const std::vector<double>& w_left,
                    const std::vector<double>& w_right, double mu_left,
                    double mu_right, std::size_t at) const {
  std::size_t n = parent->idx.size();

  if (H == 1) {
    double* dst = block.data() + at;

    for (std::size_t k = 0; k < n; k++) {
      dst[k] = base[k] + w_left[k] * mu_left + w_right[k] * mu_right;
    }

    return;
  }

  const arma::mat& e = *eta;

  for (std::size_t k = 0; k < n; k++) {
    double* dst = block.data() + (at + k) * H;
    const double* src = e.colptr(parent->idx[k]);
    std::copy(src, src + H, dst);
    dst[h] = base[k] + w_left[k] * mu_left + w_right[k] * mu_right;
  }
}

// The node's observations listed twice, so that a pair of leaf values can be
// evaluated in one block call rather than two.
void Context::make_pair_idx(const Node* node) const {
  std::size_t n = node->idx.size();
  pair_idx.resize(2 * n);
  std::copy(node->idx.begin(), node->idx.end(), pair_idx.begin());
  std::copy(node->idx.begin(), node->idx.end(), pair_idx.begin() + n);
}

double Context::log_f(const Node* node, const std::vector<double>& base,
                      double mu) const {
  double out = R::dnorm4(mu, 0.0, sigma_mu, 1);
  std::size_t n = node->idx.size();

  if (blocked) {
    reserve(n, 1);
    fill(node, base, mu, 0);
    family->logdens_block(node->idx.data(), static_cast<int>(n), block.data(),
                          values.data());

    for (std::size_t k = 0; k < n; k++) {
      out += values[k];
    }

    return out;
  }

  const arma::mat& e = *eta;
  const double* wt = node->weights();

  if (wt == nullptr) {
    for (std::size_t k = 0; k < n; k++) {
      int i = node->idx[k];
      std::copy(e.colptr(i), e.colptr(i) + H, scratch.begin());
      scratch[h] = base[k] + mu;
      out += family->logdens(i, scratch.data());
    }

    return out;
  }

  for (std::size_t k = 0; k < n; k++) {
    int i = node->idx[k];
    std::copy(e.colptr(i), e.colptr(i) + H, scratch.begin());
    scratch[h] = base[k] + wt[k] * mu;
    out += family->logdens(i, scratch.data());
  }

  return out;
}

void Context::score_info(const Node* node, const std::vector<double>& base,
                         double mu, double* score, double* info) const {
  double prec = 1.0 / (sigma_mu * sigma_mu);
  double s_acc = -mu * prec;
  double j_acc = prec;
  std::size_t n = node->idx.size();

  if (blocked) {
    reserve(n, 1);
    fill(node, base, mu, 0);
    family->score_info_block(node->idx.data(), static_cast<int>(n),
                             block.data(), h, d1.data(), d2.data());

    const double* wt = node->weights();

    for (std::size_t k = 0; k < n; k++) {
      double wk = wt == nullptr ? 1.0 : wt[k];
      s_acc += wk * d1[k];
      j_acc += wk * wk * d2[k];
    }

    *score = s_acc;
    *info = j_acc;
    return;
  }

  const arma::mat& e = *eta;
  const double* wt = node->weights();

  for (std::size_t k = 0; k < n; k++) {
    int i = node->idx[k];
    double wk = wt == nullptr ? 1.0 : wt[k];
    std::copy(e.colptr(i), e.colptr(i) + H, scratch.begin());
    scratch[h] = base[k] + wk * mu;
    double a;
    double b;
    family->score_info(i, scratch.data(), h, &a, &b);
    s_acc += wk * a;
    j_acc += wk * wk * b;
  }

  *score = s_acc;
  *info = j_acc;
}

void Context::log_f_pair(const Node* node, const std::vector<double>& base,
                         double mu_a, double mu_b, double* out_a,
                         double* out_b) const {
  double a_acc = R::dnorm4(mu_a, 0.0, sigma_mu, 1);
  double b_acc = R::dnorm4(mu_b, 0.0, sigma_mu, 1);
  std::size_t n = node->idx.size();

  if (blocked) {
    reserve(n, 2);
    fill(node, base, mu_a, 0);
    fill(node, base, mu_b, n);
    make_pair_idx(node);
    family->logdens_block(pair_idx.data(), static_cast<int>(2 * n),
                          block.data(), values.data());

    for (std::size_t k = 0; k < n; k++) {
      a_acc += values[k];
    }

    for (std::size_t k = 0; k < n; k++) {
      b_acc += values[n + k];
    }

    *out_a = a_acc;
    *out_b = b_acc;
    return;
  }

  const arma::mat& e = *eta;

  const double* wt = node->weights();

  for (std::size_t k = 0; k < n; k++) {
    int i = node->idx[k];
    std::copy(e.colptr(i), e.colptr(i) + H, scratch.begin());
    double centered = base[k];
    double weight = wt == nullptr ? 1.0 : wt[k];
    scratch[h] = centered + weight * mu_a;
    a_acc += family->logdens(i, scratch.data());
    scratch[h] = centered + weight * mu_b;
    b_acc += family->logdens(i, scratch.data());
  }

  *out_a = a_acc;
  *out_b = b_acc;
}

void Context::score_info_at(const Node* node, double mu_ref, double mu,
                            double* score, double* info) const {
  double prec = 1.0 / (sigma_mu * sigma_mu);
  double s_acc = -mu * prec;
  double j_acc = prec;
  double shift = mu - mu_ref;
  std::size_t n = node->idx.size();

  if (blocked) {
    reserve(n, 1);
    fill_at(node, shift, 0);
    family->score_info_block(node->idx.data(), static_cast<int>(n),
                             block.data(), h, d1.data(), d2.data());

    const double* wt = node->weights();

    for (std::size_t k = 0; k < n; k++) {
      double weight = wt == nullptr ? 1.0 : wt[k];
      s_acc += weight * d1[k];
      j_acc += weight * weight * d2[k];
    }

    *score = s_acc;
    *info = j_acc;
    return;
  }

  const arma::mat& e = *eta;

  const double* wt = node->weights();

  for (std::size_t k = 0; k < n; k++) {
    int i = node->idx[k];
    std::copy(e.colptr(i), e.colptr(i) + H, scratch.begin());
    double weight = wt == nullptr ? 1.0 : wt[k];
    scratch[h] += weight * shift;
    double a;
    double b;
    family->score_info(i, scratch.data(), h, &a, &b);
    s_acc += weight * a;
    j_acc += weight * weight * b;
  }

  *score = s_acc;
  *info = j_acc;
}

void Context::log_f_pair_at(const Node* node, double mu_ref, double mu_a,
                            double mu_b, double* out_a, double* out_b) const {
  double a_acc = R::dnorm4(mu_a, 0.0, sigma_mu, 1);
  double b_acc = R::dnorm4(mu_b, 0.0, sigma_mu, 1);
  double shift_a = mu_a - mu_ref;
  double shift_b = mu_b - mu_ref;
  std::size_t n = node->idx.size();

  if (blocked) {
    reserve(n, 2);
    fill_at(node, shift_a, 0);
    fill_at(node, shift_b, n);
    make_pair_idx(node);
    family->logdens_block(pair_idx.data(), static_cast<int>(2 * n),
                          block.data(), values.data());

    for (std::size_t k = 0; k < n; k++) {
      a_acc += values[k];
    }

    for (std::size_t k = 0; k < n; k++) {
      b_acc += values[n + k];
    }

    *out_a = a_acc;
    *out_b = b_acc;
    return;
  }

  const arma::mat& e = *eta;

  const double* wt = node->weights();

  for (std::size_t k = 0; k < n; k++) {
    int i = node->idx[k];
    std::copy(e.colptr(i), e.colptr(i) + H, scratch.begin());
    double weight = wt == nullptr ? 1.0 : wt[k];
    double centered = scratch[h];
    scratch[h] = centered + weight * shift_a;
    a_acc += family->logdens(i, scratch.data());
    scratch[h] = centered + weight * shift_b;
    b_acc += family->logdens(i, scratch.data());
  }

  *out_a = a_acc;
  *out_b = b_acc;
}

double Context::log_f2(const Node* parent, const std::vector<double>& base,
                       const std::vector<double>& w_left,
                       const std::vector<double>& w_right,
                       double mu_left, double mu_right) const {
  double out = R::dnorm4(mu_left, 0.0, sigma_mu, 1) +
    R::dnorm4(mu_right, 0.0, sigma_mu, 1);
  std::size_t n = parent->idx.size();

  if (blocked) {
    reserve(n, 1);
    fill2(parent, base, w_left, w_right, mu_left, mu_right, 0);
    family->logdens_block(parent->idx.data(), static_cast<int>(n), block.data(),
                          values.data());

    for (std::size_t k = 0; k < n; k++) {
      out += values[k];
    }

    return out;
  }

  const arma::mat& e = *eta;

  for (std::size_t k = 0; k < n; k++) {
    int i = parent->idx[k];
    std::copy(e.colptr(i), e.colptr(i) + H, scratch.begin());
    scratch[h] = base[k] + w_left[k] * mu_left + w_right[k] * mu_right;
    out += family->logdens(i, scratch.data());
  }

  return out;
}

void Context::score_info2(const Node* parent, const std::vector<double>& base,
                          const std::vector<double>& w_left,
                          const std::vector<double>& w_right,
                          double mu_left, double mu_right,
                          double* g, double* info) const {
  double prec = 1.0 / (sigma_mu * sigma_mu);
  g[0] = -mu_left * prec;
  g[1] = -mu_right * prec;
  info[0] = prec;
  info[1] = 0.0;
  info[2] = prec;
  std::size_t n = parent->idx.size();

  if (blocked) {
    reserve(n, 1);
    fill2(parent, base, w_left, w_right, mu_left, mu_right, 0);
    family->score_info_block(parent->idx.data(), static_cast<int>(n),
                             block.data(), h, d1.data(), d2.data());

    for (std::size_t k = 0; k < n; k++) {
      double wl = w_left[k];
      double wr = w_right[k];
      g[0] += wl * d1[k];
      g[1] += wr * d1[k];
      info[0] += wl * wl * d2[k];
      info[1] += wl * wr * d2[k];
      info[2] += wr * wr * d2[k];
    }

    return;
  }

  const arma::mat& e = *eta;

  for (std::size_t k = 0; k < n; k++) {
    int i = parent->idx[k];
    std::copy(e.colptr(i), e.colptr(i) + H, scratch.begin());
    scratch[h] = base[k] + w_left[k] * mu_left + w_right[k] * mu_right;
    double a;
    double b;
    family->score_info(i, scratch.data(), h, &a, &b);
    double wl = w_left[k];
    double wr = w_right[k];
    g[0] += wl * a;
    g[1] += wr * a;
    info[0] += wl * wl * b;
    info[1] += wl * wr * b;
    info[2] += wr * wr * b;
  }
}

// Fisher scoring from a fixed start. The start is fixed, and the tolerance
// with it, so that the fit is a deterministic function of the state and the
// reverse move reconstructs the same proposal -- which is what reversibility
// needs. A quadratic target is reached in one step from anywhere, so the loop is
// cut short for one, exactly rather than approximately.
// The log target, its score and its information at one point, in a single pass.
// The blocked path keeps its two separate calls: it exists for families that
// reach back into R, and none of those is quadratic, so this is never their hot
// path.
// One pass giving the log target, its score and its information. The family
// accumulates the three sums itself, which lets a compiled family inline its own
// arithmetic into the loop rather than being called through the vtable four
// times per observation. The prior is added here, since it is the same Gaussian
// whatever the family.
void Context::log_f_score_info(const Node* node, const std::vector<double>& base,
                               double mu, double* f, double* s,
                               double* j) const {
  int n = static_cast<int>(node->idx.size());

  double fa;
  double sa;
  double ja;

  if (blocked) {
    reserve(n, 1);
    fill(node, base, mu, 0);
    accumulate1(node, n, &fa, &sa, &ja);
  }
  else {
    accumulate1_from(node, n, base, mu, &fa, &sa, &ja);
  }

  double prec = 1.0 / (sigma_mu * sigma_mu);
  *f = fa + R::dnorm4(mu, 0.0, sigma_mu, 1);
  *s = sa - mu * prec;
  *j = ja + prec;
}

// Dispatch to the family's own accumulator, or to the dynamically dispatched
// one when the caller has asked for it.
void Context::accumulate1(const Node* node, int n, double* f, double* s,
                          double* j) const {
  if (generic_accumulate) {
    family->accumulate1_generic(node->idx.data(), node->weights(), n,
                                block.data(), h, f, s, j);
    return;
  }

  family->accumulate1(node->idx.data(), node->weights(), n, block.data(), h, f,
                      s, j);
}

void Context::accumulate1_from(const Node* node, int n,
                               const std::vector<double>& base, double mu,
                               double* f, double* s, double* j) const {
  if (generic_accumulate) {
    family->accumulate1_from_generic(node->idx.data(), node->weights(), n, *eta,
                                     base.data(), mu, h, f, s, j);
    return;
  }

  family->accumulate1_from(node->idx.data(), node->weights(), n, *eta,
                           base.data(), mu, h, f, s, j);
}

void Context::accumulate2_from(const Node* parent, int n,
                               const std::vector<double>& base,
                               const std::vector<double>& w_left,
                               const std::vector<double>& w_right,
                               double mu_left, double mu_right, double* f,
                               double* g, double* info) const {
  if (generic_accumulate) {
    family->accumulate2_from_generic(parent->idx.data(), w_left.data(),
                                     w_right.data(), n, *eta, base.data(),
                                     mu_left, mu_right, h, f, g, info);
    return;
  }

  family->accumulate2_from(parent->idx.data(), w_left.data(), w_right.data(), n,
                           *eta, base.data(), mu_left, mu_right, h, f, g, info);
}

void Context::accumulate1_at(const Node* node, int n, double shift, double* f,
                             double* s, double* j) const {
  if (generic_accumulate) {
    family->accumulate1_at_generic(node->idx.data(), node->weights(), n, *eta, h,
                                   shift, f, s, j);
    return;
  }

  family->accumulate1_at(node->idx.data(), node->weights(), n, *eta, h, shift, f,
                         s, j);
}

void Context::accumulate2(const Node* parent, int n,
                          const std::vector<double>& w_left,
                          const std::vector<double>& w_right, double* f,
                          double* g, double* info) const {
  if (generic_accumulate) {
    family->accumulate2_generic(parent->idx.data(), w_left.data(),
                                w_right.data(), n, block.data(), h, f, g, info);
    return;
  }

  family->accumulate2(parent->idx.data(), w_left.data(), w_right.data(), n,
                      block.data(), h, f, g, info);
}

void Context::log_f_score_info_at(const Node* node, double mu_ref, double mu,
                                  double* f, double* s, double* j) const {
  int n = static_cast<int>(node->idx.size());

  double fa;
  double sa;
  double ja;

  if (blocked) {
    // A family whose log density is an R function has to be handed a whole leaf
    // at a time, so the block is what it wants and there is nothing to fuse.
    reserve(n, 1);
    fill_at(node, mu - mu_ref, 0);
    accumulate1(node, n, &fa, &sa, &ja);
  }
  else {
    accumulate1_at(node, n, mu - mu_ref, &fa, &sa, &ja);
  }

  double prec = 1.0 / (sigma_mu * sigma_mu);
  *f = fa + R::dnorm4(mu, 0.0, sigma_mu, 1);
  *s = sa - mu * prec;
  *j = ja + prec;
}

void Context::log_f_score_info2(const Node* parent,
                                const std::vector<double>& base,
                                const std::vector<double>& w_left,
                                const std::vector<double>& w_right,
                                double mu_left, double mu_right, double* f,
                                double* g, double* info) const {
  int n = static_cast<int>(parent->idx.size());

  double prec = 1.0 / (sigma_mu * sigma_mu);
  g[0] = -mu_left * prec;
  g[1] = -mu_right * prec;
  info[0] = prec;
  info[1] = 0.0;
  info[2] = prec;

  double fa;

  if (blocked) {
    reserve(n, 1);
    fill2(parent, base, w_left, w_right, mu_left, mu_right, 0);
    accumulate2(parent, n, w_left, w_right, &fa, g, info);
  }
  else {
    accumulate2_from(parent, n, base, w_left, w_right, mu_left, mu_right, &fa, g,
                     info);
  }

  *f = fa + R::dnorm4(mu_left, 0.0, sigma_mu, 1) +
    R::dnorm4(mu_right, 0.0, sigma_mu, 1);
}

Laplace1 fit_laplace1(Context& ctx, const Node* node,
                      const std::vector<double>& base) {
  double m = SCORE_INIT;
  double u;
  double j;
  ctx.score_info(node, base, m, &u, &j);
  int max_steps = ctx.quadratic ? 1 : MAX_SCORE_STEPS;

  for (int step = 0; step < max_steps; step++) {
    if (!(j > 0.0) || !std::isfinite(u)) {
      break;
    }
    // A quadratic target is reached exactly by one step from anywhere, so there
    // is no tolerance to meet and stopping short of the step would only return
    // a worse proposal.
    if (!ctx.quadratic && std::fabs(u) <= SCORE_TOL * std::sqrt(j)) {
      break;
    }
    double proposed = m + capped_step(u, j, ctx.quadratic);
    if (!std::isfinite(proposed)) {
      break;
    }
    m = proposed;
    ctx.score_info(node, base, m, &u, &j);
  }

  Laplace1 out;
  out.mean = m;
  // The prior contributes 1 / sigma_mu^2 to the information, so the curvature
  // is bounded below and the proposal is never wider than the prior.
  out.sd = (j > 0.0 && std::isfinite(j)) ? std::pow(j, -0.5) : ctx.sigma_mu;
  return out;
}

Laplace1 fit_laplace1_at(Context& ctx, const Node* node, double mu_ref) {
  double m = SCORE_INIT;
  double u;
  double j;
  ctx.score_info_at(node, mu_ref, m, &u, &j);
  int max_steps = ctx.quadratic ? 1 : MAX_SCORE_STEPS;

  for (int step = 0; step < max_steps; step++) {
    if (!(j > 0.0) || !std::isfinite(u)) {
      break;
    }
    if (!ctx.quadratic && std::fabs(u) <= SCORE_TOL * std::sqrt(j)) {
      break;
    }
    double proposed = m + capped_step(u, j, ctx.quadratic);
    if (!std::isfinite(proposed)) {
      break;
    }
    m = proposed;
    ctx.score_info_at(node, mu_ref, m, &u, &j);
  }

  Laplace1 out;
  out.mean = m;
  out.sd = (j > 0.0 && std::isfinite(j)) ? std::pow(j, -0.5) : ctx.sigma_mu;
  return out;
}

Laplace2 fit_laplace2(Context& ctx, const Node* parent,
                      const std::vector<double>& base,
                      const std::vector<double>& w_left,
                      const std::vector<double>& w_right) {
  double m[2] = {SCORE_INIT, SCORE_INIT};
  double g[2];
  double info[3];
  int max_steps = ctx.quadratic ? 1 : MAX_SCORE_STEPS;

  for (int step = 0; step < max_steps; step++) {
    ctx.score_info2(parent, base, w_left, w_right, m[0], m[1], g, info);
    double det = info[0] * info[2] - info[1] * info[1];
    if (!(det > 0.0) || !std::isfinite(det)) {
      break;
    }
    double d0 = (info[2] * g[0] - info[1] * g[1]) / det;
    double d1 = (info[0] * g[1] - info[1] * g[0]) / det;
    // The squared score in standard-error units, which is the multivariate
    // form of the relative stopping rule used for the univariate fit.
    double scaled = g[0] * d0 + g[1] * d1;
    if (!std::isfinite(scaled)) {
      break;
    }
    if (!ctx.quadratic && scaled <= SCORE_TOL * SCORE_TOL) {
      break;
    }
    if (!std::isfinite(m[0] + d0) || !std::isfinite(m[1] + d1)) {
      break;
    }
    m[0] += d0;
    m[1] += d1;
  }

  ctx.score_info2(parent, base, w_left, w_right, m[0], m[1], g, info);

  Laplace2 out;
  out.mean[0] = m[0];
  out.mean[1] = m[1];
  out.prec[0] = info[0];
  out.prec[1] = info[1];
  out.prec[2] = info[2];

  double det = info[0] * info[2] - info[1] * info[1];
  if (!(det > 0.0) || !std::isfinite(det)) {
    double prec = 1.0 / (ctx.sigma_mu * ctx.sigma_mu);
    out.prec[0] = prec;
    out.prec[1] = 0.0;
    out.prec[2] = prec;
    det = prec * prec;
  }
  out.log_det_prec = std::log(det);
  return out;
}

// --- the quadratic targets --------------------------------------------------

namespace {

// The univariate fit, given the score and information at `ref`. Degenerate
// curvature falls back to the prior width, as the Fisher-scoring loop does.
Laplace1 laplace_from(double ref, double d1, double d2, double sigma_mu) {
  Laplace1 out;

  if (d2 > 0.0 && std::isfinite(d2) && std::isfinite(d1)) {
    out.mean = ref + d1 / d2;
    out.sd = std::pow(d2, -0.5);
  }
  else {
    out.mean = ref;
    out.sd = sigma_mu;
  }

  return out;
}

Laplace2 laplace2_from(const double* g, const double* info, double sigma_mu) {
  Laplace2 out;
  double det = info[0] * info[2] - info[1] * info[1];

  if (det > 0.0 && std::isfinite(det) && std::isfinite(g[0]) &&
      std::isfinite(g[1])) {
    out.mean[0] = (info[2] * g[0] - info[1] * g[1]) / det;
    out.mean[1] = (info[0] * g[1] - info[1] * g[0]) / det;
    out.prec[0] = info[0];
    out.prec[1] = info[1];
    out.prec[2] = info[2];
  }
  else {
    double prec = 1.0 / (sigma_mu * sigma_mu);
    out.mean[0] = 0.0;
    out.mean[1] = 0.0;
    out.prec[0] = prec;
    out.prec[1] = 0.0;
    out.prec[2] = prec;
    det = prec * prec;
  }

  out.log_det_prec = std::log(det);
  return out;
}

} // namespace

namespace {

// The tolerance for the scalar Newton iteration the exponential form uses. It
// runs on three numbers rather than on the data, so there is no reason to stop
// early: iterating to the machine's limit keeps the fit a function of the
// target's coefficients alone, which is what lets the forward and reverse moves
// reconstruct the same proposal.
const double EXP_TOL = 1e-12;
const int EXP_STEPS = 200;

// Whether the exponential form is usable here. It needs every observation in the
// node to carry weight one, which is what a hard rule gives and a soft one does
// not: a soft rule hands each observation its own exponent exp(s * w * mu), and a
// sum of those is not determined by three numbers.
bool exponential_usable(const Node* node, TargetForm form) {
  if (form != TARGET_EXP_UP && form != TARGET_EXP_DOWN) {
    return false;
  }

  return !node->tree->hypers->soft;
}

// Solve score(mu) = 0 by Newton's method, from a fixed start so that the answer
// depends on the target and not on where the expansion was taken.
double exponential_mode(double a, double b, double rate, double prec,
                        double* info_out) {
  double m = SCORE_INIT;
  double rate2 = rate * rate;

  for (int step = 0; step < EXP_STEPS; step++) {
    double ex = b * std::exp(rate * m);
    double score = a + rate * ex - m * prec;
    double info = prec - rate2 * ex;

    if (!(info > 0.0) || !std::isfinite(score) || !std::isfinite(info)) {
      break;
    }

    if (std::fabs(score) <= EXP_TOL * std::sqrt(info)) {
      *info_out = info;
      return m;
    }

    double proposed = m + score / info;

    if (!std::isfinite(proposed)) {
      break;
    }

    m = proposed;
  }

  double ex = b * std::exp(rate * m);
  *info_out = prec - rate2 * ex;
  return m;
}

} // namespace

Target1::Target1(Context& ctx, const Node* node, const std::vector<double>* base,
                 double ref)
  : ctx_(&ctx), node_(node), base_(base), ref_(ref), mode_(PASSES), f_(0.0),
    d1_(0.0), d2_(0.0), rate_(0.0), a_(0.0), b_(0.0), c_(0.0) {
  if (!ctx.exact_quadratic) {
    return;
  }

  bool quadratic = ctx.quadratic;
  bool exponential = exponential_usable(node, ctx.form);

  if (!quadratic && !exponential) {
    return;
  }

  if (base_ != nullptr) {
    ctx_->log_f_score_info(node_, *base_, ref_, &f_, &d1_, &d2_);
  }
  else {
    ctx_->log_f_score_info_at(node_, ref_, ref_, &f_, &d1_, &d2_);
  }

  if (quadratic) {
    mode_ = CLOSED;
    return;
  }

  // Strip the leaf prior, which is Gaussian and so known exactly, to leave the
  // likelihood's own value, slope and curvature at ref_. Then read off the
  // three coefficients of c + a * mu + b * exp(r * mu), whose curvature is
  // r^2 * b * exp(r * mu).
  double prec = 1.0 / (ctx.sigma_mu * ctx.sigma_mu);
  double value = f_ - R::dnorm4(ref_, 0.0, ctx.sigma_mu, 1);
  double slope = d1_ + ref_ * prec;
  double curve = prec - d2_;

  // b exp(r ref) = curve / r^2 and a = slope - curve / r, which at the rates of
  // +1 and -1 the Poisson and gamma have reduces to the sign arithmetic this
  // used to do, since there 1 / r == r and r^2 == 1.
  rate_ = ctx_->family->exp_rate(ctx_->h);
  double rate2 = rate_ * rate_;
  a_ = slope - curve / rate_;
  b_ = curve * std::exp(-rate_ * ref_) / rate2;
  c_ = value - a_ * ref_ - curve / rate2;

  // A non-finite coefficient means the predictor has run somewhere the
  // exponential cannot be represented; fall back rather than propagate it.
  if (!std::isfinite(a_) || !std::isfinite(b_) || !std::isfinite(c_)) {
    return;
  }

  mode_ = EXPONENTIAL;
}

double Target1::exp_score(double mu) const {
  return a_ + rate_ * b_ * std::exp(rate_ * mu) -
    mu / (ctx_->sigma_mu * ctx_->sigma_mu);
}

double Target1::exp_info(double mu) const {
  return 1.0 / (ctx_->sigma_mu * ctx_->sigma_mu) -
    rate_ * rate_ * b_ * std::exp(rate_ * mu);
}

double Target1::log_f(double mu) const {
  if (mode_ == CLOSED) {
    double d = mu - ref_;
    return f_ + d1_ * d - 0.5 * d2_ * d * d;
  }

  if (mode_ == EXPONENTIAL) {
    return c_ + a_ * mu + b_ * std::exp(rate_ * mu) +
      R::dnorm4(mu, 0.0, ctx_->sigma_mu, 1);
  }

  if (base_ != nullptr) {
    return ctx_->log_f(node_, *base_, mu);
  }

  double a;
  double b;
  ctx_->log_f_pair_at(node_, ref_, mu, ref_, &a, &b);
  return a;
}

Laplace1 Target1::laplace() const {
  if (mode_ == CLOSED) {
    return laplace_from(ref_, d1_, d2_, ctx_->sigma_mu);
  }

  if (mode_ == EXPONENTIAL) {
    double prec = 1.0 / (ctx_->sigma_mu * ctx_->sigma_mu);
    double info;
    double m = exponential_mode(a_, b_, rate_, prec, &info);

    Laplace1 out;

    if (info > 0.0 && std::isfinite(info) && std::isfinite(m)) {
      out.mean = m;
      out.sd = std::pow(info, -0.5);
    }
    else {
      out.mean = 0.0;
      out.sd = ctx_->sigma_mu;
    }

    return out;
  }

  if (base_ != nullptr) {
    return fit_laplace1(*ctx_, node_, *base_);
  }

  return fit_laplace1_at(*ctx_, node_, ref_);
}

Target2::Target2(Context& ctx, const Node* parent,
                 const std::vector<double>& base,
                 const std::vector<double>& w_left,
                 const std::vector<double>& w_right)
  : ctx_(&ctx), parent_(parent), base_(&base), w_left_(&w_left),
    w_right_(&w_right), mode_(PASSES), f_(0.0), rate_(0.0), c_(0.0) {
  g_[0] = 0.0;
  g_[1] = 0.0;
  info_[0] = 0.0;
  info_[1] = 0.0;
  info_[2] = 0.0;
  a_[0] = 0.0;
  a_[1] = 0.0;
  b_[0] = 0.0;
  b_[1] = 0.0;

  if (!ctx.exact_quadratic) {
    return;
  }

  bool quadratic = ctx.quadratic;
  bool exponential = exponential_usable(parent, ctx.form);

  if (!quadratic && !exponential) {
    return;
  }

  // Expanded about zero, which is where a fresh pair of children starts.
  ctx_->log_f_score_info2(parent_, *base_, *w_left_, *w_right_, 0.0, 0.0, &f_,
                          g_, info_);

  if (quadratic) {
    mode_ = CLOSED;
    return;
  }

  // A hard rule sends each observation to exactly one child, so no observation
  // contributes to the cross term and the two children separate. If that is not
  // what came back, the assumption behind the form does not hold here.
  if (info_[1] != 0.0) {
    return;
  }

  double prec = 1.0 / (ctx.sigma_mu * ctx.sigma_mu);
  rate_ = ctx_->family->exp_rate(ctx_->h);
  double rate2 = rate_ * rate_;
  double curve_left = prec - info_[0];
  double curve_right = prec - info_[2];
  b_[0] = curve_left / rate2;
  b_[1] = curve_right / rate2;
  a_[0] = g_[0] - curve_left / rate_;
  a_[1] = g_[1] - curve_right / rate_;

  // Only the two constants' sum is ever needed, which is as well: one pass
  // returns the two children's log densities already added together.
  c_ = f_ - b_[0] - b_[1] - 2.0 * R::dnorm4(0.0, 0.0, ctx.sigma_mu, 1);

  if (!std::isfinite(a_[0]) || !std::isfinite(a_[1]) ||
      !std::isfinite(b_[0]) || !std::isfinite(b_[1]) || !std::isfinite(c_)) {
    return;
  }

  mode_ = EXPONENTIAL;
}

double Target2::log_f(double mu_left, double mu_right) const {
  if (mode_ == CLOSED) {
    return f_ + g_[0] * mu_left + g_[1] * mu_right -
      0.5 * (info_[0] * mu_left * mu_left +
             2.0 * info_[1] * mu_left * mu_right +
             info_[2] * mu_right * mu_right);
  }

  if (mode_ == EXPONENTIAL) {
    return c_ + a_[0] * mu_left + a_[1] * mu_right +
      b_[0] * std::exp(rate_ * mu_left) + b_[1] * std::exp(rate_ * mu_right) +
      R::dnorm4(mu_left, 0.0, ctx_->sigma_mu, 1) +
      R::dnorm4(mu_right, 0.0, ctx_->sigma_mu, 1);
  }

  return ctx_->log_f2(parent_, *base_, *w_left_, *w_right_, mu_left, mu_right);
}

Laplace2 Target2::laplace() const {
  if (mode_ == CLOSED) {
    return laplace2_from(g_, info_, ctx_->sigma_mu);
  }

  if (mode_ == EXPONENTIAL) {
    double prec = 1.0 / (ctx_->sigma_mu * ctx_->sigma_mu);
    double info[3];
    double g[2];
    info[1] = 0.0;

    // Two one-dimensional problems, since the cross curvature is zero.
    for (int side = 0; side < 2; side++) {
      double curvature;
      double m = exponential_mode(a_[side], b_[side], rate_, prec, &curvature);
      double slot = side == 0 ? 0.0 : 2.0;
      info[static_cast<int>(slot)] = curvature;
      // Expressed as a score at the origin, so that laplace2_from() recovers
      // the same mode it would have found by solving the system.
      g[side] = curvature * m;
    }

    return laplace2_from(g, info, ctx_->sigma_mu);
  }

  return fit_laplace2(*ctx_, parent_, *base_, *w_left_, *w_right_);
}

namespace {

void apply_leaf_delta(Node* leaf, arma::mat& eta, int h, double new_mu) {
  double delta = new_mu - leaf->mu;
  if (delta == 0.0) {
    return;
  }
  std::size_t n = leaf->idx.size();
  const double* wt = leaf->weights();

  if (wt == nullptr) {
    for (std::size_t k = 0; k < n; k++) {
      eta(h, leaf->idx[k]) += delta;
    }
  }
  else {
    for (std::size_t k = 0; k < n; k++) {
      eta(h, leaf->idx[k]) += wt[k] * delta;
    }
  }

  leaf->mu = new_mu;
}

// Rewrite the predictor for the observations a node supports, which is how a
// completed birth, death or change move is committed.
void commit_children(const Node* parent, arma::mat& eta, int h,
                     const std::vector<double>& base,
                     const std::vector<double>& w_left,
                     const std::vector<double>& w_right,
                     double mu_left, double mu_right) {
  std::size_t n = parent->idx.size();
  for (std::size_t k = 0; k < n; k++) {
    eta(h, parent->idx[k]) = base[k] + w_left[k] * mu_left +
      w_right[k] * mu_right;
  }
}

void commit_single(const Node* node, arma::mat& eta, int h,
                   const std::vector<double>& base, double mu) {
  std::size_t n = node->idx.size();
  const double* wt = node->weights();

  if (wt == nullptr) {
    for (std::size_t k = 0; k < n; k++) {
      eta(h, node->idx[k]) = base[k] + mu;
    }

    return;
  }

  for (std::size_t k = 0; k < n; k++) {
    eta(h, node->idx[k]) = base[k] + wt[k] * mu;
  }
}

void node_birth(Tree* tree, Context& ctx, const Hypers& hypers) {

  std::vector<Node*>& leaf_list = ctx.buf_leaves;
  leaf_list.clear();
  leaves(tree->root, leaf_list);
  Node* leaf = rand_node(leaf_list);

  std::vector<double>& base = ctx.buf_base;
  ctx.make_base(leaf, base);

  double rho_d = grow_prob(&hypers, leaf->depth);

  // A forest pinned at depth zero carries a branching probability of zero, so
  // the move is impossible and the ratio below would be -Inf. Returning here
  // says so plainly and skips a Laplace fit that could only be rejected.
  if (!(rho_d > 0.0)) {
    return;
  }
  double rho_d1 = grow_prob(&hypers, leaf->depth + 1);

  // One object for the parent's target: the log density before the move and the
  // fit the reverse death move would propose from both come out of it, in one
  // pass when the family is quadratic.
  Target1 parent_target(ctx, leaf, &base, leaf->mu);
  double log_f_before = parent_target.log_f(leaf->mu);
  Laplace1 merge_fit = parent_target.laplace();
  double log_g_reverse = merge_fit.log_dens(leaf->mu);

  double p_forward = std::log(p_birth_move(tree->root)) -
    std::log(static_cast<double>(leaf_list.size()));

  // The children start at zero, so the predictor still carries the parent's own
  // contribution. `base` already has that removed and the split does not change
  // which observations reach this node, so it is the correct base for the pair.
  // The weights come out of the split itself: dividing the support and recording
  // what each side got are the same pass, and were two.
  std::vector<double>& w_left = ctx.buf_left;
  std::vector<double>& w_right = ctx.buf_right;
  leaf->birth_leaves(&w_left, &w_right);

  Target2 child_target(ctx, leaf, base, w_left, w_right);
  Laplace2 fit = child_target.laplace();
  double proposal[2];
  fit.draw(proposal);

  double log_f_after = child_target.log_f(proposal[0], proposal[1]);
  double log_g_forward = fit.log_dens(proposal[0], proposal[1]);

  double p_backward = std::log(p_death_move(tree->root)) -
    std::log(static_cast<double>(num_not_grand_branches(tree->root)));

  double log_ratio = log_f_after - log_f_before +
    std::log(rho_d) + 2.0 * std::log1p(-rho_d1) - std::log1p(-rho_d) +
    p_backward - p_forward + log_g_reverse - log_g_forward;

  if (std::isfinite(log_ratio) && std::log(unif_rand()) < log_ratio) {
    leaf->left->mu = proposal[0];
    leaf->right->mu = proposal[1];
    commit_children(leaf, *ctx.eta, ctx.h, base, w_left, w_right, proposal[0],
                    proposal[1]);
  }
  else {
    leaf->delete_leaves();
    leaf->var = 0;
    leaf->group = 0;
    leaf->na_rule = NA_LEFT;
  }
}

void node_death(Tree* tree, Context& ctx, const Hypers& hypers) {

  std::vector<Node*>& ngb = ctx.buf_branches;
  ngb.clear();
  not_grand_branches(tree->root, ngb);
  if (ngb.empty()) {
    return;
  }
  Node* branch = rand_node(ngb);
  Node* left = branch->left;
  Node* right = branch->right;

  double rho_d = grow_prob(&hypers, branch->depth);
  double rho_d1 = grow_prob(&hypers, branch->depth + 1);

  // base is the predictor with this whole subtree removed, which serves as the
  // base both for the pair of children and for the merged leaf: a merged leaf
  // has exactly the branch's own support and weights.
  std::vector<double>& base = ctx.buf_base;
  std::vector<double>& w_left = ctx.buf_left;
  std::vector<double>& w_right = ctx.buf_right;
  ctx.make_base_children(branch, base, w_left, w_right);

  Target2 child_target(ctx, branch, base, w_left, w_right);
  double log_f_before = child_target.log_f(left->mu, right->mu);
  // The reverse birth move would propose the pair from this same fit.
  Laplace2 split_fit = child_target.laplace();
  double log_g_reverse = split_fit.log_dens(left->mu, right->mu);

  int num_leaf_before = num_leaves(tree->root);
  double p_forward = std::log(p_death_move(tree->root)) -
    std::log(static_cast<double>(ngb.size()));
  double p_backward = std::log(p_birth_after_death(branch)) -
    std::log(static_cast<double>(num_leaf_before - 1));

  // Detach rather than delete, so that a rejection can put them back.
  branch->left = nullptr;
  branch->right = nullptr;
  branch->is_leaf = true;

  Target1 merged_target(ctx, branch, &base, 0.0);
  Laplace1 merge_fit = merged_target.laplace();
  double mu_new = merge_fit.draw();
  double log_f_after = merged_target.log_f(mu_new);
  double log_g_forward = merge_fit.log_dens(mu_new);

  double log_ratio = log_f_after - log_f_before +
    std::log1p(-rho_d) - std::log(rho_d) - 2.0 * std::log1p(-rho_d1) +
    p_backward - p_forward + log_g_reverse - log_g_forward;

  if (std::isfinite(log_ratio) && std::log(unif_rand()) < log_ratio) {
    tree->give_node(left);
    tree->give_node(right);
    branch->mu = mu_new;
    branch->var = 0;
    branch->group = 0;
    branch->na_rule = NA_LEFT;
    commit_single(branch, *ctx.eta, ctx.h, base, mu_new);
  }
  else {
    branch->left = left;
    branch->right = right;
    branch->is_leaf = false;
  }
}

void change_rule(Tree* tree, Context& ctx) {

  std::vector<Node*>& ngb = ctx.buf_branches;
  ngb.clear();
  not_grand_branches(tree->root, ngb);
  if (ngb.empty()) {
    return;
  }
  Node* branch = rand_node(ngb);
  Node* left = branch->left;
  Node* right = branch->right;

  std::vector<double>& base = ctx.buf_base;
  std::vector<double>& w_left = ctx.buf_left;
  std::vector<double>& w_right = ctx.buf_right;
  ctx.make_base_children(branch, base, w_left, w_right);

  Target2 old_target(ctx, branch, base, w_left, w_right);
  double log_f_before = old_target.log_f(left->mu, right->mu);
  Laplace2 old_fit = old_target.laplace();
  double log_g_reverse = old_fit.log_dens(left->mu, right->mu);

  Node::Rule rule_old = branch->rule();

  // The predictor still reflects the old rule, so the base computed above is
  // the one to carry forward; only the child weights change, and they come out
  // of the split that the new rule performs anyway.
  std::vector<double>& wl_new = ctx.buf_left2;
  std::vector<double>& wr_new = ctx.buf_right2;
  branch->resample_rule(&wl_new, &wr_new);

  Target2 new_target(ctx, branch, base, wl_new, wr_new);
  Laplace2 new_fit = new_target.laplace();
  double proposal[2];
  new_fit.draw(proposal);
  double log_f_after = new_target.log_f(proposal[0], proposal[1]);
  double log_g_forward = new_fit.log_dens(proposal[0], proposal[1]);

  // Both children keep their depth, so the tree-shape prior cancels, and the
  // rule is drawn from its own prior, so the rule prior cancels as well. The
  // move is restricted to branches whose children are both leaves, so no
  // descendant's admissible cutpoint range is disturbed.
  double log_ratio = log_f_after - log_f_before + log_g_reverse -
    log_g_forward;

  if (std::isfinite(log_ratio) && std::log(unif_rand()) < log_ratio) {
    left->mu = proposal[0];
    right->mu = proposal[1];
    commit_children(branch, *ctx.eta, ctx.h, base, wl_new, wr_new, proposal[0],
                    proposal[1]);
  }
  else {
    branch->set_rule(rule_old);
    branch->split_support();
  }
}

// Independence Metropolis refresh of every leaf value, reusing the same Laplace
// fit as the tree moves. With soft rules the leaves are dependent, so this is a
// Metropolis-within-Gibbs sweep and the predictor is updated after each leaf.
//
// When the target is quadratic the fitted normal is the conditional posterior
// itself, so the acceptance probability is one and the draw is a Gibbs step. The
// two log-density evaluations that would compute that ratio are the single
// largest block of work in a Gaussian fit -- about two of the four passes over
// every observation that this sweep otherwise costs -- and skipping them is
// exact, not an approximation.
} // namespace

// One scalar parameter that enters the predictor with weight `wt` for the
// observations the node lists, drawn from its conditional posterior.
//
// This is the leaf refresh, and it is *also* the whole of a group-level random
// intercept: a random intercept is a scalar with a Gaussian prior entering the
// predictor with weight one for the observations in its group, which is a leaf
// with the gate removed. Exposing it means the quadratic closed form, the
// exponential form and the general Laplace-and-Metropolis path all serve both
// without a second implementation. `ctx.sigma_mu` carries whichever prior scale
// applies -- the forest's for a leaf, the group's for a random intercept.
void update_scalar(Node* node, Context& ctx) {
  double mu_ref = node->mu;

  if (ctx.quadratic) {
    // One pass: the fit is the conditional posterior, so there is no ratio to
    // evaluate and nothing else to ask of the target.
    Target1 target(ctx, node, nullptr, mu_ref);
    apply_leaf_delta(node, *ctx.eta, ctx.h, target.laplace().draw());
    return;
  }

  Laplace1 fit = fit_laplace1_at(ctx, node, mu_ref);

  double mean_forward = damped_mean(fit, mu_ref);
  double mu_new = mean_forward + fit.sd * norm_rand();

  double log_f_new;
  double log_f_old;
  ctx.log_f_pair_at(node, mu_ref, mu_new, mu_ref, &log_f_new, &log_f_old);

  // The reverse proposal is the same fit damped toward the proposed value.
  double mean_reverse = damped_mean(fit, mu_new);
  Laplace1 forward = {mean_forward, fit.sd};
  Laplace1 reverse = {mean_reverse, fit.sd};

  double log_ratio = log_f_new - log_f_old + reverse.log_dens(mu_ref) -
    forward.log_dens(mu_new);

  if (std::isfinite(log_ratio) && std::log(unif_rand()) < log_ratio) {
    apply_leaf_delta(node, *ctx.eta, ctx.h, mu_new);
  }
}

namespace {

void update_leaf_params(Tree* tree, Context& ctx) {
  std::vector<Node*>& leaf_list = ctx.buf_leaves;
  leaf_list.clear();
  leaves(tree->root, leaf_list);

  for (std::size_t l = 0; l < leaf_list.size(); l++) {
    update_scalar(leaf_list[l], ctx);
  }
}

// Metropolis update of one tree's gate bandwidth. Changing it moves every gate
// on every path at once, so the whole tree's contribution to the predictor has
// to be removed, rebuilt and compared.
// Optimal acceptance rate for a one-dimensional random walk.
const double BANDWIDTH_TARGET = 0.44;

void update_bandwidth(Tree* tree, Context& ctx, const Hypers& hypers) {

  // A tree with no splits has no gate, so the bandwidth does not enter the
  // likelihood at all and its full conditional is exactly the prior. Drawing
  // from it directly is both an exact Gibbs step and free, where the Metropolis
  // step below would rebuild every membership weight in the tree and evaluate
  // the whole likelihood to decide a move that cannot change the fit. Trees are
  // small under this prior, so this is a common case rather than an edge one.
  if (tree->root->is_leaf) {
    tree->bandwidth = exp_rand() * hypers.bandwidth_scale;
    return;
  }

  double bandwidth_old = tree->bandwidth;

  // Multiplicative random walk, whose asymmetry contributes log(new / old).
  double bandwidth_new = bandwidth_old *
    std::exp(tree->log_step * (2.0 * unif_rand() - 1.0));

  const int n_obs = ctx.family->N;
  const int stride = ctx.H;
  std::vector<double>& base = ctx.buf_bw_base;
  std::vector<double>& proposed = ctx.buf_bw_new;
  base.resize(n_obs);
  proposed.resize(n_obs);

  // The predictor with this tree taken out. Read straight from the row rather
  // than through an arma temporary: this is one traversal instead of the
  // strided extraction plus two vector copies it replaces, and the current
  // predictor never has to be saved, because it is `base` plus the tree under
  // the old bandwidth and rejecting simply leaves the row alone.
  const double* e = ctx.eta->memptr() + ctx.h;
  for (int i = 0; i < n_obs; i++) {
    base[i] = e[static_cast<std::size_t>(i) * stride];
  }

  arma::rowvec base_view(base.data(), n_obs, false, true);
  accumulate(tree->root, base_view, -1.0);

  ctx.bw_support.clear();
  save_support(tree->root, ctx.bw_support);

  tree->bandwidth = bandwidth_new;
  rebuild_support(tree->root);

  std::copy(base.begin(), base.end(), proposed.begin());
  arma::rowvec prop_view(proposed.data(), n_obs, false, true);
  accumulate(tree->root, prop_view, 1.0);

  // One pass for the difference the ratio needs, rather than two full
  // likelihood evaluations for two numbers that are then subtracted.
  double log_ratio = ctx.family->loglik_delta(*ctx.eta, ctx.h, proposed.data()) +
    Rf_dexp(bandwidth_new, hypers.bandwidth_scale, 1) -
    Rf_dexp(bandwidth_old, hypers.bandwidth_scale, 1) +
    std::log(bandwidth_new) - std::log(bandwidth_old);

  bool accept = std::isfinite(log_ratio) && std::log(unif_rand()) < log_ratio;

  if (hypers.adapt) {
    // Robbins-Monro on the log step with a vanishing gain. Adaptation runs only
    // during warmup, so the kernel used for the retained draws is fixed.
    tree->attempts++;
    double gain = 1.0 / std::sqrt(1.0 + tree->attempts / 50.0);
    tree->log_step += gain * ((accept ? 1.0 : 0.0) - BANDWIDTH_TARGET);
    if (tree->log_step < std::log(1.02)) {
      tree->log_step = std::log(1.02);
    }
    if (tree->log_step > std::log(100.0)) {
      tree->log_step = std::log(100.0);
    }
  }

  if (accept) {
    double* out = ctx.eta->memptr() + ctx.h;
    for (int i = 0; i < n_obs; i++) {
      out[static_cast<std::size_t>(i) * stride] = proposed[i];
    }
    return;
  }

  // The predictor was never written, so putting the old bandwidth and the
  // supports it implies back is the whole of the rollback -- and the supports
  // come from the snapshot rather than from evaluating every gate again.
  tree->bandwidth = bandwidth_old;
  std::size_t pos = 0;
  restore_support(tree->root, ctx.bw_support, pos);
}

void update_sigma_mu(Hypers& hypers, std::vector<Tree*>& forest) {
  std::vector<double> params;
  for (std::size_t t = 0; t < forest.size(); t++) {
    collect_leaf_params(forest[t]->root, params);
  }
  if (params.empty()) {
    return;
  }
  arma::vec mu(params.data(), params.size(), false);
  double prec_old = std::pow(hypers.sigma_mu, -2.0);
  double prec_new = half_cauchy_update_precision_mh(mu, prec_old,
                                                    hypers.scale_sigma_mu);
  hypers.sigma_mu = std::pow(prec_new, -0.5);
}

} // namespace

void update_forest(std::vector<Tree*>& forest, Context& ctx, Hypers& hypers) {

  for (std::size_t t = 0; t < forest.size(); t++) {
    Tree* tree = forest[t];
    ctx.sigma_mu = hypers.sigma_mu;

    if (tree->root->is_leaf) {
      node_birth(tree, ctx, hypers);
    }
    else {
      double u = unif_rand();
      if (u < 0.5 * P_BIRTH_DEATH) {
        node_birth(tree, ctx, hypers);
      }
      else if (u < P_BIRTH_DEATH) {
        node_death(tree, ctx, hypers);
      }
      else {
        change_rule(tree, ctx);
      }
    }

    update_leaf_params(tree, ctx);

    tree->sweeps++;

    if (hypers.soft && hypers.update_bandwidth &&
        tree->sweeps % hypers.bandwidth_every == 0) {
      update_bandwidth(tree, ctx, hypers);
    }
  }

  if (hypers.update_sigma_mu) {
    update_sigma_mu(hypers, forest);
    ctx.sigma_mu = hypers.sigma_mu;
  }

  if (hypers.update_s) {
    arma::uvec counts = arma::zeros<arma::uvec>(hypers.num_groups());
    for (std::size_t t = 0; t < forest.size(); t++) {
      get_var_counts(forest[t]->root, counts);
    }
    hypers.update_s_param(counts);
    if (hypers.update_alpha) {
      hypers.update_alpha_param();
    }
  }
}

} // namespace bartisan
