#include "hypers.h"
#include "slice.h"

namespace bartisan {

Hypers::Hypers(const arma::sp_mat& group_probs, double sigma_mu_, double gamma_,
               double beta_, double alpha_, double alpha_scale_,
               double alpha_shape_1_, double alpha_shape_2_,
               bool update_sigma_mu_, bool update_s_, bool update_alpha_,
               bool soft_, double bandwidth_scale_, bool update_bandwidth_,
               int bandwidth_every_, int gate_,
               const arma::vec& split_prior_,
               const arma::vec& split_mask_) {

  group_probs_ = group_probs;
  num_groups_ = static_cast<int>(group_probs_.n_cols);

  bool fixed_s = split_prior_.n_elem == static_cast<arma::uword>(num_groups_);

  // A mask of the wrong length, or one that allows nothing, is read as allowing
  // everything: a forest that may split on no group at all is a pinned forest,
  // which never splits anyway, and giving it an empty support would leave the
  // Dirichlet below dividing by zero.
  if (split_mask_.n_elem == static_cast<arma::uword>(num_groups_)) {
    allowed_ = arma::find(split_mask_ > 0.0);
  }

  if (allowed_.n_elem == 0) {
    allowed_ = arma::regspace<arma::uvec>(0, num_groups_ - 1);
  }

  if (fixed_s) {
    s_ = split_prior_ / arma::accu(split_prior_);
  }
  else {
    // Uniform over the groups this forest may use, zero on the rest.
    s_ = arma::zeros<arma::vec>(num_groups_);
    s_.elem(allowed_).fill(1.0 / static_cast<double>(allowed_.n_elem));
  }
  log_s_ = arma::log(s_);

  gamma = gamma_;
  beta = beta_;

  alpha = alpha_;
  alpha_scale = alpha_scale_ > 0.0
    ? alpha_scale_
    : static_cast<double>(allowed_.n_elem);
  alpha_shape_1 = alpha_shape_1_;
  alpha_shape_2 = alpha_shape_2_;
  update_s = update_s_;
  update_alpha = update_alpha_;

  sigma_mu = sigma_mu_;
  scale_sigma_mu = sigma_mu_;
  update_sigma_mu = update_sigma_mu_;

  soft = soft_;
  bandwidth_scale = bandwidth_scale_;
  update_bandwidth = update_bandwidth_;
  bandwidth_every = bandwidth_every_;
  gate = gate_;
  adapt = false;

  // A single group carries no sparsity information, so leave s at its prior.
  // Counted over the groups this forest may actually use: a forest held to one
  // predictor has nothing to select between, whatever the model as a whole has.
  if (allowed_.n_elem < 2) {
    update_s = false;
    update_alpha = false;
  }

  // Weights the caller supplied are a statement about the predictors, not a
  // starting point for one, so nothing draws over them. `bartisan_control()`
  // already turns the sparsity prior off when `split_prior` is given; this makes
  // the guarantee hold whatever reaches the constructor.
  if (fixed_s) {
    update_s = false;
    update_alpha = false;
  }
}

arma::uvec Hypers::sample_var() const {
  arma::uvec group_var(2);
  int group = sample_class(s_);
  group_var(0) = static_cast<arma::uword>(group);
  group_var(1) = static_cast<arma::uword>(sample_class_col(group_probs_, group));
  return group_var;
}

// Sample s from its Dirichlet full conditional. The shape parameters are often
// small, so draw log-gamma variates directly and normalize with log-sum-exp
// rather than forming gamma draws that would underflow to zero.
void Hypers::update_s_param(const arma::uvec& counts) {
  // Over the allowed groups only, so that a group this forest may not split on
  // stays at zero instead of being handed prior mass it could never earn back
  // through a count. The concentration is alpha spread over the groups in play,
  // which is what `update_alpha_param()` assumes as well.
  arma::uword k = allowed_.n_elem;

  arma::vec shape_up = alpha / static_cast<double>(k) *
    arma::ones<arma::vec>(k);
  shape_up += arma::conv_to<arma::vec>::from(counts.elem(allowed_));

  arma::vec logs(k);
  for (arma::uword i = 0; i < k; i++) {
    logs(i) = rlgam(shape_up(i));
  }
  logs -= log_sum_exp(logs);

  log_s_.fill(R_NegInf);
  log_s_.elem(allowed_) = logs;

  s_ = arma::zeros<arma::vec>(num_groups_);
  s_.elem(allowed_) = arma::exp(logs);
}

void Hypers::copy_s_from(const Hypers& other) {
  s_ = other.s_;
  log_s_ = other.log_s_;
  alpha = other.alpha;
}

bool Hypers::same_allowed(const Hypers& other) const {
  return allowed_.n_elem == other.allowed_.n_elem &&
    arma::all(allowed_ == other.allowed_);
}

void Hypers::update_alpha_param() {
  // The allowed groups are the ones the Dirichlet is over, so they are the ones
  // its concentration is estimated from. Averaging over all of them would
  // average in the negative infinities sitting on the masked entries.
  double mean_log_s = arma::mean(log_s_.elem(allowed_));
  double p = static_cast<double>(allowed_.n_elem);
  double scale = alpha_scale;
  double shape_1 = alpha_shape_1;
  double shape_2 = alpha_shape_2;

  auto logf = [mean_log_s, p, scale, shape_1, shape_2](double rho) {
    if (!(rho > 0.0) || !(rho < 1.0)) {
      return R_NegInf;
    }
    double a = scale * rho / (1.0 - rho);
    double log_prior = (shape_1 - 1.0) * std::log(rho) +
      (shape_2 - 1.0) * std::log1p(-rho);
    return a * mean_log_s + R::lgammafn(a) - p * R::lgammafn(a / p) + log_prior;
  };

  double rho_current = alpha / (alpha + alpha_scale);
  double rho_up = slice_sampler(rho_current, logf, 0.1, 0.0, 1.0);
  alpha = alpha_scale * rho_up / (1.0 - rho_up);
}

} // namespace bartisan
