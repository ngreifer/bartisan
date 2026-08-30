#include "hypers.h"
#include "slice.h"

namespace bartisan {

Hypers::Hypers(const arma::sp_mat& group_probs, double sigma_mu_, double gamma_,
               double beta_, double alpha_, double alpha_scale_,
               double alpha_shape_1_, double alpha_shape_2_,
               bool update_sigma_mu_, bool update_s_, bool update_alpha_,
               bool soft_, double bandwidth_scale_, bool update_bandwidth_,
               int bandwidth_every_, int gate_) {

  group_probs_ = group_probs;
  num_groups_ = static_cast<int>(group_probs_.n_cols);
  s_ = arma::ones<arma::vec>(num_groups_) / static_cast<double>(num_groups_);
  log_s_ = arma::log(s_);

  gamma = gamma_;
  beta = beta_;

  alpha = alpha_;
  alpha_scale = alpha_scale_ > 0.0 ? alpha_scale_
                                  : static_cast<double>(num_groups_);
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
  if (num_groups_ < 2) {
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
  arma::vec shape_up = alpha / static_cast<double>(num_groups_) *
    arma::ones<arma::vec>(num_groups_);
  shape_up += arma::conv_to<arma::vec>::from(counts);

  arma::vec logs(num_groups_);
  for (int i = 0; i < num_groups_; i++) {
    logs(i) = rlgam(shape_up(i));
  }
  logs -= log_sum_exp(logs);

  log_s_ = logs;
  s_ = arma::exp(logs);
}

void Hypers::update_alpha_param() {
  double mean_log_s = arma::mean(log_s_);
  double p = static_cast<double>(num_groups_);
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
