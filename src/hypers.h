#ifndef GENBART_HYPERS_H
#define GENBART_HYPERS_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace genbart {

// Hyperparameters shared by every tree of one forest. A model with H additive
// predictors carries H of these, so that each predictor gets its own sparsity
// pattern and its own leaf scale.
struct Hypers {

  // Branching probability at depth d is gamma * (1 + d)^(-beta).
  double gamma;
  double beta;

  // Dirichlet(alpha / G, ..., alpha / G) prior on the splitting proportions s
  // over the G predictor groups, with a Beta prior on alpha / (alpha + scale).
  double alpha;
  double alpha_scale;
  double alpha_shape_1;
  double alpha_shape_2;
  bool update_s;
  bool update_alpha;

  // Half-Cauchy(scale_sigma_mu) prior on the leaf standard deviation.
  double sigma_mu;
  double scale_sigma_mu;
  bool update_sigma_mu;

  // Soft decision rules. The bandwidth is on the scale of the predictors, which
  // are mapped to [0, 1] before fitting, and is given the exponential prior of
  // Linero and Yang (2018) with mean bandwidth_scale.
  bool soft;
  double bandwidth_scale;
  bool update_bandwidth;

  // How many sweeps between bandwidth attempts for a given tree. The bandwidth
  // is one scalar per tree with an adaptive random-walk proposal, so attempting
  // it on every sweep buys mixing it does not need while costing a full rebuild
  // of the tree's supports each time -- measured as the single largest item in a
  // soft-rule fit.
  int bandwidth_every;

  // Which gate the soft rules use; see GateShape in node.h. Ignored when the
  // rules are hard.
  int gate;

  // True during warmup only, while the bandwidth proposal is being tuned.
  bool adapt;

  Hypers(const arma::sp_mat& group_probs, double sigma_mu_, double gamma_,
         double beta_, double alpha_, double alpha_scale_, double alpha_shape_1_,
         double alpha_shape_2_, bool update_sigma_mu_, bool update_s_,
         bool update_alpha_, bool soft_, double bandwidth_scale_,
         bool update_bandwidth_, int bandwidth_every_, int gate_);

  // Draw a (group, variable) pair: the group from s, then the variable from
  // that group's column of group_probs. Grouping lets the dummy columns of one
  // factor share a single sparsity weight.
  arma::uvec sample_var() const;

  void update_alpha_param();
  void update_s_param(const arma::uvec& counts);

  int num_groups() const { return num_groups_; }
  int num_vars() const { return static_cast<int>(group_probs_.n_rows); }
  const arma::vec& log_s() const { return log_s_; }

private:
  arma::vec s_;
  arma::vec log_s_;
  int num_groups_;
  arma::sp_mat group_probs_;
};

inline double grow_prob(const Hypers* hypers, int depth) {
  return hypers->gamma * std::pow(1.0 + static_cast<double>(depth), -hypers->beta);
}

} // namespace genbart

#endif
