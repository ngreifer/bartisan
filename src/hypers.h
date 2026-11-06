#ifndef GENBART_HYPERS_H
#define GENBART_HYPERS_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace bartisan {

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

  // The flags below were built to locate an apparent narrowing of soft-rule
  // posteriors in simulation-based calibration, which on fresh seed blocks
  // turned out to be one block's chance rather than a defect (`_dev/TASKS.md`,
  // 2026-11-05). They stay as tools for testing the sampler.
  //
  // Internal, off by default, and read only by `update_forest()`: propose the
  // bandwidth with the tree's leaves integrated out where the target is
  // quadratic, and with `collapsed_check` verify each marginal likelihood
  // against the family's exact log density. `the$collapsed_bandwidth` in R sets
  // it.
  bool collapsed_bandwidth = false;
  bool collapsed_check = false;
  // Internal and off by default: attempt only births and deaths, never a change,
  // with the move-choice terms of their acceptance ratios to match. Birth and
  // death alone are a valid sampler. `the$no_change_move` in R sets it.
  bool no_change_move = false;
  // Internal and off by default: no birth, death or change at all, so every
  // tree keeps the structure it started with while its leaves and bandwidth are
  // drawn. With the trees started at the generator's own (`fixed_trees` in the
  // control list), an SBC of this sampler tests everything but the structure
  // moves (`_dev/TASKS.md`, 2026-11-05). Set by `_dev/sbc.R` appending
  // `fixed_structure` to the control list; there is no `the$` route.
  bool fixed_structure = false;

  // Which gate the soft rules use; see GateShape in node.h. Ignored when the
  // rules are hard.
  int gate;

  // True during warmup only, while the bandwidth proposal is being tuned.
  bool adapt;

  // `split_prior_` empty leaves the splitting proportions at the uniform prior
  // and lets `update_s_` decide whether they are drawn. Non-empty fixes them at
  // the caller's weights, which is what `split_prior` in `bartisan_control()`
  // does, and then nothing draws them.
  //
  // `split_mask_` is a different statement and must not be confused with it: a
  // zero says this forest may never split on that group, because its own formula
  // does not name it, and says nothing about the weights of the groups that are
  // left. The proportions are still drawn when asked for, over the allowed
  // groups alone. Empty allows every group.
  Hypers(const arma::sp_mat& group_probs, double sigma_mu_, double gamma_,
         double beta_, double alpha_, double alpha_scale_, double alpha_shape_1_,
         double alpha_shape_2_, bool update_sigma_mu_, bool update_s_,
         bool update_alpha_, bool soft_, double bandwidth_scale_,
         bool update_bandwidth_, int bandwidth_every_, int gate_,
         const arma::vec& split_prior_ = arma::vec(),
         const arma::vec& split_mask_ = arma::vec());

  // Draw a (group, variable) pair: the group from s, then the variable from
  // that group's column of group_probs. Grouping lets the dummy columns of one
  // factor share a single sparsity weight.
  arma::uvec sample_var() const;

  void update_alpha_param();
  void update_s_param(const arma::uvec& counts);

  // Take the proportions and the concentration from another forest's draw,
  // which is how `share_sparsity` gives every forest the same ones: the counts
  // are pooled, one forest draws, and the rest copy. Only meaningful when the
  // two carry the same allowed set, which `same_allowed()` is for.
  void copy_s_from(const Hypers& other);
  bool same_allowed(const Hypers& other) const;

  int num_groups() const { return num_groups_; }
  int num_vars() const { return static_cast<int>(group_probs_.n_rows); }
  const arma::vec& log_s() const { return log_s_; }

private:
  arma::vec s_;
  arma::vec log_s_;
  int num_groups_;

  // The groups this forest may split on, as indices into s_. Every group unless
  // a mask said otherwise. The Dirichlet is over these and s_ is zero off them,
  // so `sample_var()` cannot propose one that is not here.
  arma::uvec allowed_;

  arma::sp_mat group_probs_;
};

inline double grow_prob(const Hypers* hypers, int depth) {
  return hypers->gamma * std::pow(1.0 + static_cast<double>(depth), -hypers->beta);
}

} // namespace bartisan

#endif
