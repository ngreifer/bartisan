#ifndef GENBART_RANDOM_H
#define GENBART_RANDOM_H

#include <RcppArmadillo.h>
#include <memory>
#include <string>
#include <vector>
#include "hypers.h"
#include "mcmc.h"
#include "node.h"

namespace bartisan {

// A group-level random intercept, for one grouping factor and one additive
// predictor.
//
// The whole of it reuses the leaf machinery. A random intercept is a scalar with
// a Gaussian prior that enters the predictor with weight one for the
// observations in its level -- which is a leaf whose gate has been removed. So
// each level is represented by a Node carrying that level's observations, and
// updated by update_scalar(), which is the same function the leaf refresh calls.
// The quadratic closed form, the exponential form and the general
// Laplace-and-Metropolis path therefore all apply without a second
// implementation of any of them.
//
// The Hypers here exists only to carry `soft = false` to the exponential-form
// check: a random intercept's weights really are one whatever the decision rules
// are, so the exponential form is available even in a soft-rule fit, and it would
// not be if the level nodes pointed at the forest's own hyperparameters.
struct RandomTerm {
  std::string label;
  int num_levels;

  // The scale of the intercepts, and its half-Cauchy prior.
  double tau;
  double tau_scale;
  bool update_tau;

  std::unique_ptr<Hypers> hypers;
  std::unique_ptr<Tree> tree;
  std::vector<Node*> levels;

  RandomTerm(const std::string& label_, const arma::ivec& level_of,
             int num_levels_, double tau_start, double tau_scale_,
             bool update_tau_, const arma::mat* X);

  ~RandomTerm();

  // Draw every level's intercept, then the scale they share.
  void update(Context& ctx);

  arma::vec values() const;
};

// Every random-effect term, for every additive predictor. `terms[h][r]` is the
// r-th grouping factor's intercepts for predictor h: each predictor gets its own
// set, so a family with several of them -- a zero-inflated count model, say --
// has a group effect on each part rather than one shared between them.
struct RandomEffects {
  std::vector<std::vector<std::unique_ptr<RandomTerm> > > terms;

  bool empty() const { return terms.empty(); }

  // The intercepts' contribution to the predictor, added into eta. Called once
  // at the start so that eta and the intercepts agree before the first sweep.
  void accumulate_into(arma::mat& eta) const;

  void update(Context& ctx, int h);
};

// Build the whole structure from what the R side passes down. `spec` is a list
// with one element per grouping term, each a list of `label`, `levels` (the
// zero-based level of every observation) and `num_levels`.
std::unique_ptr<RandomEffects> make_random_effects(const Rcpp::List& spec, int H,
                                                   double tau_scale,
                                                   bool update_tau,
                                                   const arma::mat* X);

} // namespace bartisan

#endif
