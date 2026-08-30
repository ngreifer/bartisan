#ifndef GENBART_MCMC_H
#define GENBART_MCMC_H

#include <RcppArmadillo.h>
#include <vector>
#include "family.h"
#include "node.h"

namespace bartisan {

// Everything the leaf-level Laplace fits need: the family, the current
// predictors for all observations, and which additive predictor is being
// updated.
struct Context {
  Family* family;
  arma::mat* eta;   // H x N, so that eta->colptr(i) is observation i
  int H;
  int h;
  double sigma_mu;

  // The shape of the family's log density in the predictor currently being
  // updated, and whether that shape is the quadratic one. Both are set whenever
  // `h` changes; see Family::target_form().
  TargetForm form;
  bool quadratic;

  // Whether to take the closed forms a quadratic target allows. Settable so the
  // test suite can check that they reproduce the general path.
  bool exact_quadratic;

  // Whether to accumulate a leaf's sums through the vtable rather than through
  // the family's own statically dispatched loop. Settable for the same reason:
  // so the test suite can check the two agree.
  bool generic_accumulate;

  // Whether to hand the family a whole leaf at a time instead of one
  // observation at a time. Set for families that pay a fixed cost per call --
  // one whose log density is an R function -- and settable directly so that the
  // test suite can check the two paths agree.
  bool blocked;

  Context(Family* family_, arma::mat* eta_, bool blocked_ = false)
    : family(family_), eta(eta_), H(family_->H), h(0), sigma_mu(1.0),
      form(family_->target_form(0)),
      quadratic(family_->is_quadratic(0)), exact_quadratic(true),
      generic_accumulate(false),
      blocked(blocked_ || family_->wants_block()) {
    scratch.resize(family_->H);
  }

  // Log of Linero's F: the leaf prior times the likelihood of the observations
  // this node supports, as a function of the node's parameter. Observations
  // outside the node's support contribute identically before and after any
  // move here and so cancel from every acceptance ratio.
  double log_f(const Node* node, const std::vector<double>& base,
               double mu) const;

  // Score and information of the same target, both including the prior, in a
  // single pass over the node's support. Fisher scoring always wants them at
  // the same point, so computing them separately doubled both the loop and the
  // number of family evaluations.
  void score_info(const Node* node, const std::vector<double>& base, double mu,
                  double* score, double* info) const;

  // The target at two values of the leaf parameter, in one pass. The leaf
  // Metropolis step needs exactly this pair.
  void log_f_pair(const Node* node, const std::vector<double>& base,
                  double mu_a, double mu_b, double* out_a, double* out_b) const;

  // The same two quantities without materializing a base array. The predictor
  // already contains this leaf's contribution at mu_ref, so the shifted value
  // is eta + wt * (mu - mu_ref). This is the hot path -- every leaf of every
  // tree on every sweep -- and it saves one pass and one temporary per visit.
  void score_info_at(const Node* node, double mu_ref, double mu, double* score,
                     double* info) const;
  void log_f_pair_at(const Node* node, double mu_ref, double mu_a, double mu_b,
                     double* out_a, double* out_b) const;

  // The two-child versions, where the children's weights are the parent's
  // weights split by the gate. The cross term of the information vanishes for
  // hard rules, which recovers a pair of independent fits.
  double log_f2(const Node* parent, const std::vector<double>& base,
                const std::vector<double>& w_left,
                const std::vector<double>& w_right,
                double mu_left, double mu_right) const;

  void score_info2(const Node* parent, const std::vector<double>& base,
                   const std::vector<double>& w_left,
                   const std::vector<double>& w_right,
                   double mu_left, double mu_right,
                   double* g, double* info) const;

  // The log target, its score and its information at one point, in a single
  // pass. Everything a quadratic target needs comes from one of these.
  void log_f_score_info(const Node* node, const std::vector<double>& base,
                        double mu, double* f, double* s, double* j) const;
  void log_f_score_info_at(const Node* node, double mu_ref, double mu, double* f,
                           double* s, double* j) const;
  void log_f_score_info2(const Node* parent, const std::vector<double>& base,
                         const std::vector<double>& w_left,
                         const std::vector<double>& w_right, double mu_left,
                         double mu_right, double* f, double* g,
                         double* info) const;

  // The predictor for this node's supported observations with the node's own
  // contribution removed.
  void make_base(const Node* node, std::vector<double>& base) const;

  // The same, removing the contribution of both children of a branch, together
  // with the weights the branch's current rule assigns to them. Valid only
  // while the predictor still reflects that rule and those leaf values.
  void make_base_children(const Node* parent, std::vector<double>& base,
                          std::vector<double>& w_left,
                          std::vector<double>& w_right) const;

  // Reused across moves so that a sweep does not allocate once per proposal.
  std::vector<double> buf_base;
  std::vector<double> buf_left;
  std::vector<double> buf_right;

  // A second pair, for the change move, which needs the weights the old rule
  // gave and the weights the new one gives at the same time. These were local
  // vectors, allocated and freed on every change proposal.
  std::vector<double> buf_left2;
  std::vector<double> buf_right2;

  // Working space for the bandwidth move, which needs two full-length vectors:
  // the predictor with the tree removed, and the same with the tree put back
  // under the proposed bandwidth.
  std::vector<double> buf_bw_base;
  std::vector<double> buf_bw_new;

  // The supports as they stood before a bandwidth proposal, so that rejecting it
  // -- which happens rather more often than not -- costs a copy rather than a
  // second evaluation of every gate in the tree.
  SupportStore bw_support;

  // The node lists a move picks from. Returning these by value allocated and
  // freed a vector on every proposal of every tree of every sweep; the buffers
  // are cleared and refilled instead. Only one is live at a time, but they are
  // kept separate so that stays true by construction rather than by reading.
  std::vector<Node*> buf_leaves;
  std::vector<Node*> buf_branches;

private:
  mutable std::vector<double> scratch;

  // Working space for the block evaluations. `block` holds the predictors of a
  // node's support in the family's H-by-n layout, sized for two leaf values at
  // once so that the paired evaluations need a single call.
  mutable std::vector<double> block;
  mutable std::vector<double> values;
  mutable std::vector<double> d1;
  mutable std::vector<double> d2;
  mutable std::vector<int> pair_idx;

  // Copy the node's current predictors into `block`, overwriting component h
  // with the leaf contribution implied by `mu`. `at` is the observation the
  // block starts at: 0 for a single value, n for the second half of a paired
  // block.
  void fill(const Node* node, const std::vector<double>& base, double mu,
            std::size_t at) const;
  void fill_at(const Node* node, double shift, std::size_t at) const;
  void fill2(const Node* parent, const std::vector<double>& base,
             const std::vector<double>& w_left,
             const std::vector<double>& w_right, double mu_left,
             double mu_right, std::size_t at) const;
  void make_pair_idx(const Node* node) const;

  // Hand the filled block to the family, statically or through the vtable.
  void accumulate1(const Node* node, int n, double* f, double* s,
                   double* j) const;

  // The same three sums without materializing the block: the family reads each
  // observation's predictors from `eta` and shifts component h by `wt * shift`.
  void accumulate1_at(const Node* node, int n, double shift, double* f,
                      double* s, double* j) const;

  // The same, forming component h from a base the caller has already computed:
  // the node's predictor with its own contribution taken out.
  void accumulate1_from(const Node* node, int n, const std::vector<double>& base,
                        double mu, double* f, double* s, double* j) const;

  void accumulate2_from(const Node* parent, int n,
                        const std::vector<double>& base,
                        const std::vector<double>& w_left,
                        const std::vector<double>& w_right, double mu_left,
                        double mu_right, double* f, double* g,
                        double* info) const;

  void accumulate2(const Node* parent, int n,
                   const std::vector<double>& w_left,
                   const std::vector<double>& w_right, double* f, double* g,
                   double* info) const;

  // Size the block buffers for `copies` evaluations of an n-observation node.
  void reserve(std::size_t n, std::size_t copies) const;
};

// A univariate Gaussian approximation to a node's conditional posterior.
struct Laplace1 {
  double mean;
  double sd;
  double log_dens(double x) const { return R::dnorm4(x, mean, sd, 1); }
  double draw() const { return mean + sd * norm_rand(); }
};

// A bivariate Gaussian approximation, stored by its precision so that no
// inversion is needed to evaluate the density.
struct Laplace2 {
  double mean[2];
  double prec[3];   // (0,0), (0,1), (1,1)
  double log_det_prec;

  double log_dens(double x, double y) const;
  void draw(double* out) const;
};

Laplace1 fit_laplace1(Context& ctx, const Node* node,
                      const std::vector<double>& base);

Laplace2 fit_laplace2(Context& ctx, const Node* parent,
                      const std::vector<double>& base,
                      const std::vector<double>& w_left,
                      const std::vector<double>& w_right);

// The log target as a function of one leaf value, over the observations a node
// supports.
//
// A leaf value enters the additive predictor linearly, so when the family's log
// density is quadratic in the predictor the target is quadratic in the leaf
// value -- and so is the leaf prior, which is Gaussian. One pass over the node
// then determines the whole function: the expansion below is an identity, not a
// Taylor approximation. The Laplace "approximation" is the conditional posterior
// exactly, and the log target at any value is arithmetic.
//
// That is where a conjugate sampler's advantage comes from, and it is worth
// naming precisely: a birth move needs the log target at two values and two
// Laplace fits, which is six passes over the node, and for a quadratic family
// two suffice. Where the family is not quadratic this class forwards to exactly
// what the sampler did before -- a pass per value and a Fisher-scoring loop per
// fit -- so nothing about those families changes.
class Target1 {
public:
  // `base` is the predictor with this node's contribution removed, or null to
  // take the node's current contribution from the predictor and shift it, in
  // which case `ref` is the leaf value the predictor already holds.
  Target1(Context& ctx, const Node* node, const std::vector<double>* base,
          double ref);

  double log_f(double mu) const;
  Laplace1 laplace() const;

private:
  // How much of the work the shape of the target lets us skip. CLOSED is the
  // quadratic case, where the whole function follows from one pass and the fit
  // is a formula. EXPONENTIAL also follows from one pass, but the mode has to be
  // found by iteration -- on three scalars, so no further passes. PASSES is
  // everything else, which is what the sampler did throughout before any of
  // this.
  enum Mode { PASSES, CLOSED, EXPONENTIAL };

  double exp_score(double mu) const;
  double exp_info(double mu) const;

  Context* ctx_;
  const Node* node_;
  const std::vector<double>* base_;
  double ref_;
  Mode mode_;

  // The quadratic case: value, score and information at ref_.
  double f_;
  double d1_;
  double d2_;

  // The exponential case: the likelihood written as c + a * mu + b * exp(r * mu),
  // with the coefficients recovered from the same single pass. Reference-free,
  // so the fit does not depend on where the expansion was taken.
  double rate_;   // the rate r in b * exp(r * mu)
  double a_;
  double b_;
  double c_;
};

// The same for the pair of leaf values a branch's children carry, which enter
// the predictor linearly through the gate weights.
class Target2 {
public:
  Target2(Context& ctx, const Node* parent, const std::vector<double>& base,
          const std::vector<double>& w_left,
          const std::vector<double>& w_right);

  double log_f(double mu_left, double mu_right) const;
  Laplace2 laplace() const;

private:
  enum Mode { PASSES, CLOSED, EXPONENTIAL };

  Context* ctx_;
  const Node* parent_;
  const std::vector<double>* base_;
  const std::vector<double>* w_left_;
  const std::vector<double>* w_right_;
  Mode mode_;
  double f_;
  double g_[2];
  double info_[3];

  // The exponential case. With hard rules each observation reaches exactly one
  // child, so the target separates: the cross curvature is exactly zero and the
  // two children are two independent one-dimensional problems.
  double rate_;   // the rate r in b * exp(r * mu)
  double a_[2];
  double b_[2];
  double c_;
};

// One scalar parameter entering the predictor over the observations a node
// lists, drawn from its conditional posterior. This is the leaf refresh, and a
// group-level random intercept is the same problem with the gate removed; see
// the comment on the definition.
void update_scalar(Node* node, Context& ctx);

// One sweep over the trees of one forest, then the forest's hyperparameters.
void update_forest(std::vector<Tree*>& forest, Context& ctx, Hypers& hypers);

} // namespace bartisan

#endif
