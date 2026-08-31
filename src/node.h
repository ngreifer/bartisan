#ifndef GENBART_NODE_H
#define GENBART_NODE_H

#include <RcppArmadillo.h>
#include <vector>
#include <cstdint>
#include "hypers.h"

namespace bartisan {

struct Node;

// What a splitting rule does with a missing value. This is missingness
// incorporated in attributes (Twala, Jones & Hand 2008; for BART, Kapelner &
// Bleich 2015): the rule itself carries the answer, and it is drawn from its
// prior alongside the variable and the cutpoint rather than being imputed or
// decided in advance. NA_ONLY is the third of Twala's rules, splitting on
// whether the value is there at all and ignoring the cutpoint.
enum NaRule { NA_LEFT = 0, NA_RIGHT = 1, NA_ONLY = 2 };

const int NUM_NA_RULES = 3;

// The shape of a soft decision rule's gate, as a function of the standardized
// distance from the cutpoint. Both are cumulative distribution functions, so
// both are monotone, both hit 1/2 at the cutpoint, and under both the two
// children's weights sum to the parent's.
//
// GATE_LOGISTIC is the logistic CDF of Linero and Yang (2018). It has unbounded
// support, and at the default bandwidth that is not a technicality: dropping a
// weight needs it below 1e-10, which for the logistic needs the observation
// about 23 bandwidths from the cutpoint -- further than the whole unit interval.
// So every observation reaches every leaf, a soft tree costs about two and a
// half times a hard one, and each gate costs an exp().
//
// GATE_SMOOTHSTEP is the Beta(2, 2) CDF, the smoothstep polynomial
// t^2 (3 - 2t), on a bounded interval around the cutpoint. Outside it the
// answer is exactly zero or one: the observation takes one side only, the
// subtree on the other side is never visited, and no transcendental function is
// evaluated. It is continuously differentiable, which is what soft rules are
// for, but not twice.
// Each gate is the CDF of a symmetric Beta kernel, so the choice is a choice of
// how many derivatives the fit has and how wide the kernel must be to smooth by
// a given amount:
//
//   smoothstep    Beta(2, 2)   one derivative     support 2.24 sd
//   smootherstep  Beta(3, 3)   two derivatives    support 2.65 sd
//   logistic                   infinitely many    unbounded
//
// Measured, the support width turns out to matter far less than expected and the
// choice is close to a free one: at a bandwidth wide enough that a compact gate
// truncates nothing at all, it is still 1.45 times faster than the logistic. What
// the compact gates save is the exp(), not the work on the far side of the
// cutpoint. So the polynomial gates all cost about the same and no gate is faster
// than smoothstep for a reason worth having; see TASKS.md for the measurements,
// and for two candidates rejected -- a raised cosine, which is in smoothstep's
// smoothness class but needs a wider kernel and a cos(), and the Beta(1, 1)
// linear ramp, which was the narrowest and the fastest by about a twentieth but
// gives up a differentiable fit, which is what soft rules are for.
enum GateShape {
  GATE_LOGISTIC = 0,
  GATE_SMOOTHSTEP = 1,
  GATE_SMOOTHERSTEP = 2
};

// Half-width of each compactly supported gate, as a multiple of the bandwidth.
// These are what make `bandwidth` mean the same amount of smoothing under every
// gate: the logistic with scale s has standard deviation s * pi / sqrt(3), and
// the Beta(a, a) density on [-b, b] has standard deviation b / sqrt(2a + 1), so
// the kernels have the same spread when b = s * pi * sqrt((2a + 1) / 3).
const double SMOOTHSTEP_HALF_WIDTH = 4.055935661788187;    // a = 2: sqrt(5/3)
const double SMOOTHERSTEP_HALF_WIDTH = 4.797928400597175;  // a = 3: sqrt(7/3)

// Per-tree state. The bandwidth of the soft decision rules is a tree-level
// parameter, as in Linero and Yang (2018), so it lives here rather than on
// each node.
struct Tree {
  Node* root;
  double bandwidth;
  Hypers* hypers;
  const arma::mat* X;

  // Which columns of X contain a missing value. A rule on a column that has
  // none needs no answer for the missing case, so it is not given one and the
  // prior is exactly what it was before missing data was supported.
  const arma::uvec* has_na;

  // The level of each observation within each predictor group whose columns are
  // mutually exclusive indicators, so that a rule on such a group can name a
  // *subset* of its levels rather than a threshold on one indicator. `codes` is
  // observations by categorical groups and holds -1 where the group is missing;
  // `cat_col` maps a group to its column of `codes`, or -1; `n_levels` maps a
  // group to its number of levels, or 0.
  //
  // A threshold on one indicator can only peel a single level off the rest, so
  // it reaches 2^K - K of the B_K partitions of K levels and leaves the bulk
  // undivided. Deshpande (2024) is the reference; a subset rule reaches all of
  // them.
  const arma::imat* codes;
  const arma::ivec* cat_col;
  const arma::ivec* n_levels;

  int levels_of(int group) const {
    return n_levels == nullptr ? 0 : (*n_levels)(group);
  }

  int code_col(int group) const { return (*cat_col)(group); }

  // State for the adaptive bandwidth proposal. The step is the half-width of
  // the multiplicative random walk on the log scale, tuned during warmup
  // towards the acceptance rate that is optimal for a one-dimensional random
  // walk, then frozen so the sampling phase uses a fixed kernel.
  double log_step;
  long attempts;

  // Sweeps this tree has seen, so that the bandwidth move can be attempted on
  // every `bandwidth_every`-th one.
  long sweeps;

  // Detached nodes, kept for reuse. A birth allocates two children and is
  // rejected about two thirds of the time, so without this the sampler spends a
  // measurable share of its time in the allocator -- and a recycled node's index
  // and weight vectors keep their capacity, so a support does not have to grow
  // them again either.
  std::vector<Node*> pool;

  Tree(Hypers* hypers_, const arma::mat* X_, const arma::uvec* has_na_,
       const arma::imat* codes_ = nullptr, const arma::ivec* cat_col_ = nullptr,
       const arma::ivec* n_levels_ = nullptr);
  ~Tree();

  // A fresh child of `parent`, from the pool if one is waiting. Identical in
  // every field to `new Node(parent)`; both go through Node::init_as_child(), so
  // the two cannot drift apart.
  Node* take_node(Node* parent);

  // Hand a detached leaf back. It must have no children of its own.
  void give_node(Node* node);

  bool splits_on_missing(int var) const {
    return has_na != nullptr && (*has_na)(var) > 0;
  }
};

// A node of a decision tree. Each node caches the observations that reach it
// together with their membership weights, which is what lets the soft and hard
// cases share one code path: for a hard tree every weight is one, and for a
// soft tree the weight is the product of the gates along the path from the
// root.
//
// Weights below WEIGHT_TOL are dropped and the subtree below them is not
// explored. Those observations contribute the same amount to the predictor
// before and after any move at this node, so they cancel from every acceptance
// ratio; dropping them is what keeps soft trees close to the cost of hard ones.
struct Node {

  static constexpr double WEIGHT_TOL = 1e-10;

  Node* left;
  Node* right;
  Node* parent;
  Tree* tree;

  bool is_root;
  bool is_leaf;
  int depth;

  int var;
  int group;
  double val;
  double lower;
  double upper;
  NaRule na_rule;

  // Non-empty exactly when this node's rule names a subset of a categorical
  // group's levels, in which case a set bit is a level that goes left and
  // `val`, `lower` and `upper` are unused. Empty for a numeric rule, so the
  // common case allocates nothing, and a node from the pool keeps its capacity.
  std::vector<std::uint32_t> mask;

  bool is_categorical() const { return !mask.empty(); }

  double mu;

  // Center and standard deviation of the Laplace approximation to this node's
  // conditional posterior.
  double mu_star;
  double v_star;

  std::vector<int> idx;

  // The membership weight of each observation in `idx`: the product of the gates
  // along the path from the root. **Empty when the rules are hard**, because
  // then a rule sends an observation entirely one way and every weight is
  // exactly one. Read it through weights(), which returns a null pointer in that
  // case, and treat null as all ones -- multiplying by 1.0 is exact, so the two
  // paths agree to the last bit, and the loops that dominate the sampler stop
  // loading a vector of ones and multiplying by it.
  std::vector<double> wt;

  const double* weights() const { return wt.empty() ? nullptr : wt.data(); }

  Node(Tree* tree_, int n_obs);
  Node(Node* parent_);
  ~Node();

  // Every field a child node starts life with, in one place, so that a recycled
  // node and a newly allocated one are indistinguishable.
  void init_as_child(Node* parent_);

  bool is_left() const;
  void get_limits();

  // Everything a splitting rule consists of, in one object. The change move has
  // to be able to put the old rule back after a rejection, and doing that from a
  // hand-written list of fields is what silently dropped `na_rule` the first
  // time the rule grew a field. The list lives here, next to the declarations,
  // and nowhere else.
  // A whole rule, for the change move to put back when its proposal is
  // rejected. The level set belongs in here as much as the cutpoint does:
  // restoring `var` and `group` from a categorical rule while leaving the mask
  // cleared, or the other way round, leaves a node whose rule says it is one
  // kind and whose `var` indexes the other kind's matrix.
  struct Rule {
    int var;
    int group;
    double val;
    double lower;
    double upper;
    NaRule na_rule;
    std::vector<std::uint32_t> mask;
  };

  Rule rule() const {
    return Rule{var, group, val, lower, upper, na_rule, mask};
  }

  void set_rule(const Rule& r) {
    var = r.var;
    group = r.group;
    val = r.val;
    lower = r.lower;
    upper = r.upper;
    na_rule = r.na_rule;
    mask = r.mask;
  }

  // Draw a splitting rule from the prior: a variable, what to do with its
  // missing values if it has any, and a cutpoint.
  void draw_rule();
  void draw_categorical_rule(int levels);
  void available_levels(int group_id, int levels,
                        std::vector<std::uint32_t>& out) const;

  // Draw a splitting rule from the prior and build the two children, dividing
  // this node's support between them.
  void birth_leaves(std::vector<double>* w_left = nullptr,
                    std::vector<double>* w_right = nullptr);
  void delete_leaves();

  // Draw a fresh splitting rule from the prior for an existing branch, keeping
  // its children in place. Used by the change move.
  void resample_rule(std::vector<double>* w_left = nullptr,
                     std::vector<double>* w_right = nullptr);

  // Recompute the children's support from this node's, for use after a rule or
  // bandwidth change.
  //
  // When `w_left` and `w_right` are given they are filled with the weights the
  // rule assigns to each of *this* node's observations, in this node's own order
  // -- which is what a two-child target needs, and what a separate pass used to
  // recompute by evaluating every gate a second time.
  void split_support(std::vector<double>* w_left = nullptr,
                     std::vector<double>* w_right = nullptr);

  double gate(int i) const;
};

// The probability that an observation at a branch goes left. The hard rule is
// x <= val, and the soft rule is its logistic relaxation, which recovers the
// hard rule as the bandwidth goes to zero.
//
// A missing value has no position to compare, so the rule decides it outright.
// That keeps the two children's weights summing to the parent's, which is what
// the whole scheme rests on, and it makes a missing value's path through the
// tree a hard one even when the rules are soft -- correctly, since there is
// nothing about being absent to smooth over.
// A level set as a bitmask, 32 levels to a word.
inline int mask_words(int levels) { return (levels + 31) / 32; }

inline bool mask_test(const std::vector<std::uint32_t>& mask, int level) {
  return (mask[static_cast<std::size_t>(level) >> 5] >>
          (level & 31)) & 1u;
}

inline void mask_set(std::vector<std::uint32_t>& mask, int level) {
  mask[static_cast<std::size_t>(level) >> 5] |= 1u << (level & 31);
}

// Whether a categorical rule sends this level's observations left. A rule on a
// set of levels is always hard, even in a soft tree: a gate is a smooth function
// of the distance from a cutpoint and there is no distance between two levels of
// a factor. On a 0/1 indicator column the old distance-based gate did apply, and
// at the default bandwidth it left one level with a fractional membership weight
// for 81% of cutpoints, which is not something anyone asked for.
inline double left_prob_categorical(int level,
                                    const std::vector<std::uint32_t>& mask,
                                    int na_rule) {
  if (level < 0) {
    return na_rule == NA_RIGHT ? 0.0 : 1.0;
  }

  if (na_rule == NA_ONLY) {
    return 0.0;
  }

  return mask_test(mask, level) ? 1.0 : 0.0;
}

inline double left_prob(double x, double val, double bandwidth, bool soft,
                        int na_rule, int gate = GATE_LOGISTIC) {
  if (std::isnan(x)) {
    return na_rule == NA_RIGHT ? 0.0 : 1.0;
  }

  // The rule asks whether the value is missing, and this one is not.
  if (na_rule == NA_ONLY) {
    return 0.0;
  }

  if (!soft) {
    return x <= val ? 1.0 : 0.0;
  }

  if (gate != GATE_LOGISTIC) {
    double half = bandwidth;

    half *= gate == GATE_SMOOTHSTEP ? SMOOTHSTEP_HALF_WIDTH
                                    : SMOOTHERSTEP_HALF_WIDTH;

    double t = 0.5 + 0.5 * (val - x) / half;

    if (t <= 0.0) {
      return 0.0;
    }

    if (t >= 1.0) {
      return 1.0;
    }

    if (gate == GATE_SMOOTHSTEP) {
      return t * t * (3.0 - 2.0 * t);
    }

    return t * t * t * (10.0 + t * (6.0 * t - 15.0));
  }

  return expit((val - x) / bandwidth);
}

// Recompute every node's support below this one, for use after the bandwidth
// changes and every gate along every path moves at once.
void rebuild_support(Node* node);

// A flat snapshot of the supports below a node, so that a rejected bandwidth
// move can be undone by copying rather than by evaluating every gate in the tree
// a second time. Nodes are visited in pre-order, and the tree's shape does not
// change in between, so one offset per node is enough to find each again.
struct SupportStore {
  std::vector<int> idx;
  std::vector<double> wt;
  std::vector<std::size_t> offset;

  void clear() {
    idx.clear();
    wt.clear();
    offset.clear();
  }
};

void save_support(const Node* node, SupportStore& store);
void restore_support(Node* node, const SupportStore& store, std::size_t& pos);

std::vector<Node*> leaves(Node* node);
void leaves(Node* node, std::vector<Node*>& out);
std::vector<Node*> not_grand_branches(Node* node);
void not_grand_branches(Node* node, std::vector<Node*>& out);
int num_leaves(Node* node);

// The count alone, for the acceptance ratios that need only how many there
// would be. Building the list to ask for its length allocated a vector on every
// birth.
int num_not_grand_branches(Node* node);

Node* rand_node(const std::vector<Node*>& nodes);

// Accumulate this tree's contribution to the predictor for every observation.
void accumulate(Node* node, arma::rowvec& eta, double sign);

// Evaluate the tree at new data.
void predict_accumulate(Node* node, const arma::mat& X, arma::vec& out,
                        double weight, int i);
// A tree's prediction for one dataset used to live here as `predict_tree()` and
// `predict_accumulate()`. Both were dead -- nothing in the package or the tests
// called them -- and both read a rule as a threshold on a column of X, which a
// categorical rule is not, so they would have silently mishandled one for
// whoever called them next. Prediction goes through the encoded forest in
// `bartisan_predict()`.

void get_var_counts(Node* node, arma::uvec& counts);
void collect_leaf_params(Node* node, std::vector<double>& out);

} // namespace bartisan

#endif
