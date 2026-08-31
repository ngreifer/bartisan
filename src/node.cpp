#include "node.h"

namespace bartisan {

Tree::Tree(Hypers* hypers_, const arma::mat* X_, const arma::uvec* has_na_,
           const arma::imat* codes_, const arma::ivec* cat_col_,
           const arma::ivec* n_levels_)
  : hypers(hypers_), X(X_), has_na(has_na_), codes(codes_), cat_col(cat_col_),
    n_levels(n_levels_) {
  bandwidth = hypers_->bandwidth_scale;
  log_step = std::log(5.0);
  attempts = 0;
  sweeps = 0;
  root = new Node(this, static_cast<int>(X_->n_rows));
}

Tree::~Tree() {
  delete root;

  for (std::size_t k = 0; k < pool.size(); k++) {
    delete pool[k];
  }
}

Node* Tree::take_node(Node* parent) {
  if (pool.empty()) {
    return new Node(parent);
  }

  Node* out = pool.back();
  pool.pop_back();
  out->init_as_child(parent);
  return out;
}

void Tree::give_node(Node* node) {
  if (node == nullptr) {
    return;
  }

  // Only a leaf can be recycled; anything else would leak its subtree.
  if (!node->is_leaf) {
    delete node;
    return;
  }

  // The vectors keep their capacity, which is half the point.
  node->idx.clear();
  node->wt.clear();
  pool.push_back(node);
}

Node::Node(Tree* tree_, int n_obs) {
  tree = tree_;
  left = nullptr;
  right = nullptr;
  parent = nullptr;
  is_root = true;
  is_leaf = true;
  depth = 0;
  var = 0;
  group = 0;
  val = 0.0;
  lower = 0.0;
  upper = 1.0;
  na_rule = NA_LEFT;
  mask.clear();
  mu = 0.0;
  mu_star = 0.0;
  v_star = 1.0;

  idx.reserve(n_obs);

  if (tree->hypers->soft) {
    wt.assign(n_obs, 1.0);
  }

  for (int i = 0; i < n_obs; i++) {
    idx.push_back(i);
  }
}

Node::Node(Node* parent_) {
  init_as_child(parent_);
}

void Node::init_as_child(Node* parent_) {
  tree = parent_->tree;
  left = nullptr;
  right = nullptr;
  parent = parent_;
  is_root = false;
  is_leaf = true;
  depth = parent_->depth + 1;
  var = 0;
  group = 0;
  val = 0.0;
  lower = 0.0;
  upper = 1.0;
  na_rule = NA_LEFT;
  mask.clear();
  mu = 0.0;
  mu_star = 0.0;
  v_star = 1.0;
}

Node::~Node() {
  delete left;
  delete right;
}

bool Node::is_left() const {
  if (parent == nullptr) {
    return true;
  }
  return this == parent->left;
}

// The admissible range for this node's cutpoint is bounded by the nearest
// ancestor that splits on the same variable. An ancestor that splits on whether
// the variable is missing says nothing about where its observed values lie, so
// it constrains nothing and is passed over.
void Node::get_limits() {
  lower = 0.0;
  upper = 1.0;
  const Node* y = this;
  bool keep_going = !y->is_root;
  while (keep_going) {
    bool from_left = y->is_left();
    y = y->parent;
    keep_going = !y->is_root;
    // A categorical ancestor constrains a set of levels, not an interval, and
    // its `var` is a column of level codes rather than of X, so it can collide
    // with a numeric column index. Both have to be numeric rules for the
    // ancestor to say anything about where this cutpoint may fall.
    if (y->var == var && !y->is_categorical() && y->na_rule != NA_ONLY) {
      keep_going = false;
      if (from_left) {
        upper = y->val;
        lower = y->lower;
      }
      else {
        upper = y->upper;
        lower = y->val;
      }
    }
  }
}

double Node::gate(int i) const {
  if (is_categorical()) {
    return left_prob_categorical((*tree->codes)(i, var), mask, na_rule);
  }

  return left_prob((*tree->X)(i, var), val, tree->bandwidth,
                   tree->hypers->soft, na_rule, tree->hypers->gate);
}

// The levels of this node's categorical group that can still reach it. Every
// ancestor splitting on the same group narrows it: going left keeps that
// ancestor's set, going right keeps its complement. Masks are absolute level
// sets, so intersecting over all such ancestors in any order gives the same
// answer as walking down from the root.
void Node::available_levels(int group_id, int levels,
                            std::vector<std::uint32_t>& out) const {
  int words = mask_words(levels);
  out.assign(words, 0u);

  for (int k = 0; k < levels; k++) {
    mask_set(out, k);
  }

  const Node* y = this;

  while (!y->is_root) {
    bool from_left = y->is_left();
    y = y->parent;

    if (y->group != group_id || !y->is_categorical() ||
        y->na_rule == NA_ONLY) {
      continue;
    }

    for (int w = 0; w < words; w++) {
      out[w] &= from_left ? y->mask[w] : ~y->mask[w];
    }
  }
}

// The variable comes from the sparsity weights and the cutpoint is uniform on
// the range the ancestors leave open, exactly as before. A variable with missing
// values additionally needs the rule to say what becomes of them, and that is
// drawn uniformly over the three options. A variable with no missing values is
// not given the extra draw at all, so the prior -- and the sequence of random
// numbers -- is unchanged for complete data.
//
// Drawing this from the prior is what keeps it out of every acceptance ratio:
// birth and change both propose a rule from its prior, so the prior and the
// proposal density carry the same factor and it cancels, just as the choice of
// variable and cutpoint already do.
void Node::draw_rule() {
  arma::uvec group_var = tree->hypers->sample_var();
  group = static_cast<int>(group_var(0));
  var = static_cast<int>(group_var(1));

  int levels = tree->levels_of(group);

  if (levels > 0) {
    draw_categorical_rule(levels);
    return;
  }

  mask.clear();
  na_rule = NA_LEFT;

  if (tree->splits_on_missing(var)) {
    na_rule = static_cast<NaRule>(sample_class(NUM_NA_RULES));
  }

  get_limits();
  val = lower + (upper - lower) * unif_rand();
}

// A subset of the levels still available at this node, drawn by sending each of
// them left with probability one half and rejecting the two draws that would
// leave a child empty. That rejection is part of the prior rather than a repair
// of it: the prior on a categorical rule *is* the uniform distribution over the
// 2^m - 2 non-degenerate subsets of the m available levels, and rejection
// sampling draws from exactly that, so the rule still cancels out of every
// acceptance ratio the way the variable and the cutpoint do.
//
// `var` becomes the group's column of level codes rather than a column of X,
// which is what `gate()` and the encoded tree both index by.
void Node::draw_categorical_rule(int levels) {
  var = tree->code_col(group);

  na_rule = NA_LEFT;

  if (tree->codes != nullptr && tree->cat_col != nullptr) {
    // A group with a missing value needs the rule to say what becomes of it,
    // drawn from its prior exactly as for a numeric column.
    bool any_missing = false;

    for (arma::uword i = 0; i < tree->codes->n_rows && !any_missing; i++) {
      any_missing = (*tree->codes)(i, var) < 0;
    }

    if (any_missing) {
      na_rule = static_cast<NaRule>(sample_class(NUM_NA_RULES));
    }
  }

  std::vector<std::uint32_t> avail;
  available_levels(group, levels, avail);

  std::vector<int> open;
  open.reserve(levels);

  for (int k = 0; k < levels; k++) {
    if (mask_test(avail, k)) {
      open.push_back(k);
    }
  }

  int words = mask_words(levels);

  // Fewer than two levels left: nothing to divide. The rule is recorded as a
  // categorical one that sends everything left, and the empty child makes the
  // move fail the same way an exhausted numeric column does.
  if (open.size() < 2u) {
    mask = avail;
    if (mask.empty()) {
      mask.assign(words, 0u);
    }
    return;
  }

  mask.assign(words, 0u);

  while (true) {
    int taken = 0;

    for (int w = 0; w < words; w++) {
      mask[w] = 0u;
    }

    for (std::size_t j = 0; j < open.size(); j++) {
      if (unif_rand() < 0.5) {
        mask_set(mask, open[j]);
        taken++;
      }
    }

    if (taken > 0 && taken < static_cast<int>(open.size())) {
      return;
    }
  }
}

void Node::birth_leaves(std::vector<double>* w_left,
                        std::vector<double>* w_right) {
  if (!is_leaf) {
    return;
  }

  is_leaf = false;
  draw_rule();

  left = tree->take_node(this);
  right = tree->take_node(this);
  split_support(w_left, w_right);
}

void Node::resample_rule(std::vector<double>* w_left,
                         std::vector<double>* w_right) {
  if (is_leaf) {
    return;
  }
  draw_rule();
  split_support(w_left, w_right);
}

void save_support(const Node* node, SupportStore& store) {
  if (node->is_leaf) {
    return;
  }

  const Node* kids[2] = {node->left, node->right};

  for (int c = 0; c < 2; c++) {
    const Node* kid = kids[c];
    store.offset.push_back(store.idx.size());
    store.idx.insert(store.idx.end(), kid->idx.begin(), kid->idx.end());
    store.wt.insert(store.wt.end(), kid->wt.begin(), kid->wt.end());
  }

  save_support(node->left, store);
  save_support(node->right, store);
}

void restore_support(Node* node, const SupportStore& store, std::size_t& pos) {
  if (node->is_leaf) {
    return;
  }

  Node* kids[2] = {node->left, node->right};

  for (int c = 0; c < 2; c++) {
    Node* kid = kids[c];
    std::size_t from = store.offset[pos];
    std::size_t to = pos + 1 < store.offset.size() ? store.offset[pos + 1]
                                                   : store.idx.size();
    pos++;

    kid->idx.assign(store.idx.begin() + from, store.idx.begin() + to);

    // A hard tree stores no weights at all, so there is nothing to put back.
    if (store.wt.empty()) {
      kid->wt.clear();
      continue;
    }

    kid->wt.assign(store.wt.begin() + from, store.wt.begin() + to);
  }

  restore_support(node->left, store, pos);
  restore_support(node->right, store, pos);
}

void rebuild_support(Node* node) {
  if (node->is_leaf) {
    return;
  }
  node->split_support();
  rebuild_support(node->left);
  rebuild_support(node->right);
}

void Node::delete_leaves() {
  tree->give_node(left);
  tree->give_node(right);
  left = nullptr;
  right = nullptr;
  is_leaf = true;
}

void Node::split_support(std::vector<double>* w_left,
                         std::vector<double>* w_right) {
  // The children's vectors are sized to the parent's support -- the most either
  // could take -- filled by index, and trimmed. `push_back` tested capacity on
  // every element of every one of the four vectors, in the innermost loop of the
  // whole sampler; the value-initialization that `resize` does to the unused tail
  // is a memset and turns out to cost a fraction of that. Measured at 20% of a
  // soft-rule fit, which is more than the two structural fusions in this round
  // put together.
  bool soft = tree->hypers->soft;
  std::size_t n = idx.size();
  bool record = w_left != nullptr;

  if (record) {
    w_left->resize(n);
    w_right->resize(n);
  }

  left->idx.resize(n);
  right->idx.resize(n);
  std::size_t nl = 0;
  std::size_t nr = 0;

  if (!soft) {
    left->wt.clear();
    right->wt.clear();

    for (std::size_t k = 0; k < n; k++) {
      int i = idx[k];
      double g = gate(i);

      if (record) {
        (*w_left)[k] = g;
        (*w_right)[k] = 1.0 - g;
      }

      if (g > 0.5) {
        left->idx[nl++] = i;
      }
      else {
        right->idx[nr++] = i;
      }
    }

    left->idx.resize(nl);
    right->idx.resize(nr);
    return;
  }

  left->wt.resize(n);
  right->wt.resize(n);

  for (std::size_t k = 0; k < n; k++) {
    int i = idx[k];
    double g = gate(i);
    double wl = wt[k] * g;
    double wr = wt[k] - wl;

    if (record) {
      (*w_left)[k] = wl;
      (*w_right)[k] = wr;
    }

    if (wl > WEIGHT_TOL) {
      left->idx[nl] = i;
      left->wt[nl] = wl;
      nl++;
    }
    if (wr > WEIGHT_TOL) {
      right->idx[nr] = i;
      right->wt[nr] = wr;
      nr++;
    }
  }

  left->idx.resize(nl);
  left->wt.resize(nl);
  right->idx.resize(nr);
  right->wt.resize(nr);
}

void leaves(Node* node, std::vector<Node*>& out) {
  if (node->is_leaf) {
    out.push_back(node);
  }
  else {
    leaves(node->left, out);
    leaves(node->right, out);
  }
}

std::vector<Node*> leaves(Node* node) {
  std::vector<Node*> out;
  leaves(node, out);
  return out;
}

int num_leaves(Node* node) {
  if (node->is_leaf) {
    return 1;
  }
  return num_leaves(node->left) + num_leaves(node->right);
}

void not_grand_branches(Node* node, std::vector<Node*>& out) {
  if (node->is_leaf) {
    return;
  }
  if (node->left->is_leaf && node->right->is_leaf) {
    out.push_back(node);
  }
  else {
    not_grand_branches(node->left, out);
    not_grand_branches(node->right, out);
  }
}

int num_not_grand_branches(Node* node) {
  if (node->is_leaf) {
    return 0;
  }
  if (node->left->is_leaf && node->right->is_leaf) {
    return 1;
  }
  return num_not_grand_branches(node->left) +
    num_not_grand_branches(node->right);
}

std::vector<Node*> not_grand_branches(Node* node) {
  std::vector<Node*> out;
  not_grand_branches(node, out);
  return out;
}

Node* rand_node(const std::vector<Node*>& nodes) {
  return nodes[sample_class(static_cast<int>(nodes.size()))];
}

void accumulate(Node* node, arma::rowvec& eta, double sign) {
  if (node->is_leaf) {
    std::size_t n = node->idx.size();
    const double* wt = node->weights();
    double contribution = sign * node->mu;

    if (wt == nullptr) {
      for (std::size_t k = 0; k < n; k++) {
        eta(node->idx[k]) += contribution;
      }
      return;
    }

    // The product is formed in the order the unweighted version above formed it
    // in, so that a soft tree's arithmetic is untouched by this.
    for (std::size_t k = 0; k < n; k++) {
      eta(node->idx[k]) += sign * wt[k] * node->mu;
    }
    return;
  }
  accumulate(node->left, eta, sign);
  accumulate(node->right, eta, sign);
}

void get_var_counts(Node* node, arma::uvec& counts) {
  if (node->is_leaf) {
    return;
  }
  counts(node->group) += 1;
  get_var_counts(node->left, counts);
  get_var_counts(node->right, counts);
}

void collect_leaf_params(Node* node, std::vector<double>& out) {
  if (node->is_leaf) {
    out.push_back(node->mu);
    return;
  }
  collect_leaf_params(node->left, out);
  collect_leaf_params(node->right, out);
}

} // namespace bartisan
