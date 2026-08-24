#include "node.h"

namespace genbart {

Tree::Tree(Hypers* hypers_, const arma::mat* X_, const arma::uvec* has_na_)
  : hypers(hypers_), X(X_), has_na(has_na_) {
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
    if (y->var == var && y->na_rule != NA_ONLY) {
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
  return left_prob((*tree->X)(i, var), val, tree->bandwidth,
                   tree->hypers->soft, na_rule, tree->hypers->gate);
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

  na_rule = NA_LEFT;

  if (tree->splits_on_missing(var)) {
    na_rule = static_cast<NaRule>(sample_class(NUM_NA_RULES));
  }

  get_limits();
  val = lower + (upper - lower) * unif_rand();
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

void predict_accumulate(Node* node, const arma::mat& X, arma::vec& out,
                        double weight, int i) {
  if (weight <= Node::WEIGHT_TOL) {
    return;
  }
  if (node->is_leaf) {
    out(i) += weight * node->mu;
    return;
  }
  double g = left_prob(X(i, node->var), node->val, node->tree->bandwidth,
                       node->tree->hypers->soft, node->na_rule,
                       node->tree->hypers->gate);
  predict_accumulate(node->left, X, out, weight * g, i);
  predict_accumulate(node->right, X, out, weight * (1.0 - g), i);
}

arma::vec predict_tree(Tree* tree, const arma::mat& X) {
  arma::vec out = arma::zeros<arma::vec>(X.n_rows);
  int n = static_cast<int>(X.n_rows);
  for (int i = 0; i < n; i++) {
    predict_accumulate(tree->root, X, out, 1.0, i);
  }
  return out;
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

} // namespace genbart
