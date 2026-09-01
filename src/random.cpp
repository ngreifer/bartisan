#include "random.h"

namespace bartisan {

namespace {

// Hyperparameters that exist only so a level node has somewhere to read
// `soft = false` from. Nothing else on them is used: a level node is never
// split, so no branching probability, sparsity weight or bandwidth is ever
// consulted.
std::unique_ptr<Hypers> level_hypers(double tau) {
  arma::sp_mat one(1, 1);
  one(0, 0) = 1.0;

  return std::unique_ptr<Hypers>(new Hypers(one, tau, 0.95, 2.0, 1.0, 1.0, 0.5,
                                            1.0, false, false, false,
                                            /* soft = */ false, 0.1, false, 1,
                                            GATE_LOGISTIC));
}

} // namespace

RandomTerm::RandomTerm(const std::string& label_, const arma::ivec& level_of,
                       int num_levels_, double tau_start, double tau_scale_,
                       bool update_tau_, const arma::mat* X)
  : label(label_), num_levels(num_levels_), tau(tau_start),
    tau_scale(tau_scale_), update_tau(update_tau_) {

  hypers = level_hypers(tau_start);
  tree.reset(new Tree(hypers.get(), X, nullptr));

  // How many observations each level has, so the index vectors can be sized
  // once rather than grown.
  std::vector<int> counts(num_levels, 0);

  for (arma::uword i = 0; i < level_of.n_elem; i++) {
    int l = level_of(i);
    if (l >= 0 && l < num_levels) {
      counts[l]++;
    }
  }

  levels.reserve(num_levels);

  for (int l = 0; l < num_levels; l++) {
    Node* node = new Node(tree.get(), 0);
    node->idx.reserve(counts[l]);
    // A level's weights are all one, and an empty weight vector is how that is
    // said; see Node::weights().
    node->wt.clear();
    node->mu = 0.0;
    levels.push_back(node);
  }

  for (arma::uword i = 0; i < level_of.n_elem; i++) {
    int l = level_of(i);
    if (l >= 0 && l < num_levels) {
      levels[l]->idx.push_back(static_cast<int>(i));
    }
  }
}

RandomTerm::~RandomTerm() {
  for (std::size_t l = 0; l < levels.size(); l++) {
    delete levels[l];
  }
}

arma::vec RandomTerm::values() const {
  arma::vec out(num_levels);

  for (int l = 0; l < num_levels; l++) {
    out(l) = levels[l]->mu;
  }

  return out;
}

void RandomTerm::update(Context& ctx) {
  // The prior scale the level updates see is this term's, not the forest's.
  double forest_scale = ctx.sigma_mu;
  ctx.sigma_mu = tau;

  for (int l = 0; l < num_levels; l++) {
    // A level with no observations has only its prior, which is a direct draw.
    if (levels[l]->idx.empty()) {
      levels[l]->mu = norm_rand() * tau;
      continue;
    }

    update_scalar(levels[l], ctx);
  }

  ctx.sigma_mu = forest_scale;

  if (!update_tau || num_levels < 2) {
    return;
  }

  // The same half-Cauchy update the leaf scale uses, with the intercepts in
  // place of the leaf values -- which is what they are: a set of draws from a
  // common mean-zero normal whose scale is wanted.
  arma::vec b = values();
  double prec_old = std::pow(tau, -2.0);
  double prec_new = half_cauchy_update_precision_mh(b, prec_old, tau_scale);
  tau = std::pow(prec_new, -0.5);
  hypers->sigma_mu = tau;
}

void RandomEffects::accumulate_into(arma::mat& eta) const {
  for (std::size_t h = 0; h < terms.size(); h++) {
    for (std::size_t r = 0; r < terms[h].size(); r++) {
      const RandomTerm& term = *terms[h][r];

      for (int l = 0; l < term.num_levels; l++) {
        const Node* node = term.levels[l];
        double value = node->mu;

        for (std::size_t k = 0; k < node->idx.size(); k++) {
          eta(h, node->idx[k]) += value;
        }
      }
    }
  }
}

void RandomEffects::update(Context& ctx, int h) {
  if (h < 0 || static_cast<std::size_t>(h) >= terms.size()) {
    return;
  }

  for (std::size_t r = 0; r < terms[h].size(); r++) {
    terms[h][r]->update(ctx);
  }
}

std::unique_ptr<RandomEffects> make_random_effects(const Rcpp::List& spec, int H,
                                                   double tau_scale,
                                                   bool update_tau,
                                                   const arma::mat* X,
                                                   const arma::ivec& vc_column) {
  std::unique_ptr<RandomEffects> out(new RandomEffects());

  if (spec.size() == 0) {
    return out;
  }

  out->terms.resize(H);

  for (int h = 0; h < H; h++) {
    // A group intercept belongs to the control function, not to a coefficient.
    // Giving one to a coefficient's forest would make the coefficient itself
    // vary by group, which is a random slope -- and `split_random()` refuses
    // `(x | g)` in as many words, so producing one here without being asked
    // would contradict that refusal silently.
    if (static_cast<int>(vc_column.n_elem) > h && vc_column(h) >= 0) {
      continue;
    }

    for (int r = 0; r < spec.size(); r++) {
      Rcpp::List one = spec[r];
      std::string label = Rcpp::as<std::string>(one["label"]);
      arma::ivec level_of = Rcpp::as<arma::ivec>(one["levels"]);
      int num_levels = Rcpp::as<int>(one["num_levels"]);

      out->terms[h].push_back(std::unique_ptr<RandomTerm>(
        new RandomTerm(label, level_of, num_levels, tau_scale, tau_scale,
                       update_tau, X)));
    }
  }

  return out;
}

} // namespace bartisan
