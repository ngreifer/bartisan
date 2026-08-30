#include <RcppArmadillo.h>
#include <memory>
#include <vector>
#include "family.h"
#include "hypers.h"
#include "mcmc.h"
#include "node.h"
#include "random.h"

// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace bartisan;

namespace {

// Trees are stored for later prediction as a preorder sequence of five-number
// records: whether the node is a leaf, its splitting variable, its cutpoint,
// what its rule does with a missing value, and its leaf value. A branch's two
// subtrees follow it immediately, so decoding is a single pass with a cursor.
const int RECORD_SIZE = 5;

// `leaf_shift` is subtracted from every leaf value as the tree is written out.
// A tree's membership weights sum to one for every observation, so subtracting
// the same amount from all of its leaves moves that tree's contribution to the
// predictor by exactly that amount -- which is how a whole forest is recentered
// by 1/num_trees of the shift per tree. See Family::report_shift().
void encode_tree(const Node* node, std::vector<double>& out,
                 double leaf_shift) {
  out.push_back(node->is_leaf ? 1.0 : 0.0);
  out.push_back(static_cast<double>(node->var));
  out.push_back(node->val);
  out.push_back(static_cast<double>(node->na_rule));
  out.push_back(node->mu - leaf_shift);
  if (!node->is_leaf) {
    encode_tree(node->left, out, leaf_shift);
    encode_tree(node->right, out, leaf_shift);
  }
}

// Both subtrees are always walked, even when a hard rule gives one of them zero
// weight, so that the cursor advances past the whole tree exactly once.
double eval_tree(const double* record, int& pos, const arma::mat& X, int i,
                 double weight, double bandwidth, bool soft, int gate) {
  bool leaf = record[pos] > 0.5;
  int var = static_cast<int>(record[pos + 1]);
  double val = record[pos + 2];
  int na_rule = static_cast<int>(record[pos + 3]);
  double mu = record[pos + 4];
  pos += RECORD_SIZE;

  if (leaf) {
    return weight * mu;
  }

  double g = left_prob(X(i, var), val, bandwidth, soft, na_rule, gate);
  double left = eval_tree(record, pos, X, i, weight * g, bandwidth, soft, gate);
  double right = eval_tree(record, pos, X, i, weight * (1.0 - g), bandwidth,
                           soft, gate);
  return left + right;
}

// The sampler passes forests around as vectors of raw pointers, so ownership
// has to be released by hand. Doing it in a destructor rather than at the end of
// the fit covers the paths that do not reach the end: the interrupt the sampler
// checks for every iteration, and the error a family whose log density is an R
// function can raise.
struct ForestGuard {
  std::vector<std::vector<Tree*>>& forests;

  ~ForestGuard() {
    for (std::size_t h = 0; h < forests.size(); h++) {
      for (std::size_t t = 0; t < forests[h].size(); t++) {
        delete forests[h][t];
      }
      forests[h].clear();
    }
  }
};

template <typename T>
List wrap_matrices(const std::vector<T>& x) {
  List out(x.size());
  for (std::size_t h = 0; h < x.size(); h++) {
    out[h] = wrap(x[h]);
  }
  return out;
}

} // namespace

//' Fit a generalized BART model
//'
//' The workhorse behind [bartisan()]. Not intended to be called directly: it
//' assumes the design matrix has already been mapped to the unit interval and
//' the response already coerced to the form the requested family expects.
//'
//' @param X design matrix with entries in `[0, 1]`, possibly with `NA`.
//' @param has_na indicator per column of `X` of whether it contains a missing
//'   value. A rule on a column with none is not given a missing-value branch, so
//'   complete data reproduces the sampler exactly as it was.
//' @param y response, coerced by the calling family.
//' @param weights prior weights.
//' @param offset an `H` by `N` matrix of fixed contributions to the additive
//'   predictors.
//' @param group_probs sparse matrix whose columns are predictor groups.
//' @param family_name,link,family_opts the family specification.
//' @param control a list of sampler and prior settings.
//' @return A list of posterior draws and the encoded forests.
//' @keywords internal
// [[Rcpp::export(.bartisan_fit)]]
List bartisan_fit(const arma::mat& X, const arma::uvec& has_na,
                 const arma::vec& y,
                 const arma::vec& weights, const arma::mat& offset,
                 const arma::sp_mat& group_probs, std::string family_name,
                 std::string link, List family_opts, List control,
                 List random_spec) {

  int n = static_cast<int>(X.n_rows);

  if (has_na.n_elem != X.n_cols) {
    stop("`has_na` must have one value per column of `X`.");
  }

  std::unique_ptr<Family> family(make_family(family_name, link, y, weights,
                                             family_opts));

  // The probit likelihood is the margin of a Gaussian one, and working with the
  // Gaussian makes the target quadratic in the predictor -- which is worth a
  // large constant factor. Only for Bernoulli data, where a single latent per
  // observation carries the whole likelihood.
  std::vector<std::string> augment;

  if (control.containsElementNamed("augment")) {
    augment = as<std::vector<std::string>>(control["augment"]);
  }

  if (!augment.empty()) {
    Family* rewritten = augmented_family(family_name, link, y, weights,
                                         family_opts, augment);

    if (rewritten != nullptr) {
      family.reset(rewritten);
    }
  }

  int H = family->H;

  if (static_cast<int>(offset.n_rows) != H ||
      static_cast<int>(offset.n_cols) != n) {
    stop("`offset` must be a %d by %d matrix.", H, n);
  }

  // One tree count per additive predictor. A single value is recycled by
  // `bartisan()` before it gets here, so this is always length H. The forests
  // are stored back to back rather than as a rectangle, so a per-forest offset
  // is what indexes them.
  std::vector<int> num_trees = as<std::vector<int>>(control["num_trees"]);

  if (static_cast<int>(num_trees.size()) != H) {
    stop("`control$num_trees` must have one value per additive predictor (%d).",
         H);
  }

  // The trailing forests may be nuisance parameters rather than predictors:
  // pinned at depth zero so that each is a single scalar, and reported in `aux`
  // rather than in `eta`. `n_report` is how many forests the caller asked for.
  int n_pinned = family->num_pinned();
  int n_report = H - n_pinned;

  std::vector<int> tree_offset(H + 1, 0);
  for (int h = 0; h < H; h++) {
    tree_offset[h + 1] = tree_offset[h] + num_trees[h];
  }
  int num_burn = as<int>(control["num_burn"]);
  int num_thin = as<int>(control["num_thin"]);
  int num_save = as<int>(control["num_save"]);
  int num_print = as<int>(control["num_print"]);
  bool verbose = as<bool>(control["verbose"]);
  bool soft = as<bool>(control["soft"]);
  int gate = as<int>(control["gate"]);
  bool update_sigma_mu = as<bool>(control["update_sigma_mu"]);

  // One leaf scale per additive predictor, since a location-scale model puts
  // its two predictors on quite different scales.
  arma::vec sigma_mu_target = as<arma::vec>(control["sigma_mu"]);

  if (static_cast<int>(sigma_mu_target.n_elem) != H) {
    stop("`control$sigma_mu` must have one value per additive predictor (%d).",
         H);
  }

  // Linero (2025), Remark 2: raising sigma_mu from near zero over the first
  // part of burn-in is described as essential, because starting it at a large
  // value can leave the sampler stuck away from the stationary distribution.
  double ramp_fraction = as<double>(control["sigma_mu_ramp"]);
  int num_ramp = static_cast<int>(std::floor(ramp_fraction * num_burn));

  std::vector<std::unique_ptr<Hypers>> hypers;
  std::vector<std::vector<Tree*>> forests(H);
  ForestGuard guard = {forests};

  for (int h = 0; h < H; h++) {
    // A pinned forest gets a branching probability of zero, which makes every
    // birth proposal impossible, and a fixed leaf scale, because one leaf cannot
    // identify a scale of its own -- left to draw it, `sigma_mu` wanders over an
    // order of magnitude.
    bool pinned = h >= n_report;
    hypers.push_back(std::unique_ptr<Hypers>(new Hypers(
      group_probs, sigma_mu_target(h),
      pinned ? 0.0 : as<double>(control["gamma"]),
      as<double>(control["beta"]), as<double>(control["alpha"]),
      as<double>(control["alpha_scale"]), as<double>(control["alpha_shape_1"]),
      as<double>(control["alpha_shape_2"]), pinned ? false : update_sigma_mu,
      as<bool>(control["update_s"]), as<bool>(control["update_alpha"]), soft,
      as<double>(control["bandwidth"]),
      as<bool>(control["update_bandwidth"]),
      as<int>(control["bandwidth_every"]), gate)));

    for (int t = 0; t < num_trees[h]; t++) {
      forests[h].push_back(new Tree(hypers[h].get(), &X, &has_na));
    }
  }

  arma::mat eta = offset;
  // Forcing the blocked evaluation path lets the test suite check that it
  // reproduces the per-observation path exactly.
  bool block_eval = control.containsElementNamed("block_eval") &&
    as<bool>(control["block_eval"]);

  Context ctx(family.get(), &eta, block_eval);

  // The closed forms a quadratic target allows can be switched off, so that the
  // test suite can check they reproduce the general path.
  if (control.containsElementNamed("exact_quadratic")) {
    ctx.exact_quadratic = as<bool>(control["exact_quadratic"]);
  }

  if (control.containsElementNamed("generic_accumulate")) {
    ctx.generic_accumulate = as<bool>(control["generic_accumulate"]);
  }

  int num_groups = static_cast<int>(group_probs.n_cols);
  std::vector<std::string> aux_names = family->aux_names();
  int num_aux = static_cast<int>(aux_names.size());

  std::vector<arma::mat> eta_out(n_report,
                                 arma::mat(num_save, n, arma::fill::zeros));
  std::vector<arma::umat> counts(n_report, arma::umat(num_save, num_groups,
                                                     arma::fill::zeros));
  arma::mat sigma_mu_out(num_save, n_report, arma::fill::zeros);
  arma::mat aux_out(num_save, std::max(num_aux, 1), arma::fill::zeros);
  arma::mat bandwidth_out(num_save, std::max(tree_offset[n_report], 1),
                          arma::fill::zeros);
  arma::vec loglik_out(num_save, arma::fill::zeros);

  std::vector<double> forest_flat;
  std::vector<int> tree_start;
  tree_start.push_back(0);

  // A family whose own state has a different size at every draw -- a Dirichlet
  // process mixture -- reports it the same way, flat with per-draw offsets.
  std::vector<double> mixture_flat;
  std::vector<int> mixture_start;
  mixture_start.push_back(0);

  // Group-level random intercepts, one set per additive predictor. The prior
  // scale starts at, and is given a half-Cauchy prior with, the same value the
  // leaf scale uses, so a group effect and a tree's leaf are shrunk on the same
  // footing.
  std::unique_ptr<RandomEffects> ranef =
    make_random_effects(random_spec, H, sigma_mu_target(0),
                        as<bool>(control["update_tau"]), &X);

  int num_ranef = 0;

  if (!ranef->empty()) {
    for (std::size_t r = 0; r < ranef->terms[0].size(); r++) {
      num_ranef += ranef->terms[0][r]->num_levels;
    }
  }

  std::vector<arma::mat> ranef_out(
    ranef->empty() ? 0 : H, arma::mat(num_save, std::max(num_ranef, 1),
                                      arma::fill::zeros));
  std::vector<arma::mat> tau_out(
    ranef->empty() ? 0 : H,
    arma::mat(num_save, ranef->empty() ? 1 : ranef->terms[0].size(),
              arma::fill::zeros));

  auto sweep = [&]() {
    for (int h = 0; h < H; h++) {
      ctx.h = h;
      ctx.form = family->target_form(h);
      ctx.quadratic = family->is_quadratic(h);
      family->before_forest(h, eta);
      update_forest(forests[h], ctx, *hypers[h]);
      // After the forest, so that the intercepts are drawn against the
      // predictor the trees have just settled on rather than the one before it.
      ranef->update(ctx, h);
    }
    family->update_aux(eta);
  };

  if (verbose) {
    Rcout << "Running " << num_burn << " warmup and " << num_save * num_thin
          << " sampling iterations\n";
  }

  for (int h = 0; h < H; h++) {
    hypers[h]->adapt = true;
  }

  for (int iter = 0; iter < num_burn; iter++) {
    if (iter < num_ramp) {
      double fraction = static_cast<double>(iter + 1) /
        static_cast<double>(num_ramp);
      for (int h = 0; h < H; h++) {
        hypers[h]->sigma_mu = sigma_mu_target(h) * fraction;
        hypers[h]->update_sigma_mu = false;
      }
      for (int h = n_report; h < H; h++) {
        // A pinned forest is not ramped: its scale is a prior the caller set,
        // not something the sampler is working its way towards.
        hypers[h]->sigma_mu = sigma_mu_target(h);
      }
    }

    sweep();
    Rcpp::checkUserInterrupt();

    if (verbose && num_print > 0 && (iter + 1) % num_print == 0) {
      Rcout << "\rWarmup " << iter + 1 << " / " << num_burn << "        ";
    }
  }
  if (verbose && num_burn > 0) {
    Rcout << "\n";
  }

  // Put the leaf scale at its target and hand it back to its own update, and
  // freeze the tuned proposals, before any draw is retained. Restoring here
  // rather than inside the loop matters when the ramp spans the whole of warmup,
  // where a restore conditioned on reaching iteration num_ramp never fires and
  // would leave the leaf scale pinned for the entire sampling phase.
  for (int h = 0; h < H; h++) {
    hypers[h]->sigma_mu = sigma_mu_target(h);
    hypers[h]->update_sigma_mu = h < n_report ? update_sigma_mu : false;
    hypers[h]->adapt = false;
  }

  for (int iter = 0; iter < num_save; iter++) {
    for (int j = 0; j < num_thin; j++) {
      sweep();
    }
    Rcpp::checkUserInterrupt();

    // The chart the draw is recorded in may differ from the one the sampler
    // works in; for every family but the ordinal ones this is zero.
    arma::vec shift = family->report_shift(eta);

    for (int h = 0; h < n_report; h++) {
      eta_out[h].row(iter) = eta.row(h) - shift(h);
      sigma_mu_out(iter, h) = hypers[h]->sigma_mu;

      arma::uvec var_counts = arma::zeros<arma::uvec>(num_groups);
      double per_tree = static_cast<double>(num_trees[h]);

      for (int t = 0; t < num_trees[h]; t++) {
        get_var_counts(forests[h][t]->root, var_counts);
        bandwidth_out(iter, tree_offset[h] + t) = forests[h][t]->bandwidth;
        encode_tree(forests[h][t]->root, forest_flat, shift(h) / per_tree);
        tree_start.push_back(static_cast<int>(forest_flat.size()));
      }
      counts[h].row(iter) = var_counts.t();
    }

    if (num_aux > 0) {
      aux_out.row(iter) = family->aux_values_shifted(shift).t();
    }

    {
      arma::vec drawn = family->mixture_flat();

      for (arma::uword k = 0; k < drawn.n_elem; k++) {
        mixture_flat.push_back(drawn(k));
      }

      mixture_start.push_back(static_cast<int>(mixture_flat.size()));
    }

    for (std::size_t h = 0; h < ranef_out.size(); h++) {
      int at = 0;

      for (std::size_t r = 0; r < ranef->terms[h].size(); r++) {
        arma::vec b = ranef->terms[h][r]->values();
        ranef_out[h](iter, arma::span(at, at + b.n_elem - 1)) = b.t();
        tau_out[h](iter, r) = ranef->terms[h][r]->tau;
        at += static_cast<int>(b.n_elem);
      }
    }
    loglik_out(iter) = family->reported_loglik(eta);

    if (verbose && num_print > 0 && (iter + 1) % num_print == 0) {
      Rcout << "\rSampling " << iter + 1 << " / " << num_save << "        ";
    }
  }
  if (verbose && num_save > 0) {
    Rcout << "\n";
  }

  List out;
  out["eta"] = wrap_matrices(eta_out);
  out["counts"] = wrap_matrices(counts);
  out["sigma_mu"] = sigma_mu_out;
  out["bandwidth"] = bandwidth_out;
  out["loglik"] = loglik_out;
  out["forest_flat"] = forest_flat;
  out["tree_start"] = tree_start;

  if (!mixture_flat.empty()) {
    out["mixture_flat"] = mixture_flat;
    out["mixture_start"] = mixture_start;
  }
  out["num_forest"] = n_report;
  out["num_trees"] = std::vector<int>(num_trees.begin(),
                                      num_trees.begin() + n_report);

  if (num_aux > 0) {
    out["aux"] = aux_out;
    out["aux_names"] = aux_names;
  }

  if (!ranef->empty()) {
    out["ranef"] = wrap_matrices(ranef_out);
    out["tau"] = wrap_matrices(tau_out);
  }

  return out;
}

//' Evaluate stored forests at new data
//'
//' @param X design matrix with entries in `[0, 1]`.
//' @param forest_flat,tree_start the encoded forests returned by
//'   `.bartisan_fit()`.
//' @param bandwidth a matrix of per-tree bandwidths.
//' @param num_forest,num_trees,num_save dimensions of the stored chain.
//' @param soft whether the decision rules are soft.
//' @param gate which gate the soft rules use; see `GateShape` in `node.h`.
//' @param iterations the zero-based saved iterations to evaluate.
//' @return A list of `num_forest` matrices of additive predictors.
//' @keywords internal
// [[Rcpp::export(.bartisan_predict)]]
List bartisan_predict(const arma::mat& X, const std::vector<double>& forest_flat,
                     const std::vector<int>& tree_start,
                     const arma::mat& bandwidth, int num_forest,
                     const std::vector<int>& num_trees,
                     int num_save, bool soft, int gate,
                     const std::vector<int>& iterations) {

  int n = static_cast<int>(X.n_rows);
  int num_iter = static_cast<int>(iterations.size());

  if (static_cast<int>(num_trees.size()) != num_forest) {
    stop("`num_trees` must have one value per forest.");
  }

  std::vector<int> tree_offset(num_forest + 1, 0);
  for (int h = 0; h < num_forest; h++) {
    tree_offset[h + 1] = tree_offset[h] + num_trees[h];
  }
  int total_trees = tree_offset[num_forest];

  std::vector<arma::mat> out(num_forest, arma::mat(num_iter, n,
                                                  arma::fill::zeros));

  for (int s = 0; s < num_iter; s++) {
    int iter = iterations[s];
    if (iter < 0 || iter >= num_save) {
      stop("`iterations` must index the saved draws.");
    }
    for (int h = 0; h < num_forest; h++) {
      for (int t = 0; t < num_trees[h]; t++) {
        // Trees were written iteration-major, then forest, then tree, with the
        // forests back to back rather than as a rectangle.
        int flat_index = iter * total_trees + tree_offset[h] + t;
        int begin = tree_start[flat_index];
        double band = bandwidth(iter, tree_offset[h] + t);
        for (int i = 0; i < n; i++) {
          int pos = begin;
          out[h](s, i) += eval_tree(forest_flat.data(), pos, X, i, 1.0, band,
                                    soft, gate);
        }
      }
    }
    Rcpp::checkUserInterrupt();
  }

  List result(num_forest);
  for (int h = 0; h < num_forest; h++) {
    result[h] = wrap(out[h]);
  }
  return result;
}

//' Conditional log density of the outcome at stored posterior draws
//'
//' Evaluates the family's log density for each observation at each posterior
//' draw of the additive predictors, using the nuisance parameters drawn at the
//' same iteration. This is the likelihood contribution of an observation, so it
//' is a density for a continuous response, a probability for a discrete one,
//' and a survival probability for a censored survival time.
//'
//' @param y the outcome, coerced as the family expects.
//' @param weights prior weights.
//' @param eta_draws a list of `H` matrices of draws by observations.
//' @param family_name,link,family_opts the family specification.
//' @param aux a matrix of draws by nuisance parameters, with zero columns when
//'   the family has none.
//' @return A matrix of draws by observations.
//' @keywords internal
// [[Rcpp::export(.bartisan_logdens)]]
arma::mat bartisan_logdens(const arma::vec& y, const arma::vec& weights,
                          const List& eta_draws, std::string family_name,
                          std::string link, List family_opts,
                          const arma::mat& aux) {

  std::unique_ptr<Family> family(make_family(family_name, link, y, weights,
                                             family_opts));

  int H = family->H;
  int n = static_cast<int>(y.n_elem);

  if (static_cast<int>(eta_draws.size()) != H) {
    stop("`eta_draws` must have one matrix per additive predictor (%d).", H);
  }

  std::vector<arma::mat> eta(H);
  for (int h = 0; h < H; h++) {
    eta[h] = as<arma::mat>(eta_draws[h]);
    if (static_cast<int>(eta[h].n_cols) != n) {
      stop("`eta_draws` and `y` disagree about the number of observations.");
    }
  }

  int num_draws = static_cast<int>(eta[0].n_rows);
  bool has_aux = aux.n_cols > 0;

  if (has_aux && static_cast<int>(aux.n_rows) != num_draws) {
    stop("`aux` must have one row per draw.");
  }

  arma::mat out(num_draws, n, arma::fill::zeros);
  arma::mat draw(H, n);
  arma::vec row(n);

  for (int d = 0; d < num_draws; d++) {
    if (has_aux) {
      family->set_aux(aux.row(d).t());
    }
    for (int i = 0; i < n; i++) {
      for (int h = 0; h < H; h++) {
        draw(h, i) = eta[h](d, i);
      }
    }
    family->logdens_all(draw, row.memptr());
    out.row(d) = row.t();
    Rcpp::checkUserInterrupt();
  }

  return out;
}

//' Category probabilities of a multinomial probit fit, by simulation
//'
//' The probability that the argmax of a correlated Gaussian vector falls in each
//' category is an orthant probability with no closed form, so it is simulated.
//' Fresh draws are taken on every call, which makes the estimate unbiased; the
//' Monte Carlo error is then averaged down by the posterior draws, so a modest
//' number of replicates per draw is enough for a posterior mean.
//'
//' @param eta_draws a list of one draws-by-observations matrix per latent
//'   variable.
//' @param sigma a matrix of draws by the lower triangle of the covariance
//'   matrix, column-major within a row, as `aux` stores it.
//' @param replicates simulation replicates per draw and observation.
//' @return An array of draws by observations by categories.
//' @keywords internal
// [[Rcpp::export(.bartisan_mnp_probs)]]
arma::cube bartisan_mnp_probs(const List& eta_draws, const arma::mat& sigma,
                             int replicates) {

  int C = static_cast<int>(eta_draws.size());

  if (C < 1) {
    stop("`eta_draws` must have at least one matrix.");
  }

  std::vector<arma::mat> eta(C);

  for (int l = 0; l < C; l++) {
    eta[l] = as<arma::mat>(eta_draws[l]);
  }

  int num_draws = static_cast<int>(eta[0].n_rows);
  int n = static_cast<int>(eta[0].n_cols);

  if (static_cast<int>(sigma.n_rows) != num_draws ||
      static_cast<int>(sigma.n_cols) != C * (C + 1) / 2) {
    stop("`sigma` must have one row per draw and one column per lower-triangle "
         "entry.");
  }

  arma::cube out(num_draws, n, C + 1, arma::fill::zeros);
  arma::mat covariance(C, C);
  arma::mat factor(C, C);
  arma::vec value(C);
  arma::vec noise(C);

  for (int d = 0; d < num_draws; d++) {
    int at = 0;

    for (int l = 0; l < C; l++) {
      for (int k = 0; k <= l; k++) {
        covariance(l, k) = sigma(d, at);
        covariance(k, l) = sigma(d, at);
        at++;
      }
    }

    if (!arma::chol(factor, covariance, "lower")) {
      factor = arma::eye<arma::mat>(C, C);
    }

    for (int i = 0; i < n; i++) {
      for (int r = 0; r < replicates; r++) {
        for (int l = 0; l < C; l++) {
          noise(l) = norm_rand();
        }

        int best = -1;
        double largest = 0.0;

        for (int l = 0; l < C; l++) {
          double shift = 0.0;

          for (int k = 0; k <= l; k++) {
            shift += factor(l, k) * noise(k);
          }

          value(l) = eta[l](d, i) + shift;

          if (value(l) >= largest) {
            largest = value(l);
            best = l;
          }
        }

        out(d, i, best + 1) += 1.0;
      }
    }

    Rcpp::checkUserInterrupt();
  }

  out /= static_cast<double>(replicates);

  return out;
}

//' Score and information of a family, analytic or by differences
//'
//' Exists so that the test suite can check each family's analytic derivatives
//' against central differences of its own log density.
//'
//' @param y,weights,eta_draws,family_name,link,family_opts,aux as for
//'   `.bartisan_logdens()`.
//' @param component which additive predictor to differentiate with respect to.
//' @param by_difference use central differences instead of the analytic form.
//' @param blocked evaluate a whole draw at once through the family's block
//'   methods rather than one observation at a time. The two paths should agree;
//'   they differ for a family whose per-observation route falls back on
//'   differences while its block route does not.
//' @return A list with matrices `d1` and `info`, draws by observations.
//' @keywords internal
// [[Rcpp::export(.bartisan_derivs)]]
List bartisan_derivs(const arma::vec& y, const arma::vec& weights,
                    const List& eta_draws, std::string family_name,
                    std::string link, List family_opts, const arma::mat& aux,
                    int component, bool by_difference, bool blocked = false) {

  std::unique_ptr<Family> family(make_family(family_name, link, y, weights,
                                             family_opts));

  int H = family->H;
  int n = static_cast<int>(y.n_elem);

  std::vector<arma::mat> eta(H);
  for (int h = 0; h < H; h++) {
    eta[h] = as<arma::mat>(eta_draws[h]);
  }

  int num_draws = static_cast<int>(eta[0].n_rows);
  bool has_aux = aux.n_cols > 0;

  arma::mat d1_out(num_draws, n, arma::fill::zeros);
  arma::mat info_out(num_draws, n, arma::fill::zeros);
  std::vector<double> column(H);

  if (blocked) {
    arma::mat draw(H, n);
    arma::vec a(n);
    arma::vec b(n);
    std::vector<int> idx(n);

    for (int i = 0; i < n; i++) {
      idx[i] = i;
    }

    for (int d = 0; d < num_draws; d++) {
      if (has_aux) {
        family->set_aux(aux.row(d).t());
      }
      for (int i = 0; i < n; i++) {
        for (int h = 0; h < H; h++) {
          draw(h, i) = eta[h](d, i);
        }
      }
      family->score_info_block(idx.data(), n, draw.memptr(), component,
                               a.memptr(), b.memptr());
      // The block methods apply the prior weight; strip it so that the result
      // is comparable with the per-observation route.
      d1_out.row(d) = (a / weights).t();
      info_out.row(d) = (b / weights).t();
    }

    return List::create(_["d1"] = d1_out, _["info"] = info_out);
  }

  for (int d = 0; d < num_draws; d++) {
    if (has_aux) {
      family->set_aux(aux.row(d).t());
    }
    for (int i = 0; i < n; i++) {
      for (int h = 0; h < H; h++) {
        column[h] = eta[h](d, i);
      }
      double a;
      double b;
      if (by_difference) {
        family->score_info_by_difference(i, column.data(), component, &a, &b);
      }
      else {
        family->score_info(i, column.data(), component, &a, &b);
        // score_info() applies the prior weight; the difference version does
        // not, so strip it here to make the two comparable.
        a /= weights(i);
        b /= weights(i);
      }
      d1_out(d, i) = a;
      info_out(d, i) = b;
    }
  }

  return List::create(_["d1"] = d1_out, _["info"] = info_out);
}
