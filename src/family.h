#ifndef GENBART_FAMILY_H
#define GENBART_FAMILY_H

#include <RcppArmadillo.h>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include "utils.h"

namespace genbart {

// The shape of a family's log density as a function of one additive predictor,
// with the others held fixed. It matters because a leaf value enters the
// predictor linearly, so the log target over a leaf inherits the same shape --
// and where the shape is known, one pass over the node determines the whole
// function and everything after that is arithmetic.
enum TargetForm {
  // Nothing assumed. Every value of a leaf parameter costs its own pass.
  TARGET_GENERAL = 0,

  // c + a * eta - b * eta^2 / 2. The Gaussian response, and every likelihood
  // that genbart_control(augment = ) rewrites into a Gaussian one. Here the
  // Laplace approximation is exact: Fisher scoring lands on the mode in one
  // step from anywhere and the fitted normal is the conditional posterior.
  TARGET_QUADRATIC = 1,

  // c + a * eta + b * exp(eta): the Poisson. This is the shape Hill et al.
  // (2020, sec. 3.1.5) identify as common to the count models. It is not
  // quadratic, so the Laplace fit still has to be iterated -- but three numbers
  // from one pass determine the target everywhere, so the iteration runs on
  // scalars and costs no further passes.
  TARGET_EXP_UP = 2,

  // c + a * eta + b * exp(-eta): gamma with a log link, and the negative
  // binomial once its Poisson-gamma latent variable is introduced.
  TARGET_EXP_DOWN = 3
};

// The entire likelihood-specific surface of the sampler. Given the additive
// predictors eta, a family supplies the log density of each observation and its
// first two derivatives with respect to each predictor. Nothing else in the
// package knows which family is in use.
//
// eta is stored as an H x N matrix so that eta.colptr(i) is a contiguous block
// of H doubles holding the predictors for observation i.
//
// Derived classes implement the *_unit methods, which are the contribution of a
// single unit of prior weight. The base class applies the weights, so weights
// need only be correct in one place. For a binomial family this means y holds
// proportions and the weights hold the number of trials, matching glm().
struct Family {

  int H;
  int N;
  arma::vec y;
  arma::vec w;

  Family(const arma::vec& y_, const arma::vec& w_, int H_)
    : H(H_), N(static_cast<int>(y_.n_elem)), y(y_), w(w_),
      eta_free_total(0.0), norm_const_total(0.0) {
    scratch.resize(H_);
  }

  virtual ~Family() {}

  double logdens(int i, const double* eta) const {
    return w(i) * logdens_unit(i, eta);
  }

  double dlogdens(int i, const double* eta, int h) const {
    return w(i) * dlogdens_unit(i, eta, h);
  }

  // Minus the second derivative. Clamped at zero rather than at a positive
  // floor: a saturated observation may legitimately carry no information, and
  // the prior term added at the node keeps the total curvature positive.
  double info(int i, const double* eta, int h) const {
    return clamp_info(w(i) * info_unit(i, eta, h));
  }

  // Score and information together. Fisher scoring always wants both at the
  // same point, and computing them in one call lets the families that difference
  // the log density numerically share the evaluations: three instead of five.
  void score_info(int i, const double* eta, int h, double* d1,
                  double* d2) const {
    double a;
    double b;
    score_info_unit(i, eta, h, &a, &b);
    *d1 = w(i) * a;
    *d2 = clamp_info(w(i) * b);
  }

  virtual double logdens_unit(int i, const double* eta) const = 0;

  // The part of the log density that does not scale with the prior weight, such
  // as the binomial coefficient. It is constant in eta, so it cancels from
  // every acceptance ratio and the sampler ignores it; it is needed only when
  // reporting a properly normalized density to the user.
  virtual double log_norm_const(int i) const { return 0.0; }

  double logdens_full(int i, const double* eta) const {
    double out = logdens(i, eta) + log_norm_const(i);
    if (!eta_free.is_empty()) {
      out += w(i) * eta_free(i);
    }
    return out;
  }

  // Exposed so that a family that wraps another can borrow its eta-free part
  // rather than recomputing it.
  const arma::vec& eta_free_part() const { return eta_free; }

  // Everything in the log likelihood that does not depend on the predictors:
  // the eta-free terms and the normalizing constants, already summed over the
  // sample. A wrapping family delegates this to what it wraps.
  virtual double logdens_extra_total() const {
    return eta_free_total + norm_const_total;
  }

  // Terms of the log density that do not involve eta at all. They cancel from
  // every acceptance ratio, so `logdens_unit()` is free to omit them and the
  // sampler never pays for them; a reported density and the reported log
  // likelihood add them back. Recomputed whenever a nuisance parameter that
  // they depend on changes.
  void refresh_eta_free() {
    eta_free = compute_eta_free();
    eta_free_total = eta_free.is_empty() ? 0.0 : arma::dot(w, eta_free);

    // The normalizing constants depend only on the response and the weights, so
    // they are summed once here and carried into the reported log likelihood.
    norm_const_total = 0.0;
    for (int i = 0; i < N; i++) {
      norm_const_total += log_norm_const(i);
    }
  }

  virtual arma::vec compute_eta_free() const { return arma::vec(); }

  // Central difference, for families whose analytic derivative is awkward.
  virtual double dlogdens_unit(int i, const double* eta, int h) const {
    double d1;
    double d2;
    score_info_numeric(i, eta, h, &d1, &d2);
    return d1;
  }

  virtual double info_unit(int i, const double* eta, int h) const {
    double d1;
    double d2;
    score_info_numeric(i, eta, h, &d1, &d2);
    return d2;
  }

  // Defaults to the analytic pair. Families that fall back on differences
  // override this with score_info_numeric(), which produces both from the same
  // three log-density evaluations.
  virtual void score_info_unit(int i, const double* eta, int h, double* d1,
                               double* d2) const {
    *d1 = dlogdens_unit(i, eta, h);
    *d2 = info_unit(i, eta, h);
  }

  // Draw the nuisance parameters given the current predictors.
  virtual void update_aux(const arma::mat& eta) {}

  // Called before each forest's sweep, with the forest's index. An augmented
  // family redraws its latent variables here rather than in update_aux(),
  // because the augmentation has to condition on the predictors as they stand
  // when that forest is about to move -- which for a multinomial model means
  // after the other categories have already moved this sweep.
  virtual void before_forest(int h, const arma::mat& eta) {}

  virtual std::vector<std::string> aux_names() const {
    return std::vector<std::string>();
  }

  virtual arma::vec aux_values() const { return arma::vec(); }

  virtual void set_aux(const arma::vec& values) {}

  // Exposed only so that the test suite can check each family's analytic
  // derivatives against the central-difference versions.
  void score_info_by_difference(int i, const double* eta, int h, double* d1,
                                double* d2) const {
    score_info_numeric(i, eta, h, d1, d2);
  }

  // The shape of this family's log density in the hth additive predictor. The
  // sampler uses it to replace passes over the data with arithmetic; see
  // TargetForm above.
  virtual TargetForm target_form(int h) const { return TARGET_GENERAL; }

  bool is_quadratic(int h) const {
    return target_form(h) == TARGET_QUADRATIC;
  }

  // The sign in exp(sign * eta) for the two exponential forms, and zero for
  // the others.
  double exp_sign(int h) const {
    TargetForm form = target_form(h);

    if (form == TARGET_EXP_UP) {
      return 1.0;
    }

    if (form == TARGET_EXP_DOWN) {
      return -1.0;
    }

    return 0.0;
  }

  // Whether this family would rather be handed a whole leaf at a time. A family
  // whose cost is dominated by a fixed charge per call -- one whose log density
  // is an R function -- says yes; a compiled family says no, because
  // materializing the block costs more than the density it saves.
  virtual bool wants_block() const { return false; }

  // Batch evaluation. Every likelihood the sampler asks for is a sum over the
  // observations reaching one leaf, so a family that wants to amortize a fixed
  // cost can take the whole leaf at once. `block` holds the predictors in the
  // same H-by-n column-major layout as `eta`, so block + k * H are the
  // predictors of observation idx[k].
  //
  // The default implementations are the same loop the sampler runs inline, in
  // the same order, so the two paths agree to the last bit for any family that
  // does not override them.
  virtual void logdens_block(const int* idx, int n, const double* block,
                             double* out) const {
    for (int k = 0; k < n; k++) {
      out[k] = logdens(idx[k], block + k * H);
    }
  }

  virtual void score_info_block(const int* idx, int n, const double* block,
                                int h, double* d1, double* d2) const {
    for (int k = 0; k < n; k++) {
      score_info(idx[k], block + k * H, h, d1 + k, d2 + k);
    }
  }

  // Accumulate, over the observations a node supports, the three sums the
  // leaf-level work needs: the log density, the weighted score, and the weighted
  // information. `block` holds the predictors in the same H-by-n column-major
  // layout as `eta`, already carrying whatever leaf value is being evaluated.
  //
  // This exists to be overridden statically. The default below is the loop the
  // sampler used to run inline, and it costs four virtual calls per observation
  // that the compiler cannot inline; Concrete<Derived> replaces it with the same
  // loop against a concrete family type, which is one virtual call per leaf and
  // the family's arithmetic inlined. Once the leaf-level work has collapsed to a
  // single pass, that dispatch is most of what is left.
  virtual void accumulate1(const int* idx, const double* wt, int n,
                           const double* block, int h, double* f, double* s,
                           double* j) const {
    accumulate1_generic(idx, wt, n, block, h, f, s, j);
  }

  virtual void accumulate2(const int* idx, const double* w_left,
                           const double* w_right, int n, const double* block,
                           int h, double* f, double* g, double* info) const {
    accumulate2_generic(idx, w_left, w_right, n, block, h, f, g, info);
  }

  // The same three sums, taking each observation's predictors from `eta` and
  // shifting component h by `wt * shift`, rather than from a block the caller has
  // written them into first. The leaf's own contribution is already in `eta`, so
  // that shift is exactly what a trial value of a leaf parameter does -- which
  // makes this the hot path, and materializing the block for it was a full
  // write-then-read pass over the node that only a family reaching back into R
  // needs.
  //
  // The default is the same loop as accumulate1_generic(), reading from `eta`;
  // Concrete<Derived> replaces it with a statically dispatched one.
  virtual void accumulate1_at(const int* idx, const double* wt, int n,
                              const arma::mat& eta, int h, double shift,
                              double* f, double* s, double* j) const {
    accumulate1_at_generic(idx, wt, n, eta, h, shift, f, s, j);
  }

  virtual void accumulate1_from(const int* idx, const double* wt, int n,
                                const arma::mat& eta, const double* base,
                                double mu, int h, double* f, double* s,
                                double* j) const {
    accumulate1_from_generic(idx, wt, n, eta, base, mu, h, f, s, j);
  }

  virtual void accumulate2_from(const int* idx, const double* w_left,
                                const double* w_right, int n,
                                const arma::mat& eta, const double* base,
                                double mu_left, double mu_right, int h,
                                double* f, double* g, double* info) const {
    accumulate2_from_generic(idx, w_left, w_right, n, eta, base, mu_left,
                             mu_right, h, f, g, info);
  }

  void accumulate1_from_generic(const int* idx, const double* wt, int n,
                                const arma::mat& eta, const double* base,
                                double mu, int h, double* f, double* s,
                                double* j) const {
    double fa = 0.0;
    double sa = 0.0;
    double ja = 0.0;
    const double* eta_ptr = eta.memptr();

    for (int k = 0; k < n; k++) {
      int i = idx[k];
      double wk = wt == nullptr ? 1.0 : wt[k];
      const double* col = eta_ptr + static_cast<std::size_t>(i) * H;
      std::copy(col, col + H, scratch.begin());
      scratch[h] = base[k] + wk * mu;
      const double* e = scratch.data();
      double a;
      double b;
      score_info(i, e, h, &a, &b);
      fa += logdens(i, e);
      sa += wk * a;
      ja += wk * wk * b;
    }

    *f = fa;
    *s = sa;
    *j = ja;
  }

  void accumulate2_from_generic(const int* idx, const double* w_left,
                                const double* w_right, int n,
                                const arma::mat& eta, const double* base,
                                double mu_left, double mu_right, int h,
                                double* f, double* g, double* info) const {
    double fa = 0.0;
    const double* eta_ptr = eta.memptr();

    for (int k = 0; k < n; k++) {
      int i = idx[k];
      double wl = w_left[k];
      double wr = w_right[k];
      const double* col = eta_ptr + static_cast<std::size_t>(i) * H;
      std::copy(col, col + H, scratch.begin());
      scratch[h] = base[k] + wl * mu_left + wr * mu_right;
      const double* e = scratch.data();
      double a;
      double b;
      score_info(i, e, h, &a, &b);
      fa += logdens(i, e);
      g[0] += wl * a;
      g[1] += wr * a;
      info[0] += wl * wl * b;
      info[1] += wl * wr * b;
      info[2] += wr * wr * b;
    }

    *f = fa;
  }

  void accumulate1_at_generic(const int* idx, const double* wt, int n,
                              const arma::mat& eta, int h, double shift,
                              double* f, double* s, double* j) const {
    double fa = 0.0;
    double sa = 0.0;
    double ja = 0.0;
    const double* eta_ptr = eta.memptr();

    for (int k = 0; k < n; k++) {
      int i = idx[k];
      double wk = wt == nullptr ? 1.0 : wt[k];
      const double* col = eta_ptr + static_cast<std::size_t>(i) * H;
      std::copy(col, col + H, scratch.begin());
      scratch[h] = col[h] + wk * shift;
      const double* e = scratch.data();
      double a;
      double b;
      score_info(i, e, h, &a, &b);
      fa += logdens(i, e);
      sa += wk * a;
      ja += wk * wk * b;
    }

    *f = fa;
    *s = sa;
    *j = ja;
  }

  // The dynamically dispatched versions, kept reachable so that the test suite
  // can check a family's static override against them.
  //
  // `wt` may be null, which means every membership weight is one -- a hard rule
  // sends an observation entirely one way, so a hard tree stores none of them.
  // Templated on which case it is so the multiplication disappears at compile
  // time rather than being a branch per observation; because multiplying by 1.0
  // is exact, the two instantiations agree to the last bit.
  template <bool Weighted>
  void accumulate1_loop(const int* idx, const double* wt, int n,
                        const double* block, int h, double* f, double* s,
                        double* j) const {
    double fa = 0.0;
    double sa = 0.0;
    double ja = 0.0;

    for (int k = 0; k < n; k++) {
      int i = idx[k];
      const double* e = block + static_cast<std::size_t>(k) * H;
      double a;
      double b;
      score_info(i, e, h, &a, &b);
      fa += logdens(i, e);
      double wk = Weighted ? wt[k] : 1.0;
      sa += wk * a;
      ja += wk * wk * b;
    }

    *f = fa;
    *s = sa;
    *j = ja;
  }

  void accumulate1_generic(const int* idx, const double* wt, int n,
                           const double* block, int h, double* f, double* s,
                           double* j) const {
    if (wt == nullptr) {
      accumulate1_loop<false>(idx, wt, n, block, h, f, s, j);
      return;
    }

    accumulate1_loop<true>(idx, wt, n, block, h, f, s, j);
  }

  void accumulate2_generic(const int* idx, const double* w_left,
                           const double* w_right, int n, const double* block,
                           int h, double* f, double* g, double* info) const {
    double fa = 0.0;

    for (int k = 0; k < n; k++) {
      int i = idx[k];
      const double* e = block + static_cast<std::size_t>(k) * H;
      double a;
      double b;
      score_info(i, e, h, &a, &b);
      fa += logdens(i, e);
      double wl = w_left[k];
      double wr = w_right[k];
      g[0] += wl * a;
      g[1] += wr * a;
      info[0] += wl * wl * b;
      info[1] += wl * wr * b;
      info[2] += wr * wr * b;
    }

    *f = fa;
  }

  // A location the recorded draw should be shifted by, one per additive
  // predictor.
  //
  // This exists for the ordinal families, which are identified only up to a
  // common shift of the cutpoints and the predictor: adding a constant to every
  // cutpoint and to eta leaves every category probability alone. The sampler
  // works in the chart that pins the first cutpoint at zero, which is well
  // conditioned and keeps a two-category response identical to binary
  // regression. That is a poor chart to *report* in, because it makes the
  // cutpoints incomparable with `polr()`'s and leaves the reader to remember
  // that the predictor carries an offset. So each retained draw is recorded in
  // the chart where the predictor has mean zero over the fitted sample and every
  // cutpoint is free, which is what `polr()` reports.
  //
  // It is a change of chart, not of model: the shift is applied to the recorded
  // predictor, to the recorded cutpoints, and to the recorded leaf values, so
  // the stored forest replays to the shifted predictor and every identified
  // quantity is untouched.
  virtual arma::vec report_shift(const arma::mat& eta) const {
    return arma::zeros<arma::vec>(H);
  }

  // The nuisance parameters as they should be recorded, given that shift. Only a
  // family whose nuisance parameters live on the predictor's scale needs to do
  // anything here.
  virtual arma::vec aux_values_shifted(const arma::vec& shift) const {
    return aux_values();
  }

  // The log likelihood to report to the caller. It is the sampler's target for
  // every family but the augmented probit, where the target is the augmented
  // density and the quantity of interest is the one it is a device for.
  virtual double reported_loglik(const arma::mat& eta) const {
    return total_loglik(eta);
  }

  // The actual log likelihood, including the terms the sampler is free to drop.
  // Those terms are constant in eta, so adding them here cannot disturb any
  // acceptance ratio that compares two values of the predictor.
  double total_loglik(const arma::mat& eta) const {
    double out = logdens_extra_total();
    all_units();
    unit_values.resize(N);
    logdens_block(unit_idx.data(), N, eta.memptr(), unit_values.data());
    for (int i = 0; i < N; i++) {
      out += unit_values[i];
    }
    return out;
  }

  // The change in the log likelihood when component h of the predictor is
  // replaced by `new_h`, leaving the other components alone.
  //
  // The bandwidth move needs exactly this and nothing else. Getting it from two
  // calls to total_loglik() traversed the whole predictor matrix twice, built an
  // index vector of every observation twice, and summed a vector of N values
  // twice, to produce two numbers whose difference is all that is used. Here the
  // two evaluations of each observation happen while its predictors are in
  // cache, and the eta-free terms cancel rather than being computed.
  //
  // Blocked so that a family which overrides logdens_block() -- one whose log
  // density is an R function -- still goes through its own implementation. The
  // observations are visited in index order either way, so the sum is the one
  // the two-call version produced.
  virtual double loglik_delta(const arma::mat& eta, int h,
                              const double* new_h) const {
    const int CHUNK = 256;
    delta_idx.resize(CHUNK);
    delta_block.resize(static_cast<std::size_t>(CHUNK) * H);
    delta_old.resize(CHUNK);
    delta_new.resize(CHUNK);

    const double* e = eta.memptr();
    double out = 0.0;

    for (int start = 0; start < N; start += CHUNK) {
      int n = std::min(CHUNK, N - start);
      const double* from = e + static_cast<std::size_t>(start) * H;

      for (int k = 0; k < n; k++) {
        delta_idx[k] = start + k;
      }

      logdens_block(delta_idx.data(), n, from, delta_old.data());

      std::memcpy(delta_block.data(), from,
                  sizeof(double) * static_cast<std::size_t>(n) * H);

      for (int k = 0; k < n; k++) {
        delta_block[static_cast<std::size_t>(k) * H + h] = new_h[start + k];
      }

      logdens_block(delta_idx.data(), n, delta_block.data(), delta_new.data());

      for (int k = 0; k < n; k++) {
        out += delta_new[k] - delta_old[k];
      }
    }

    return out;
  }

  // The log density of every observation at one draw of the predictors, which
  // is what predict(type = "density") reports.
  void logdens_all(const arma::mat& eta, double* out) const {
    all_units();
    logdens_block(unit_idx.data(), N, eta.memptr(), out);
    for (int i = 0; i < N; i++) {
      out[i] += log_norm_const(i);
      if (!eta_free.is_empty()) {
        out[i] += w(i) * eta_free(i);
      }
    }
  }

protected:
  arma::vec eta_free;
  double eta_free_total;
  double norm_const_total;

  static double clamp_info(double value) {
    if (!std::isfinite(value) || value < 0.0) {
      return 0.0;
    }
    return value;
  }

  // Both derivatives from one set of three evaluations, at a step chosen for
  // the second difference. The resulting gradient carries a truncation error of
  // order step^2, which is immaterial for a proposal and, crucially, is a
  // deterministic function of the current state, so reversibility is preserved.
  void score_info_numeric(int i, const double* eta, int h, double* d1,
                          double* d2) const {
    const double step = 1e-4;
    std::copy(eta, eta + H, scratch.begin());
    double mid = logdens_unit(i, eta);
    scratch[h] = eta[h] + step;
    double up = logdens_unit(i, scratch.data());
    scratch[h] = eta[h] - step;
    double down = logdens_unit(i, scratch.data());
    *d1 = 0.5 * (up - down) / step;
    *d2 = -(up - 2.0 * mid + down) / (step * step);
  }

  mutable std::vector<double> scratch;

  // 0, 1, ..., N - 1, for the whole-sample block calls above.
  void all_units() const {
    if (static_cast<int>(unit_idx.size()) != N) {
      unit_idx.resize(N);
      for (int i = 0; i < N; i++) {
        unit_idx[i] = i;
      }
    }
  }

  mutable std::vector<int> unit_idx;
  mutable std::vector<double> unit_values;

  // Working space for loglik_delta(), sized once on first use.
  mutable std::vector<int> delta_idx;
  mutable std::vector<double> delta_block;
  mutable std::vector<double> delta_old;
  mutable std::vector<double> delta_new;
};

// Families derive through this instead of from Family directly, which is the
// curiously recurring template pattern: the derived type is the template
// argument, so `static_cast<const Derived&>(*this)` has a concrete type and
// calls on it are not virtual. The two accumulators below are then the same
// loops as Family's, with the family's own log density and derivatives inlined
// into them.
//
// The `Derived::` qualification on each call is what forces static dispatch --
// without it the calls would go back through the vtable and the whole point
// would be lost.
template <typename Derived>
struct Concrete : Family {
  Concrete(const arma::vec& y_, const arma::vec& w_, int H_)
    : Family(y_, w_, H_) {}

  // `Weighted` says whether the node carries membership weights at all: a hard
  // rule sends an observation entirely one way, so a hard tree stores none and
  // `wt` is null. `Fused` says whether the predictors are read from `eta`
  // directly, shifting component h by `wt * shift`, rather than from a block the
  // caller has already written them into. Both instantiations of each compute
  // the identical expression -- multiplying by 1.0 is exact, and the fused form
  // holds the shifted value in a register instead of storing and reloading it --
  // so all four agree to the last bit, which the test suite checks.
  // The three sums over a node's support, from a block the caller has already
  // written the predictors into.
  template <bool Weighted>
  void accumulate1_block(const int* idx, const double* wt, int n,
                         const double* block, int h, double* f, double* s,
                         double* j) const {
    const Derived& self = static_cast<const Derived&>(*this);
    double fa = 0.0;
    double sa = 0.0;
    double ja = 0.0;

    for (int k = 0; k < n; k++) {
      int i = idx[k];
      const double* e = block + static_cast<std::size_t>(k) * H;
      double a;
      double b;
      self.Derived::score_info_unit(i, e, h, &a, &b);
      double weight = w(i);
      // The same three products, in the same order, as Family::score_info() and
      // Family::logdens() would have formed them.
      double da = weight * a;
      double db = clamp_info(weight * b);
      fa += weight * self.Derived::logdens_unit(i, e);
      double wk = Weighted ? wt[k] : 1.0;
      sa += wk * da;
      ja += wk * wk * db;
    }

    *f = fa;
    *s = sa;
    *j = ja;
  }

  // The same sums without a block. `FromBase` picks the value of component h:
  // `base[k] + wt * amount` when the caller has the node's predictor with this
  // leaf's contribution already removed, and `eta[i] + wt * amount` when it does
  // not and the contribution is being shifted in place. With one additive
  // predictor the value never leaves a register; with several, the other
  // components are copied into a buffer of H rather than one of n * H.
  //
  // The expression is the one fill()/fill_at() wrote into the block, so this
  // agrees with accumulate1_block() to the last bit -- which the test suite
  // checks by running a whole fit each way.
  template <bool Weighted, bool FromBase>
  void accumulate1_fused(const int* idx, const double* wt, int n,
                         const double* eta_ptr, const double* base,
                         double amount, int h, double* f, double* s,
                         double* j) const {
    const Derived& self = static_cast<const Derived&>(*this);
    double fa = 0.0;
    double sa = 0.0;
    double ja = 0.0;
    double one;

    for (int k = 0; k < n; k++) {
      int i = idx[k];
      double wk = Weighted ? wt[k] : 1.0;
      const double* col = eta_ptr + static_cast<std::size_t>(i) * H;
      double value = FromBase ? base[k] + wk * amount : col[h] + wk * amount;
      const double* e;

      if (H == 1) {
        one = value;
        e = &one;
      }
      else {
        std::copy(col, col + H, scratch.begin());
        scratch[h] = value;
        e = scratch.data();
      }

      double a;
      double b;
      self.Derived::score_info_unit(i, e, h, &a, &b);
      double weight = w(i);
      double da = weight * a;
      double db = clamp_info(weight * b);
      fa += weight * self.Derived::logdens_unit(i, e);
      sa += wk * da;
      ja += wk * wk * db;
    }

    *f = fa;
    *s = sa;
    *j = ja;
  }

  void accumulate1(const int* idx, const double* wt, int n,
                   const double* block, int h, double* f, double* s,
                   double* j) const override {
    if (wt == nullptr) {
      accumulate1_block<false>(idx, wt, n, block, h, f, s, j);
      return;
    }

    accumulate1_block<true>(idx, wt, n, block, h, f, s, j);
  }

  void accumulate1_at(const int* idx, const double* wt, int n,
                      const arma::mat& eta, int h, double shift, double* f,
                      double* s, double* j) const override {
    if (wt == nullptr) {
      accumulate1_fused<false, false>(idx, wt, n, eta.memptr(), nullptr, shift,
                                      h, f, s, j);
      return;
    }

    accumulate1_fused<true, false>(idx, wt, n, eta.memptr(), nullptr, shift, h,
                                   f, s, j);
  }

  void accumulate1_from(const int* idx, const double* wt, int n,
                        const arma::mat& eta, const double* base, double mu,
                        int h, double* f, double* s, double* j) const override {
    if (wt == nullptr) {
      accumulate1_fused<false, true>(idx, wt, n, eta.memptr(), base, mu, h, f,
                                     s, j);
      return;
    }

    accumulate1_fused<true, true>(idx, wt, n, eta.memptr(), base, mu, h, f, s,
                                  j);
  }

  // The whole-sample log density, statically dispatched for the same reason.
  // This one is not a leaf loop: it is what the bandwidth move and the reported
  // log likelihood run over every observation, and through the virtual
  // interface it was two virtual calls per observation per tree per sweep --
  // measured at half the cost of a soft-rule fit. The product is formed in the
  // same order as Family::logdens(), so the two agree to the last bit.
  void logdens_block(const int* idx, int n, const double* block,
                     double* out) const override {
    const Derived& self = static_cast<const Derived&>(*this);

    for (int k = 0; k < n; k++) {
      int i = idx[k];
      out[k] = w(i) * self.Derived::logdens_unit(
        i, block + static_cast<std::size_t>(k) * H);
    }
  }

  // `Fused` reads the predictors from `eta` and forms component h from `base`
  // and the two children's weights, rather than from a block. This is the one
  // that matters most: the birth, death and change moves are three quarters of a
  // hard-rule fit and all three go through here, whereas the leaf refresh --
  // which is what accumulate1 serves -- is a fifth of it.
  template <bool Fused>
  void accumulate2_loop(const int* idx, const double* w_left,
                        const double* w_right, int n, const double* block,
                        const double* eta_ptr, const double* base,
                        double mu_left, double mu_right, int h, double* f,
                        double* g, double* info) const {
    const Derived& self = static_cast<const Derived&>(*this);
    double fa = 0.0;
    double one;

    for (int k = 0; k < n; k++) {
      int i = idx[k];
      double wl = w_left[k];
      double wr = w_right[k];
      const double* e;

      if (Fused) {
        double value = base[k] + wl * mu_left + wr * mu_right;

        if (H == 1) {
          one = value;
          e = &one;
        }
        else {
          const double* col = eta_ptr + static_cast<std::size_t>(i) * H;
          std::copy(col, col + H, scratch.begin());
          scratch[h] = value;
          e = scratch.data();
        }
      }
      else {
        e = block + static_cast<std::size_t>(k) * H;
      }

      double a;
      double b;
      self.Derived::score_info_unit(i, e, h, &a, &b);
      double weight = w(i);
      double da = weight * a;
      double db = clamp_info(weight * b);
      fa += weight * self.Derived::logdens_unit(i, e);
      g[0] += wl * da;
      g[1] += wr * da;
      info[0] += wl * wl * db;
      info[1] += wl * wr * db;
      info[2] += wr * wr * db;
    }

    *f = fa;
  }

  void accumulate2(const int* idx, const double* w_left, const double* w_right,
                   int n, const double* block, int h, double* f, double* g,
                   double* info) const override {
    accumulate2_loop<false>(idx, w_left, w_right, n, block, nullptr, nullptr,
                            0.0, 0.0, h, f, g, info);
  }

  void accumulate2_from(const int* idx, const double* w_left,
                        const double* w_right, int n, const arma::mat& eta,
                        const double* base, double mu_left, double mu_right,
                        int h, double* f, double* g,
                        double* info) const override {
    accumulate2_loop<true>(idx, w_left, w_right, n, nullptr, eta.memptr(), base,
                           mu_left, mu_right, h, f, g, info);
  }
};

// Build a family from the name and options passed down from R.
Family* make_family(const std::string& name, const std::string& link,
                    const arma::vec& y, const arma::vec& w,
                    const Rcpp::List& opts);

// The same likelihood written as the margin of a Gaussian one, where that is
// possible and the data allow it, which makes the target quadratic in the
// predictor and the Laplace approximation exact. Returns null when there is no
// such rewriting, so the caller keeps what make_family() gave it. Only the
// sampler uses this: everything that reports a density works with the family as
// the caller asked for it.
// `enabled` names the engine families the caller wants rewritten, since the
// rewriting is worth it for some and not others; see genbart_control().
Family* augmented_family(const std::string& name, const std::string& link,
                         const arma::vec& y, const arma::vec& w,
                         const Rcpp::List& opts,
                         const std::vector<std::string>& enabled);

} // namespace genbart

#endif
