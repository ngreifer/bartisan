#include "family.h"
#include "polyagamma.h"
#include "slice.h"

using namespace Rcpp;

namespace genbart {

namespace {

// log(F(b) - F(a)) for a < b, evaluated in whichever tail keeps the ratio of
// the two terms away from one. Used by the cumulative-link ordinal families,
// where the naive difference loses all precision once an observation sits far
// into a tail.
double log_diff_logistic(double a, double b) {
  if (!(b > a)) {
    return R_NegInf;
  }
  // expit(b) - expit(a) == expit(b) * expit(-a) * (1 - exp(-(b - a)))
  return log_expit(b) + log_expit(-a) + log1mexp(b - a);
}

// The complementary log-log link's distribution function is the smallest extreme
// value one, F(z) = 1 - exp(-exp(z)), and its survivor is exp(-exp(z)) exactly.
// That makes the interval probability a difference of two survivors, which needs
// no cancellation-prone subtraction of numbers near one:
//
//   F(b) - F(a) = exp(-exp(a)) - exp(-exp(b))
//               = exp(-exp(a)) * (1 - exp(-(exp(b) - exp(a))))
//
// so the log is -exp(a) plus log1mexp of the difference of the two exponentials.
double log_diff_cloglog(double a, double b) {
  if (!(b > a)) {
    return R_NegInf;
  }

  double ea = std::exp(a);
  double eb = std::exp(b);

  if (!std::isfinite(eb)) {
    // The upper limit has saturated: the interval is the whole upper tail.
    return -ea;
  }

  return -ea + log1mexp(eb - ea);
}

double log_diff_normal(double a, double b) {
  if (!(b > a)) {
    return R_NegInf;
  }
  if (a > 0.0) {
    // Work with the upper tails, where Phi(-a) dominates Phi(-b).
    double log_hi = R::pnorm5(-a, 0.0, 1.0, 1, 1);
    double log_lo = R::pnorm5(-b, 0.0, 1.0, 1, 1);
    return log_hi + log1mexp(log_hi - log_lo);
  }
  double log_hi = R::pnorm5(b, 0.0, 1.0, 1, 1);
  double log_lo = R::pnorm5(a, 0.0, 1.0, 1, 1);
  return log_hi + log1mexp(log_hi - log_lo);
}

// Inverse Mills ratio phi(r) / Phi(-r), computed on the log scale so that it
// stays accurate in the far upper tail where both terms underflow.
double inv_mills(double r) {
  return std::exp(R::dnorm4(r, 0.0, 1.0, 1) - R::pnorm5(-r, 0.0, 1.0, 1, 1));
}

} // namespace

// ---------------------------------------------------------------------------
// Gaussian with identity link. Conditionally conjugate, so this family is not
// the point of the package, but it is the reference case for testing.
// ---------------------------------------------------------------------------

struct GaussianFamily : Concrete<GaussianFamily> {
  double sigma;
  double sigma_hat;

  double prec;

  GaussianFamily(const arma::vec& y_, const arma::vec& w_, double sigma_hat_)
    : Concrete<GaussianFamily>(y_, w_, 1), sigma(sigma_hat_), sigma_hat(sigma_hat_) {
    prec = std::pow(sigma, -2.0);
  }

  double logdens_unit(int i, const double* eta) const override {
    double r = y(i) - eta[0];
    return -0.5 * r * r * prec;
  }

  // The normalizing constant depends on the scale but not on eta.
  arma::vec compute_eta_free() const override {
    return arma::vec(N).fill(-0.5 * LN_2PI - std::log(sigma));
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return (y(i) - eta[0]) * prec;
  }

  double info_unit(int i, const double* eta, int h) const override {
    return prec;
  }

  // A Gaussian log density is -prec * (y - eta)^2 / 2, quadratic in eta.
  TargetForm target_form(int h) const override {
    return TARGET_QUADRATIC;
  }

  // Overridden not for shared arithmetic -- there is none to share here -- but
  // because the default calls dlogdens_unit and info_unit, and those are two
  // virtual calls per observation that the compiler cannot inline. On a family
  // whose derivatives are this cheap, that dispatch is most of the cost.
  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    *d1 = (y(i) - eta[0]) * prec;
    *d2 = prec;
  }

  void update_aux(const arma::mat& eta) override {
    arma::vec r = y - eta.row(0).t();
    double sse = arma::dot(w, arma::square(r));
    double n = arma::sum(w);
    double drawn = half_cauchy_update_precision_mh(sse, n, prec, sigma_hat);
    sigma = std::pow(drawn, -0.5);
    prec = drawn;
    refresh_eta_free();
  }

  std::vector<std::string> aux_names() const override {
    return std::vector<std::string>{"sigma"};
  }

  arma::vec aux_values() const override {
    return arma::vec{sigma};
  }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem > 0) {
      sigma = values(0);
      prec = std::pow(sigma, -2.0);
      refresh_eta_free();
    }
  }
};

// ---------------------------------------------------------------------------
// Binomial. y holds proportions and the prior weights hold the number of
// trials, so Bernoulli data is the special case of unit weights.
// ---------------------------------------------------------------------------

struct BinomialFamily : Family {
  enum Link { LOGIT, PROBIT, CLOGLOG };
  Link link;

  // Where the probit derivatives switch from the natural scale to logarithms.
  // At |eta| = 5 the smaller tail is 2.9e-7, so 1 - p still carries nine
  // significant digits, and beyond that the log form takes over.
  static constexpr double PROBIT_DIRECT = 5.0;

  BinomialFamily(const arma::vec& y_, const arma::vec& w_, Link link_)
    : Family(y_, w_, 1), link(link_) {}

  // The weights hold the number of trials, so the binomial coefficient is the
  // piece of the mass function that is not simply the weight times a per-unit
  // term.
  //
  // Written with log-gamma rather than Rf_lchoose() for two reasons. The count
  // of successes is recovered as trials * proportion, and that product is not
  // reliably an integer in floating point even when the inputs are, so lchoose
  // would round it and emit a warning per observation. And a fractional weight
  // -- a Bayesian bootstrap draw, say -- is not a trial count at all, for which
  // lchoose has no meaning; the log-gamma form is the natural continuous
  // extension and returns zero for binary data at any weight, as it should.
  double log_norm_const(int i) const override {
    double trials = w(i);
    double successes = trials * y(i);

    if (!(trials > 0.0)) {
      return 0.0;
    }

    return R::lgammafn(trials + 1.0) - R::lgammafn(successes + 1.0) -
      R::lgammafn(trials - successes + 1.0);
  }

  // log Phi(eta) and log Phi(-eta) together. Only the smaller of the two is
  // taken from pnorm; the other follows from log1mexp, which is accurate
  // precisely where the probability is not small. That halves the number of
  // pnorm calls, and pnorm is by far the most expensive thing this family does:
  // on a probit fit it accounted for a factor of eleven in the cost of visiting
  // an observation, against the same fit with a logit link.
  static void log_probit_tails(double e, double* lp, double* lq) {
    if (e <= 0.0) {
      *lp = R::pnorm5(e, 0.0, 1.0, 1, 1);
      *lq = log1mexp(-(*lp));
    }
    else {
      *lq = R::pnorm5(-e, 0.0, 1.0, 1, 1);
      *lp = log1mexp(-(*lq));
    }
  }

  // log of the standard normal density, inlined rather than called: it is two
  // multiplications and a subtraction, and dnorm4 is a function call with
  // argument checking.
  static double log_dnorm(double e) {
    return -0.5 * (e * e + LN_2PI);
  }

  // Almost every binomial response is binary, and then one of the two terms is
  // multiplied by zero. Computing it anyway doubles the transcendental calls,
  // which is most of what this function costs.
  double logdens_unit(int i, const double* eta) const override {
    double e = eta[0];
    double p = y(i);
    bool one = p == 1.0;
    bool zero = p == 0.0;

    switch (link) {
    case LOGIT:
      if (one) {
        return log_expit(e);
      }
      if (zero) {
        return log1m_expit(e);
      }
      return p * log_expit(e) + (1.0 - p) * log1m_expit(e);
    case PROBIT: {
      if (one) {
        return R::pnorm5(e, 0.0, 1.0, 1, 1);
      }
      if (zero) {
        return R::pnorm5(-e, 0.0, 1.0, 1, 1);
      }
      double lp;
      double lq;
      log_probit_tails(e, &lp, &lq);
      return p * lp + (1.0 - p) * lq;
    }
    default: {
      // p(success) = 1 - exp(-exp(eta))
      double t = std::exp(e);
      if (zero) {
        return -t;
      }
      if (one) {
        return log1mexp(t);
      }
      return p * log1mexp(t) - (1.0 - p) * t;
    }
    }
  }

  // The score and the expected information share almost everything they need,
  // and the sampler always wants both at the same point, so they are computed
  // together. The default in Family computes them separately, which for probit
  // meant six calls to pnorm and dnorm where three suffice -- and, with the
  // tail trick above, where one does.
  //
  // The information is the expected Fisher information (dp/deta)^2 / (p(1 - p)),
  // strictly positive for every link, rather than the observed second
  // derivative, which for probit and cloglog can be negative.
  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double e = eta[0];

    switch (link) {
    case LOGIT: {
      double p = expit(e);
      *d1 = y(i) - p;
      *d2 = p * (1.0 - p);
      return;
    }
    case PROBIT: {
      // The inverse Mills ratios phi/Phi(eta) and phi/Phi(-eta). Away from the
      // tails they come from one call to pnorm on the natural scale and a
      // division, with no logarithms at all; only where one tail is small
      // enough for 1 - p to lose its digits is the log form needed.
      double a;
      double b;

      if (std::fabs(e) <= PROBIT_DIRECT) {
        double phi = std::exp(log_dnorm(e));
        double p = R::pnorm5(e, 0.0, 1.0, 1, 0);
        a = phi / p;
        b = phi / (1.0 - p);
      }
      else {
        double lp;
        double lq;
        log_probit_tails(e, &lp, &lq);
        double log_phi = log_dnorm(e);
        a = std::exp(log_phi - lp);
        b = std::exp(log_phi - lq);
      }

      *d1 = y(i) * a - (1.0 - y(i)) * b;
      *d2 = a * b;
      return;
    }
    default: {
      // t = exp(eta), so log(t) is eta and needs no log call.
      double t = std::exp(e);
      double l1m = log1mexp(t);
      double num = std::exp(e - t - l1m);
      *d1 = y(i) * num - (1.0 - y(i)) * t;
      *d2 = t * num;
      return;
    }
    }
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    double d1;
    double d2;
    score_info_unit(i, eta, h, &d1, &d2);
    return d1;
  }

  double info_unit(int i, const double* eta, int h) const override {
    double d1;
    double d2;
    score_info_unit(i, eta, h, &d1, &d2);
    return d2;
  }
};

// ---------------------------------------------------------------------------
// Probit regression by data augmentation (Albert and Chib 1993), for Bernoulli
// data.
//
// A latent z_i drawn from a normal with mean eta_i and unit variance, truncated
// to the positive half-line when y_i = 1 and to the negative half-line when
// y_i = 0, has the property that P(z_i > 0) = Phi(eta_i). So the probit
// likelihood is the margin of a Gaussian one, and conditional on z the model is
// a Gaussian BART with known variance.
//
// The point is not the augmentation itself but what it buys. A Gaussian log
// density is quadratic in the predictor, so the Laplace approximation stops
// being an approximation: Fisher scoring lands on the mode in one step, the
// leaf refresh becomes a Gibbs step with acceptance one, and the per-observation
// cost falls from a call to pnorm to a squared difference. This is exactly why
// dbarts fits a probit model as fast as a Gaussian one.
//
// What it costs is mixing. The augmented chain has to move z and eta in
// alternation, and for probabilities near zero or one that is known to be slow
// -- the reason this is an option rather than the default, and why the choice
// should be made on effective sample size per second rather than on time alone.
//
// Restricted to Bernoulli data. With w(i) trials the augmentation needs one
// latent per trial, and there is no single latent that carries the same
// information, so a binomial-count response uses the direct probit family.
// ---------------------------------------------------------------------------

struct ProbitAugmentedFamily : Concrete<ProbitAugmentedFamily> {
  arma::vec success;   // the observed response, 0 or 1
  arma::vec latent;    // z, redrawn every sweep

  ProbitAugmentedFamily(const arma::vec& y_, const arma::vec& w_)
    : Concrete<ProbitAugmentedFamily>(arma::vec(y_.n_elem, arma::fill::zeros), w_, 1), success(y_) {
    latent.set_size(N);
    // A deterministic start, replaced by a proper draw at the end of the first
    // sweep. Only the first few iterations of warmup ever see it.
    for (int i = 0; i < N; i++) {
      latent(i) = success(i) > 0.5 ? 0.5 : -0.5;
    }
    y = latent;
  }

  static bool applies(const arma::vec& y_, const arma::vec& w_) {
    for (arma::uword i = 0; i < y_.n_elem; i++) {
      if (w_(i) != 1.0) {
        return false;
      }
      if (y_(i) != 0.0 && y_(i) != 1.0) {
        return false;
      }
    }
    return true;
  }

  TargetForm target_form(int h) const override {
    return TARGET_QUADRATIC;
  }

  double logdens_unit(int i, const double* eta) const override {
    double r = latent(i) - eta[0];
    return -0.5 * r * r;
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    *d1 = latent(i) - eta[0];
    *d2 = 1.0;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return latent(i) - eta[0];
  }

  double info_unit(int i, const double* eta, int h) const override {
    return 1.0;
  }

  arma::vec compute_eta_free() const override {
    return arma::vec(N).fill(-0.5 * LN_2PI);
  }

  // A standard normal truncated to (lower, infinity), by inverting the upper
  // tail on the log scale so that a tail probability too small to represent
  // does not collapse the draw.
  static double truncated_normal_above(double lower) {
    double log_tail = R::pnorm5(lower, 0.0, 1.0, 0, 1);
    double u = unif_rand();

    if (u <= 0.0) {
      u = std::numeric_limits<double>::min();
    }

    return R::qnorm5(std::log(u) + log_tail, 0.0, 1.0, 0, 1);
  }

  void update_aux(const arma::mat& eta) override {
    const arma::rowvec& e = eta.row(0);

    for (int i = 0; i < N; i++) {
      // z - eta is standard normal truncated to the side y sits on.
      double shifted = (success(i) > 0.5)
        ? truncated_normal_above(-e(i))
        : -truncated_normal_above(e(i));
      latent(i) = e(i) + shifted;
    }

    y = latent;
  }

  // The sampler's target is the augmented density, which is what total_loglik
  // has to return, but the number worth reporting is the probit log likelihood
  // the augmentation is a device for.
  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      out += (success(i) > 0.5)
        ? R::pnorm5(e(i), 0.0, 1.0, 1, 1)
        : R::pnorm5(-e(i), 0.0, 1.0, 1, 1);
    }

    return out;
  }
};

// The category probability and the cutpoint update are the same problem whether
// the model is fitted directly or through its latent normal, so they live here
// rather than in either family. Category `cat` in 0, ..., num_cat - 1 occupies
// (c[cat - 1], c[cat]), with the outer edges at -Inf and +Inf.
// Which cumulative link the ordinal families use. Kept as an enum rather than a
// pair of booleans because there are now three of them and the wrong default
// would be a silent change of model.
enum OrdinalLink { ORD_LOGIT = 0, ORD_PROBIT = 1, ORD_CLOGLOG = 2 };

inline double ordinal_log_prob(int cat, double eta, const arma::vec& c,
                               int num_cat, int link) {
  bool has_lower = cat > 0;
  bool has_upper = cat < num_cat - 1;

  if (!has_lower && !has_upper) {
    return 0.0;
  }
  if (!has_lower) {
    double z = c(cat) - eta;

    if (link == ORD_LOGIT) {
      return log_expit(z);
    }
    if (link == ORD_PROBIT) {
      return R::pnorm5(z, 0.0, 1.0, 1, 1);
    }
    // log(1 - exp(-exp(z))).
    return log1mexp(std::exp(z));
  }
  if (!has_upper) {
    double z = c(cat - 1) - eta;

    if (link == ORD_LOGIT) {
      return log1m_expit(z);
    }
    if (link == ORD_PROBIT) {
      return R::pnorm5(-z, 0.0, 1.0, 1, 1);
    }
    // The complementary log-log survivor is exp(-exp(z)), so its log is exact.
    return -std::exp(z);
  }

  double lo = c(cat - 1) - eta;
  double hi = c(cat) - eta;

  if (link == ORD_LOGIT) {
    return log_diff_logistic(lo, hi);
  }
  if (link == ORD_PROBIT) {
    return log_diff_normal(lo, hi);
  }
  return log_diff_cloglog(lo, hi);
}

// cuts(0) stays at zero for identifiability; each remaining cutpoint is drawn
// between its neighbors, so the ordering constraint holds by construction.
//
// Cutpoint k enters the likelihood only through categories k and k + 1, being
// the upper limit of one and the lower limit of the other. Summing over just
// those two groups makes a sweep over every cutpoint cost O(n) in total rather
// than O(n * num_cat), because the group sizes add to n twice over. That is the
// same sparsity rms::orm() exploits in its information matrix, and it is what
// makes as many cutpoints as observations affordable.
//
// This is the *marginal* update: it uses the ordinal likelihood itself, not the
// latent normals, even when the rest of the sampler is running on them. Drawing
// a cutpoint from its full conditional given the latent variables -- uniform
// between the two order statistics that bracket it, as in Albert and Chib
// (1993) -- mixes badly once n is large, because those order statistics pin it
// to an interval of width O(1/n). Marginalizing the latent variables out of this
// one step and redrawing them straight afterwards is a partially collapsed Gibbs
// sampler in the sense of Van Dyk and Park (2008), and is the standard fix
// (Cowles 1996).
inline void update_ordinal_cuts(arma::vec& cuts, int num_cat, int link,
                                const std::vector<std::vector<int> >& by_cat,
                                const arma::vec& w, const arma::rowvec& e) {
  if (num_cat < 3) {
    return;
  }

  arma::vec trial = cuts;

  for (int k = 1; k < num_cat - 1; k++) {
    double lower = cuts(k - 1);
    double upper = (k + 1 <= num_cat - 2) ? cuts(k + 1) : 1e4;

    const std::vector<int>& below = by_cat[k];
    const std::vector<int>& above = by_cat[k + 1];

    auto logf = [&](double value) {
      trial(k) = value;
      double out = 0.0;
      for (std::size_t m = 0; m < below.size(); m++) {
        int i = below[m];
        out += w(i) * ordinal_log_prob(k, e(i), trial, num_cat, link);
      }
      for (std::size_t m = 0; m < above.size(); m++) {
        int i = above[m];
        out += w(i) * ordinal_log_prob(k + 1, e(i), trial, num_cat, link);
      }
      return out;
    };

    // With many cutpoints the admissible gap is narrow, so start the slice at
    // the scale of the gap rather than at a fixed width.
    double span = upper - lower;
    double width = span < 1.0 ? std::max(0.5 * span, 1e-8) : 0.5;

    cuts(k) = slice_sampler(cuts(k), logf, width, lower, upper);
    trial(k) = cuts(k);
  }
}

// Observation indices grouped by category, which is what the cutpoint update
// walks.
inline std::vector<std::vector<int> > group_by_category(const arma::vec& y,
                                                       int num_cat) {
  std::vector<std::vector<int> > out(num_cat);

  for (arma::uword i = 0; i < y.n_elem; i++) {
    int cat = static_cast<int>(y(i));
    if (cat >= 0 && cat < num_cat) {
      out[cat].push_back(static_cast<int>(i));
    }
  }

  return out;
}

// ---------------------------------------------------------------------------
// Ordinal probit through its latent normal (Albert and Chib 1993).
//
// A latent z_i ~ N(eta_i, 1) is introduced and the category is read off from
// which interval between the cutpoints it lands in. Conditional on z the log
// density in eta is -(z_i - eta_i)^2 / 2, so the target over a leaf is exactly
// quadratic and the sampler takes the closed form: one pass over a node instead
// of Fisher scoring plus a Metropolis ratio, and the two cumulative-normal
// evaluations per observation per pass disappear entirely. Measured on a
// thousand observations and fifty trees, that is the difference between an
// ordinal fit costing twenty times a Gaussian one and costing about twice.
//
// The cutpoints are *not* drawn from their conditional given z; see
// update_ordinal_cuts() for why, and for what is done instead.
// ---------------------------------------------------------------------------

struct OrdinalProbitAugmentedFamily : Concrete<OrdinalProbitAugmentedFamily> {
  int num_cat;
  arma::vec cuts;       // length num_cat - 1, cuts(0) fixed at 0
  bool update_cuts;
  arma::vec cat;        // the observed category, 0 ... num_cat - 1
  arma::vec latent;     // z, redrawn every sweep
  std::vector<std::vector<int> > by_cat;

  OrdinalProbitAugmentedFamily(const arma::vec& y_, const arma::vec& w_,
                               int num_cat_, const arma::vec& cuts_,
                               bool update_cuts_)
    : Concrete<OrdinalProbitAugmentedFamily>(
        arma::vec(y_.n_elem, arma::fill::zeros), w_, 1),
      num_cat(num_cat_), cuts(cuts_), update_cuts(update_cuts_), cat(y_) {
    by_cat = group_by_category(y_, num_cat_);
    latent.set_size(N);

    // A deterministic start, replaced by a proper draw at the end of the first
    // sweep: the midpoint of the observation's interval, with the open ends
    // taken one unit past the cutpoint that bounds them.
    for (int i = 0; i < N; i++) {
      int k = static_cast<int>(cat(i));
      double lo = k > 0 ? cuts(k - 1) : cuts(0) - 2.0;
      double hi = k < num_cat - 1 ? cuts(k) : cuts(num_cat - 2) + 2.0;
      latent(i) = 0.5 * (lo + hi);
    }

    y = latent;
  }

  static bool applies(const arma::vec& w_) {
    for (arma::uword i = 0; i < w_.n_elem; i++) {
      if (w_(i) != 1.0) {
        return false;
      }
    }
    return true;
  }

  TargetForm target_form(int h) const override {
    return TARGET_QUADRATIC;
  }

  double logdens_unit(int i, const double* eta) const override {
    double r = latent(i) - eta[0];
    return -0.5 * r * r;
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    *d1 = latent(i) - eta[0];
    *d2 = 1.0;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return latent(i) - eta[0];
  }

  double info_unit(int i, const double* eta, int h) const override {
    return 1.0;
  }

  arma::vec compute_eta_free() const override {
    return arma::vec(N).fill(-0.5 * LN_2PI);
  }

  void update_aux(const arma::mat& eta) override {
    const arma::rowvec& e = eta.row(0);

    // The cutpoints first, from the ordinal likelihood with the latent normals
    // integrated out, and the latent normals immediately afterwards from their
    // conditional given the drawn cutpoints.
    if (update_cuts) {
      update_ordinal_cuts(cuts, num_cat, ORD_PROBIT, by_cat, w, e);
    }

    for (int i = 0; i < N; i++) {
      int k = static_cast<int>(cat(i));
      double lo = k > 0 ? cuts(k - 1) - e(i) : R_NegInf;
      double hi = k < num_cat - 1 ? cuts(k) - e(i) : R_PosInf;
      latent(i) = e(i) + truncated_normal_between(lo, hi);
    }

    y = latent;
  }

  // The sampler's target is the augmented density; the number worth reporting is
  // the ordinal log likelihood the augmentation is a device for.
  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      out += w(i) * ordinal_log_prob(static_cast<int>(cat(i)), e(i), cuts,
                                     num_cat, ORD_PROBIT);
    }

    return out;
  }


  // Identified only up to a common shift of the cutpoints and the predictor, so
  // the draw is recorded in the chart where the predictor has mean zero over the
  // fitted sample and every cutpoint is free.
  //
  // Two categories are the exception, and deliberately so: there the model *is*
  // binary regression, the single boundary is conventionally folded into the
  // intercept, and reporting it as a free cutpoint against a centered predictor
  // would put the same fit on a different scale from `binomial()`.
  arma::vec report_shift(const arma::mat& eta) const override {
    arma::vec out(1, arma::fill::zeros);

    if (num_cat > 2) {
      out(0) = arma::mean(eta.row(0));
    }

    return out;
  }

  arma::vec aux_values_shifted(const arma::vec& shift) const override {
    return cuts - shift(0);
  }
  std::vector<std::string> aux_names() const override {
    std::vector<std::string> out;
    for (int k = 0; k < num_cat - 1; k++) {
      out.push_back("cut" + std::to_string(k + 1));
    }
    return out;
  }

  arma::vec aux_values() const override { return cuts; }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem == cuts.n_elem) {
      cuts = values;
    }
  }
};

// ---------------------------------------------------------------------------
// Ordinal logit through a normal scale mixture.
//
// The cumulative logit has no latent normal of its own -- its latent variable is
// logistic -- but a logistic variate is a normal whose precision is itself
// random, so introducing that precision buys the same thing the probit's latent
// normal buys: a target quadratic in the additive predictor.
//
// The mixing distribution is usually reached through the Kolmogorov-Smirnov
// density (Andrews and Mallows 1974; Holmes and Held 2006), which needs a
// sampler of its own. It does not have to be. Polson, Scott and Windle (2013,
// Theorem 1) with a = 1 and b = 2 reads
//
//     e^x / (1 + e^x)^2  =  (1/4) E[ exp(-w x^2 / 2) ],   w ~ PG(2, 0)
//
// and the left-hand side is exactly the standard logistic density. So a logistic
// residual is a mean-zero normal with precision w, mixed over that density -- and
// because PG(b, c) is PG(b, 0) tilted by exp(-c^2 w / 2), the conditional of w
// given a residual r is exactly PG(2, |r|). That is a Polya-Gamma draw at
// integer b, which the exact Devroye sampler in polyagamma.cpp already covers.
// Nothing approximate enters anywhere.
//
// The cutpoints are drawn as in the probit case: from the ordinal likelihood with
// the latent variables integrated out, with the latent variables redrawn
// straight afterwards.
// ---------------------------------------------------------------------------

struct OrdinalLogitAugmentedFamily : Concrete<OrdinalLogitAugmentedFamily> {
  int num_cat;
  arma::vec cuts;       // length num_cat - 1, cuts(0) fixed at 0
  bool update_cuts;
  arma::vec cat;        // the observed category, 0 ... num_cat - 1
  arma::vec latent;     // z
  arma::vec omega;      // the precision of z, one per observation
  std::vector<std::vector<int> > by_cat;

  // A precision this small makes the latent normal wider than any cutpoint gap
  // could need, and guards the division below.
  static constexpr double OMEGA_MIN = 1e-10;

  OrdinalLogitAugmentedFamily(const arma::vec& y_, const arma::vec& w_,
                              int num_cat_, const arma::vec& cuts_,
                              bool update_cuts_)
    : Concrete<OrdinalLogitAugmentedFamily>(
        arma::vec(y_.n_elem, arma::fill::zeros), w_, 1),
      num_cat(num_cat_), cuts(cuts_), update_cuts(update_cuts_), cat(y_) {
    by_cat = group_by_category(y_, num_cat_);
    latent.set_size(N);
    omega.set_size(N);

    // A deterministic start, replaced by a proper draw at the end of the first
    // sweep. The precision starts at the mean of PG(2, 0), which is 1/2.
    for (int i = 0; i < N; i++) {
      int k = static_cast<int>(cat(i));
      double lo = k > 0 ? cuts(k - 1) : cuts(0) - 2.0;
      double hi = k < num_cat - 1 ? cuts(k) : cuts(num_cat - 2) + 2.0;
      latent(i) = 0.5 * (lo + hi);
      omega(i) = 0.5;
    }

    y = latent;
  }

  static bool applies(const arma::vec& w_) {
    for (arma::uword i = 0; i < w_.n_elem; i++) {
      if (w_(i) != 1.0) {
        return false;
      }
    }
    return true;
  }

  TargetForm target_form(int h) const override {
    return TARGET_QUADRATIC;
  }

  double logdens_unit(int i, const double* eta) const override {
    double r = latent(i) - eta[0];
    return -0.5 * omega(i) * r * r;
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    *d1 = omega(i) * (latent(i) - eta[0]);
    *d2 = omega(i);
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return omega(i) * (latent(i) - eta[0]);
  }

  double info_unit(int i, const double* eta, int h) const override {
    return omega(i);
  }

  void update_aux(const arma::mat& eta) override {
    const arma::rowvec& e = eta.row(0);

    if (update_cuts) {
      update_ordinal_cuts(cuts, num_cat, ORD_LOGIT, by_cat, w, e);
    }

    for (int i = 0; i < N; i++) {
      // The precision, given the residual the previous sweep left behind.
      double r = latent(i) - e(i);
      omega(i) = std::max(rpg(2.0, std::fabs(r)), OMEGA_MIN);

      // Then the latent variable, given that precision and the cutpoints.
      double sd = 1.0 / std::sqrt(omega(i));
      int k = static_cast<int>(cat(i));
      double lo = k > 0 ? (cuts(k - 1) - e(i)) / sd : R_NegInf;
      double hi = k < num_cat - 1 ? (cuts(k) - e(i)) / sd : R_PosInf;
      latent(i) = e(i) + sd * truncated_normal_between(lo, hi);
    }

    y = latent;
  }

  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      out += w(i) * ordinal_log_prob(static_cast<int>(cat(i)), e(i), cuts,
                                     num_cat, ORD_LOGIT);
    }

    return out;
  }


  // Identified only up to a common shift of the cutpoints and the predictor, so
  // the draw is recorded in the chart where the predictor has mean zero over the
  // fitted sample and every cutpoint is free.
  //
  // Two categories are the exception, and deliberately so: there the model *is*
  // binary regression, the single boundary is conventionally folded into the
  // intercept, and reporting it as a free cutpoint against a centered predictor
  // would put the same fit on a different scale from `binomial()`.
  arma::vec report_shift(const arma::mat& eta) const override {
    arma::vec out(1, arma::fill::zeros);

    if (num_cat > 2) {
      out(0) = arma::mean(eta.row(0));
    }

    return out;
  }

  arma::vec aux_values_shifted(const arma::vec& shift) const override {
    return cuts - shift(0);
  }
  std::vector<std::string> aux_names() const override {
    std::vector<std::string> out;
    for (int k = 0; k < num_cat - 1; k++) {
      out.push_back("cut" + std::to_string(k + 1));
    }
    return out;
  }

  arma::vec aux_values() const override { return cuts; }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem == cuts.n_elem) {
      cuts = values;
    }
  }
};

// ---------------------------------------------------------------------------
// Ordinal complementary log-log through a latent exponential time.
//
// The cumulative cloglog model is the discrete proportional hazards model, and
// that is what gives it a latent variable. Writing the survivor as
//
//     P(Y > k) = exp(-exp(c_k - eta)) = exp(-L_k * exp(-eta)),   L_k = exp(c_k)
//
// says exactly that a latent T with an exponential distribution of rate
// exp(-eta) exceeds L_k. So the category is which interval between the
// transformed cutpoints T lands in, and conditional on T the log density in eta
// is
//
//     -eta - T * exp(-eta)
//
// which is TARGET_EXP_DOWN -- the same shape as the gamma family and the
// augmented negative binomial. Three numbers from one pass determine the whole
// function under hard rules, and even under soft rules the per-observation
// density is one exp() instead of a difference of two extreme-value
// distribution functions.
//
// The latent time is drawn from its exponential distribution truncated to the
// observed interval, by inverting the distribution function on a scale where
// neither endpoint has to be formed as exp(c_k): the two quantities the draw
// needs are exp(c_{k-1} - eta) and exp(c_k - eta), which are what the log
// density computes anyway.
//
// The cutpoints are drawn from the marginal ordinal likelihood with the latent
// times integrated out, for the same reason as in the other two ordinal
// augmentations: their conditional given the latent times is uniform between two
// order statistics, which pins them to an interval of width O(1/n).
// ---------------------------------------------------------------------------

struct OrdinalCloglogAugmentedFamily : Concrete<OrdinalCloglogAugmentedFamily> {
  int num_cat;
  arma::vec cuts;
  bool update_cuts;
  arma::vec cat;        // the observed category, 0 ... num_cat - 1
  arma::vec latent;     // T, redrawn every sweep
  std::vector<std::vector<int> > by_cat;

  OrdinalCloglogAugmentedFamily(const arma::vec& y_, const arma::vec& w_,
                                int num_cat_, const arma::vec& cuts_,
                                bool update_cuts_)
    : Concrete<OrdinalCloglogAugmentedFamily>(
        arma::vec(y_.n_elem, arma::fill::zeros), w_, 1),
      num_cat(num_cat_), cuts(cuts_), update_cuts(update_cuts_), cat(y_) {
    by_cat = group_by_category(y_, num_cat_);
    latent.set_size(N);

    // A deterministic start, replaced by a proper draw at the end of the first
    // sweep: the geometric midpoint of the observation's interval, with the open
    // ends taken a factor of e away from the cutpoint that bounds them.
    for (int i = 0; i < N; i++) {
      int k = static_cast<int>(cat(i));
      double lo = k > 0 ? cuts(k - 1) : cuts(0) - 1.0;
      double hi = k < num_cat - 1 ? cuts(k) : cuts(num_cat - 2) + 1.0;
      latent(i) = std::exp(0.5 * (lo + hi));
    }

    y = latent;
  }

  static bool applies(const arma::vec& w_) {
    for (arma::uword i = 0; i < w_.n_elem; i++) {
      if (w_(i) != 1.0) {
        return false;
      }
    }
    return true;
  }

  TargetForm target_form(int h) const override {
    return TARGET_EXP_DOWN;
  }

  double logdens_unit(int i, const double* eta) const override {
    return -eta[0] - latent(i) * std::exp(-eta[0]);
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double scaled = latent(i) * std::exp(-eta[0]);
    *d1 = scaled - 1.0;
    *d2 = scaled;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return latent(i) * std::exp(-eta[0]) - 1.0;
  }

  double info_unit(int i, const double* eta, int h) const override {
    return latent(i) * std::exp(-eta[0]);
  }

  void update_aux(const arma::mat& eta) override {
    const arma::rowvec& e = eta.row(0);

    if (update_cuts) {
      update_ordinal_cuts(cuts, num_cat, ORD_CLOGLOG, by_cat, w, e);
    }

    for (int i = 0; i < N; i++) {
      int k = static_cast<int>(cat(i));

      // The rate times each endpoint, which is exp(cut - eta) and never needs
      // exp(cut) on its own.
      double a = k > 0 ? std::exp(cuts(k - 1) - e(i)) : 0.0;
      double gap = k < num_cat - 1 ? std::exp(cuts(k) - e(i)) - a : R_PosInf;

      double u = unif_rand();

      if (u >= 1.0) {
        u = 1.0 - std::numeric_limits<double>::epsilon();
      }

      // Inverting the truncated exponential: with `a` and `a + gap` the rate
      // times the two endpoints, the draw is (a - log1p(-u * (1 - exp(-gap))))
      // divided by the rate. Written this way the argument of log1p stays in
      // (-1, 0] whatever the gap, including an infinite one.
      double shape = std::isinf(gap) ? a - std::log1p(-u)
                                     : a - std::log1p(u * std::expm1(-gap));

      double draw = shape * std::exp(e(i));

      latent(i) = std::isfinite(draw) && draw > 0.0 ? draw : std::exp(e(i));
    }

    y = latent;
  }

  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      out += w(i) * ordinal_log_prob(static_cast<int>(cat(i)), e(i), cuts,
                                     num_cat, ORD_CLOGLOG);
    }

    return out;
  }


  // Identified only up to a common shift of the cutpoints and the predictor, so
  // the draw is recorded in the chart where the predictor has mean zero over the
  // fitted sample and every cutpoint is free.
  //
  // Two categories are the exception, and deliberately so: there the model *is*
  // binary regression, the single boundary is conventionally folded into the
  // intercept, and reporting it as a free cutpoint against a centered predictor
  // would put the same fit on a different scale from `binomial()`.
  arma::vec report_shift(const arma::mat& eta) const override {
    arma::vec out(1, arma::fill::zeros);

    if (num_cat > 2) {
      out(0) = arma::mean(eta.row(0));
    }

    return out;
  }

  arma::vec aux_values_shifted(const arma::vec& shift) const override {
    return cuts - shift(0);
  }
  std::vector<std::string> aux_names() const override {
    std::vector<std::string> out;
    for (int k = 0; k < num_cat - 1; k++) {
      out.push_back("cut" + std::to_string(k + 1));
    }
    return out;
  }

  arma::vec aux_values() const override { return cuts; }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem == cuts.n_elem) {
      cuts = values;
    }
  }
};

// ---------------------------------------------------------------------------
// Poisson with log link.
// ---------------------------------------------------------------------------

struct PoissonFamily : Concrete<PoissonFamily> {
  PoissonFamily(const arma::vec& y_, const arma::vec& w_)
    : Concrete<PoissonFamily>(y_, w_, 1) {}

  double logdens_unit(int i, const double* eta) const override {
    return y(i) * eta[0] - std::exp(eta[0]);
  }

  arma::vec compute_eta_free() const override {
    arma::vec out(N);
    for (int i = 0; i < N; i++) {
      out(i) = -R::lgammafn(y(i) + 1.0);
    }
    return out;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return y(i) - std::exp(eta[0]);
  }

  double info_unit(int i, const double* eta, int h) const override {
    return std::exp(eta[0]);
  }

  // The mean is both the score's subtrahend and the information, so computing
  // them together saves an exponential as well as two virtual calls.
  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double mu = std::exp(eta[0]);
    *d1 = y(i) - mu;
    *d2 = mu;
  }

  // y * eta - exp(eta): linear plus an exponential.
  TargetForm target_form(int h) const override { return TARGET_EXP_UP; }
};

// ---------------------------------------------------------------------------
// Negative binomial with log link, in the mean/size parameterization so that
// mu = exp(eta) and var = mu + mu^2 / theta. theta is given a gamma prior and
// drawn by slice sampling on the log scale.
// ---------------------------------------------------------------------------

struct NegBinFamily : Family {
  double theta;
  double prior_shape;
  double prior_rate;
  bool update_theta;

  NegBinFamily(const arma::vec& y_, const arma::vec& w_, double theta_,
               double prior_shape_, double prior_rate_, bool update_theta_)
    : Family(y_, w_, 1), theta(theta_), prior_shape(prior_shape_),
      prior_rate(prior_rate_), update_theta(update_theta_) {
    log_theta = std::log(theta);
  }

  double log_theta;

  static double loglik_one(double y, double eta, double theta) {
    double mu = std::exp(eta);
    double log_theta_mu = std::log(theta + mu);
    return R::lgammafn(y + theta) - R::lgammafn(theta) - R::lgammafn(y + 1.0) +
      theta * (std::log(theta) - log_theta_mu) + y * (eta - log_theta_mu);
  }

  // Only the two terms involving mu vary with eta; the three log-gamma terms
  // depend on y and theta alone.
  double logdens_unit(int i, const double* eta) const override {
    double mu = std::exp(eta[0]);
    double log_theta_mu = std::log(theta + mu);
    return theta * (log_theta - log_theta_mu) + y(i) * (eta[0] - log_theta_mu);
  }

  arma::vec compute_eta_free() const override {
    arma::vec out(N);
    double lg_theta = R::lgammafn(theta);
    for (int i = 0; i < N; i++) {
      out(i) = R::lgammafn(y(i) + theta) - lg_theta - R::lgammafn(y(i) + 1.0);
    }
    return out;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    double mu = std::exp(eta[0]);
    return theta * (y(i) - mu) / (theta + mu);
  }

  double info_unit(int i, const double* eta, int h) const override {
    double mu = std::exp(eta[0]);
    return mu * theta / (mu + theta);
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double mu = std::exp(eta[0]);
    double ratio = theta / (theta + mu);
    *d1 = ratio * (y(i) - mu);
    *d2 = ratio * mu;
  }

  void update_aux(const arma::mat& eta) override {
    if (!update_theta) {
      return;
    }
    const arma::rowvec& e = eta.row(0);
    auto logf = [this, &e](double log_theta) {
      double th = std::exp(log_theta);
      if (!(th > 0.0) || !std::isfinite(th)) {
        return R_NegInf;
      }
      double out = prior_shape * log_theta - prior_rate * th;
      for (int i = 0; i < N; i++) {
        out += w(i) * loglik_one(y(i), e(i), th);
      }
      return out;
    };
    theta = std::exp(slice_sampler(std::log(theta), logf, 1.0, -20.0, 20.0));
    log_theta = std::log(theta);
    refresh_eta_free();
  }

  std::vector<std::string> aux_names() const override {
    return std::vector<std::string>{"theta"};
  }

  arma::vec aux_values() const override { return arma::vec{theta}; }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem > 0) {
      theta = values(0);
      log_theta = std::log(theta);
      refresh_eta_free();
    }
  }
};

// ---------------------------------------------------------------------------
// Gamma with log link, in the GLM parameterization: mean exp(eta) and constant
// shape, so that shape plays the role of the inverse dispersion.
// ---------------------------------------------------------------------------

struct GammaFamily : Concrete<GammaFamily> {
  double shape;
  double prior_shape;
  double prior_rate;
  bool update_shape;

  GammaFamily(const arma::vec& y_, const arma::vec& w_, double shape_,
              double prior_shape_, double prior_rate_, bool update_shape_)
    : Concrete<GammaFamily>(y_, w_, 1), shape(shape_), prior_shape(prior_shape_),
      prior_rate(prior_rate_), update_shape(update_shape_) {}

  static double loglik_one(double y, double eta, double shape) {
    return shape * (std::log(shape) - eta) - R::lgammafn(shape) +
      (shape - 1.0) * std::log(y) - shape * y * std::exp(-eta);
  }

  double logdens_unit(int i, const double* eta) const override {
    return -shape * eta[0] - shape * y(i) * std::exp(-eta[0]);
  }

  arma::vec compute_eta_free() const override {
    arma::vec out(N);
    double head = shape * std::log(shape) - R::lgammafn(shape);
    for (int i = 0; i < N; i++) {
      out(i) = head + (shape - 1.0) * std::log(y(i));
    }
    return out;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return shape * (y(i) * std::exp(-eta[0]) - 1.0);
  }

  // The observed second derivative, shape * y / mu, rather than the expected
  // information, which is just the shape. Most families that report the
  // expected version do so because the observed one can go negative; here it
  // cannot, since the response is strictly positive. Reporting the true
  // curvature makes the Laplace fit a genuine second-order match rather than a
  // Fisher-scoring one, and it is what the exponential form below has to
  // recover the target's coefficients from.
  double info_unit(int i, const double* eta, int h) const override {
    return shape * y(i) * std::exp(-eta[0]);
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double scaled = shape * y(i) * std::exp(-eta[0]);
    *d1 = scaled - shape;
    *d2 = scaled;
  }

  // -shape * eta - shape * y * exp(-eta): linear plus a decaying exponential.
  TargetForm target_form(int h) const override { return TARGET_EXP_DOWN; }

  void update_aux(const arma::mat& eta) override {
    if (!update_shape) {
      return;
    }
    const arma::rowvec& e = eta.row(0);
    auto logf = [this, &e](double log_shape) {
      double sh = std::exp(log_shape);
      if (!(sh > 0.0) || !std::isfinite(sh)) {
        return R_NegInf;
      }
      double out = prior_shape * log_shape - prior_rate * sh;
      for (int i = 0; i < N; i++) {
        out += w(i) * loglik_one(y(i), e(i), sh);
      }
      return out;
    };
    shape = std::exp(slice_sampler(std::log(shape), logf, 1.0, -20.0, 20.0));
    refresh_eta_free();
  }

  std::vector<std::string> aux_names() const override {
    return std::vector<std::string>{"shape"};
  }

  arma::vec aux_values() const override { return arma::vec{shape}; }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem > 0) {
      shape = values(0);
      refresh_eta_free();
    }
  }
};

// ---------------------------------------------------------------------------
// Ordinal cumulative-link model, parameterized as P(Y <= k) = F(c_k - eta) so
// that larger eta shifts mass towards higher categories, matching polr(). The
// first cutpoint is pinned at zero for identifiability, which makes the
// two-category case exactly binary regression with the matching link.
// ---------------------------------------------------------------------------

struct OrdinalFamily : Family {
  int link;
  int num_cat;
  arma::vec cuts;   // length num_cat - 1, cuts(0) fixed at 0
  bool update_cuts;

  // Observation indices grouped by category, so that a cutpoint update can
  // touch only the two groups that involve it.
  std::vector<std::vector<int>> by_cat;

  OrdinalFamily(const arma::vec& y_, const arma::vec& w_, int link_,
                int num_cat_, const arma::vec& cuts_, bool update_cuts_)
    : Family(y_, w_, 1), link(link_), num_cat(num_cat_), cuts(cuts_),
      update_cuts(update_cuts_) {
    by_cat = group_by_category(y_, num_cat_);
  }

  double log_prob(int cat, double eta, const arma::vec& c) const {
    return ordinal_log_prob(cat, eta, c, num_cat, link);
  }

  double logdens_unit(int i, const double* eta) const override {
    return log_prob(static_cast<int>(y(i)), eta[0], cuts);
  }

  // Density and its derivative for the chosen link, at z = cut - eta.
  void link_dens(double z, double* f, double* fp) const {
    if (link == ORD_LOGIT) {
      double p = expit(z);
      *f = p * (1.0 - p);
      *fp = *f * (1.0 - 2.0 * p);
    }
    else if (link == ORD_PROBIT) {
      *f = R::dnorm4(z, 0.0, 1.0, 0);
      *fp = -z * (*f);
    }
    else {
      // The smallest extreme value density, exp(z - exp(z)), whose derivative is
      // itself times (1 - exp(z)). Both underflow to zero rather than overflow,
      // which is the right behavior at the ends of the scale.
      double e = std::exp(z);
      *f = std::exp(z - e);
      *fp = *f * (1.0 - e);
    }
  }

  // With P = F(b) - F(a) and both limits shifted by -eta, dP/deta = f(a) - f(b),
  // so the score is (f(a) - f(b)) / P and the information follows from one more
  // differentiation. The interval probability of a log-concave density is
  // log-concave in eta, so the result is positive for both supported links.
  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    int cat = static_cast<int>(y(i));
    double e = eta[0];
    double prob = std::exp(log_prob(cat, e, cuts));

    if (!(prob > 1e-300) || !std::isfinite(prob)) {
      score_info_numeric(i, eta, h, d1, d2);
      return;
    }

    double f_lo = 0.0;
    double fp_lo = 0.0;
    double f_hi = 0.0;
    double fp_hi = 0.0;

    if (cat > 0) {
      link_dens(cuts(cat - 1) - e, &f_lo, &fp_lo);
    }
    if (cat < num_cat - 1) {
      link_dens(cuts(cat) - e, &f_hi, &fp_hi);
    }

    double diff = f_lo - f_hi;
    *d1 = diff / prob;
    *d2 = ((fp_lo - fp_hi) * prob + diff * diff) / (prob * prob);
  }

  void update_aux(const arma::mat& eta) override {
    if (!update_cuts) {
      return;
    }

    update_ordinal_cuts(cuts, num_cat, link, by_cat, w, eta.row(0));
  }


  // Identified only up to a common shift of the cutpoints and the predictor, so
  // the draw is recorded in the chart where the predictor has mean zero over the
  // fitted sample and every cutpoint is free.
  //
  // Two categories are the exception, and deliberately so: there the model *is*
  // binary regression, the single boundary is conventionally folded into the
  // intercept, and reporting it as a free cutpoint against a centered predictor
  // would put the same fit on a different scale from `binomial()`.
  arma::vec report_shift(const arma::mat& eta) const override {
    arma::vec out(1, arma::fill::zeros);

    if (num_cat > 2) {
      out(0) = arma::mean(eta.row(0));
    }

    return out;
  }

  arma::vec aux_values_shifted(const arma::vec& shift) const override {
    return cuts - shift(0);
  }
  std::vector<std::string> aux_names() const override {
    std::vector<std::string> out;
    for (int k = 0; k < num_cat - 1; k++) {
      out.push_back("cut" + std::to_string(k + 1));
    }
    return out;
  }

  arma::vec aux_values() const override { return cuts; }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem == cuts.n_elem) {
      cuts = values;
    }
  }
};

// ---------------------------------------------------------------------------
// Multinomial logistic, in either of two parameterizations.
//
// Reference coding carries num_cat - 1 forests and pins the first category's
// predictor at zero, so the fitted functions are log odds against that
// category and are interpretable the same way as the coefficients from
// nnet::multinom(). The cost is that the prior is not symmetric in the
// categories: which one is the reference changes the model.
//
// The symmetric coding of Murray (2021) carries one forest per category and
// leaves the model unidentified, since adding any function of x to every
// predictor leaves the probabilities alone. The proper leaf prior makes the
// posterior proper regardless, and every identified quantity -- probabilities,
// odds ratios -- is recovered by post-processing. The prior is then exchangeable
// over the categories, which is the point.
// ---------------------------------------------------------------------------

struct MultinomFamily : Family {
  int num_cat;
  bool symmetric;

  MultinomFamily(const arma::vec& y_, const arma::vec& w_, int num_cat_,
                 bool symmetric_)
    : Family(y_, w_, symmetric_ ? num_cat_ : num_cat_ - 1), num_cat(num_cat_),
      symmetric(symmetric_) {}

  // log sum_j exp(eta_j), shifting out the largest predictor first. Under
  // reference coding the omitted category contributes exp(0).
  double log_denom(const double* eta) const {
    double m = symmetric ? eta[0] : 0.0;
    for (int h = 0; h < H; h++) {
      if (eta[h] > m) {
        m = eta[h];
      }
    }
    double total = symmetric ? 0.0 : std::exp(-m);
    for (int h = 0; h < H; h++) {
      total += std::exp(eta[h] - m);
    }
    return m + std::log(total);
  }

  // Which predictor belongs to the observed category, or -1 when the category
  // is the reference and has no predictor of its own.
  int own(int i) const {
    int cat = static_cast<int>(y(i));
    return symmetric ? cat : cat - 1;
  }

  double logdens_unit(int i, const double* eta) const override {
    int j = own(i);
    return (j >= 0 ? eta[j] : 0.0) - log_denom(eta);
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    double p = std::exp(eta[h] - log_denom(eta));
    return (own(i) == h ? 1.0 : 0.0) - p;
  }

  double info_unit(int i, const double* eta, int h) const override {
    double p = std::exp(eta[h] - log_denom(eta));
    return p * (1.0 - p);
  }

  // log_denom is a pass over every category, and the default would run it
  // twice for one Fisher-scoring evaluation.
  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double p = std::exp(eta[h] - log_denom(eta));
    *d1 = (own(i) == h ? 1.0 : 0.0) - p;
    *d2 = p * (1.0 - p);
  }
};

// ---------------------------------------------------------------------------
// Accelerated failure time models for right-censored data. The response is
// log(time) and the model is log(T) = eta + sigma * epsilon, with the error
// distribution selecting the survival family: Gumbel gives Weibull, logistic
// gives log-logistic and normal gives log-normal.
// ---------------------------------------------------------------------------

struct AFTFamily : Family {
  enum Dist { WEIBULL, LOGLOGISTIC, LOGNORMAL };
  Dist dist;
  arma::vec event;
  double sigma;
  double sigma_hat;
  bool update_sigma;

  AFTFamily(const arma::vec& y_, const arma::vec& w_, Dist dist_,
            const arma::vec& event_, double sigma_hat_, bool update_sigma_)
    : Family(y_, w_, 1), dist(dist_), event(event_), sigma(sigma_hat_),
      sigma_hat(sigma_hat_), update_sigma(update_sigma_) {}

  static double loglik_one(Dist dist, double y, double eta, double delta,
                           double sigma) {
    double r = (y - eta) / sigma;
    double log_sigma = std::log(sigma);
    // Retained in full for the nuisance-parameter update, which varies sigma.
    switch (dist) {
    case WEIBULL:
      // Standard Gumbel minimum: density exp(r - exp(r)), survival exp(-exp(r))
      return delta * (r - log_sigma) - std::exp(r);
    case LOGLOGISTIC:
      return delta * (r - log_sigma) - (1.0 + delta) * std::log1p(std::exp(r));
    default:
      return delta * (R::dnorm4(r, 0.0, 1.0, 1) - log_sigma) +
        (1.0 - delta) * R::pnorm5(-r, 0.0, 1.0, 1, 1);
    }
  }

  double logdens_unit(int i, const double* eta) const override {
    return loglik_one(dist, y(i), eta[0], event(i), sigma);
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    double r = (y(i) - eta[0]) / sigma;
    double delta = event(i);
    switch (dist) {
    case WEIBULL:
      return (std::exp(r) - delta) / sigma;
    case LOGLOGISTIC:
      return ((1.0 + delta) * expit(r) - delta) / sigma;
    default: {
      if (delta > 0.0) {
        return r / sigma;
      }
      return inv_mills(r) / sigma;
    }
    }
  }

  double info_unit(int i, const double* eta, int h) const override {
    double r = (y(i) - eta[0]) / sigma;
    double delta = event(i);
    double s2 = sigma * sigma;
    switch (dist) {
    case WEIBULL:
      return std::exp(r) / s2;
    case LOGLOGISTIC: {
      double p = expit(r);
      return (1.0 + delta) * p * (1.0 - p) / s2;
    }
    default: {
      if (delta > 0.0) {
        return 1.0 / s2;
      }
      double lambda = inv_mills(r);
      return lambda * (lambda - r) / s2;
    }
    }
  }

  void update_aux(const arma::mat& eta) override {
    if (!update_sigma) {
      return;
    }
    const arma::rowvec& e = eta.row(0);
    auto logf = [this, &e](double log_sigma) {
      double s = std::exp(log_sigma);
      if (!(s > 0.0) || !std::isfinite(s)) {
        return R_NegInf;
      }
      // Half-Cauchy prior on sigma, with the Jacobian for the log scale.
      double out = Rf_dcauchy(s, 0.0, sigma_hat, 1) + log_sigma;
      for (int i = 0; i < N; i++) {
        out += w(i) * loglik_one(dist, y(i), e(i), event(i), s);
      }
      return out;
    };
    sigma = std::exp(slice_sampler(std::log(sigma), logf, 0.5, -20.0, 20.0));
  }

  std::vector<std::string> aux_names() const override {
    return std::vector<std::string>{"sigma"};
  }

  arma::vec aux_values() const override { return arma::vec{sigma}; }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem > 0) {
      sigma = values(0);
    }
  }
};

// ---------------------------------------------------------------------------
// Gaussian location-scale regression: one forest for the mean and a second for
// the log standard deviation, so the variance is an unrestricted function of
// the predictors.
// ---------------------------------------------------------------------------

struct LocationScaleFamily : Concrete<LocationScaleFamily> {

  // Quadratic in the mean but not in the log standard deviation, which is why
  // this is asked per predictor rather than per family.
  TargetForm target_form(int h) const override {
    return h == 0 ? TARGET_QUADRATIC : TARGET_GENERAL;
  }

  LocationScaleFamily(const arma::vec& y_, const arma::vec& w_)
    : Concrete<LocationScaleFamily>(y_, w_, 2) {}

  double logdens_unit(int i, const double* eta) const override {
    double r = (y(i) - eta[0]) * std::exp(-eta[1]);
    return -0.5 * LN_2PI - eta[1] - 0.5 * r * r;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    double r = (y(i) - eta[0]) * std::exp(-eta[1]);
    if (h == 0) {
      return r * std::exp(-eta[1]);
    }
    return r * r - 1.0;
  }

  double info_unit(int i, const double* eta, int h) const override {
    if (h == 0) {
      return std::exp(-2.0 * eta[1]);
    }
    // Expected information for the log scale, which is constant and avoids the
    // noise in the observed version.
    return 2.0;
  }
};

// ---------------------------------------------------------------------------
// Zero-inflated counts. Two additive predictors: the first is the log mean of
// the count component, the second the logit of the probability that an
// observation is a structural zero. Both are modeled nonparametrically, so the
// excess-zero mechanism is free to depend on the predictors.
// ---------------------------------------------------------------------------

struct ZeroInflatedFamily : Family {
  bool negbin;
  double theta;
  double prior_shape;
  double prior_rate;
  bool update_theta;

  ZeroInflatedFamily(const arma::vec& y_, const arma::vec& w_, bool negbin_,
                     double theta_, double prior_shape_, double prior_rate_,
                     bool update_theta_)
    : Family(y_, w_, 2), negbin(negbin_), theta(theta_),
      prior_shape(prior_shape_), prior_rate(prior_rate_),
      update_theta(update_theta_) {}

  // log P(count component yields zero).
  double log_p0(double eta_count, double th) const {
    double mu = std::exp(eta_count);
    if (!negbin) {
      return -mu;
    }
    return th * (std::log(th) - std::log(th + mu));
  }

  double count_logpmf(double y_i, double eta_count, double th) const {
    double mu = std::exp(eta_count);
    if (!negbin) {
      return y_i * eta_count - mu - R::lgammafn(y_i + 1.0);
    }
    double log_theta_mu = std::log(th + mu);
    return R::lgammafn(y_i + th) - R::lgammafn(th) - R::lgammafn(y_i + 1.0) +
      th * (std::log(th) - log_theta_mu) + y_i * (eta_count - log_theta_mu);
  }

  double loglik_one(double y_i, const double* eta, double th) const {
    if (y_i > 0.0) {
      return log1m_expit(eta[1]) + count_logpmf(y_i, eta[0], th);
    }
    // A zero arises either structurally or from the count component.
    return log_sum_exp(log_expit(eta[1]),
                       log1m_expit(eta[1]) + log_p0(eta[0], th));
  }

  // Only observations with a positive count reach the count log-density, and
  // for those the log-gamma terms depend on y and theta but not on eta.
  arma::vec compute_eta_free() const override {
    arma::vec out(N, arma::fill::zeros);
    double lg_theta = negbin ? R::lgammafn(theta) : 0.0;
    for (int i = 0; i < N; i++) {
      if (y(i) > 0.0) {
        out(i) = -R::lgammafn(y(i) + 1.0);
        if (negbin) {
          out(i) += R::lgammafn(y(i) + theta) - lg_theta;
        }
      }
    }
    return out;
  }

  double logdens_unit(int i, const double* eta) const override {
    double y_i = y(i);
    if (y_i <= 0.0) {
      return log_sum_exp(log_expit(eta[1]),
                         log1m_expit(eta[1]) + log_p0(eta[0], theta));
    }
    double mu = std::exp(eta[0]);
    if (!negbin) {
      return log1m_expit(eta[1]) + y_i * eta[0] - mu;
    }
    double log_theta_mu = std::log(theta + mu);
    return log1m_expit(eta[1]) +
      theta * (std::log(theta) - log_theta_mu) + y_i * (eta[0] - log_theta_mu);
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    double pi = expit(eta[1]);
    double mu = std::exp(eta[0]);

    if (y(i) > 0.0) {
      if (h == 1) {
        return -pi;
      }
      if (!negbin) {
        return y(i) - mu;
      }
      return theta * (y(i) - mu) / (theta + mu);
    }

    double log_structural = log_expit(eta[1]);
    double log_count = log1m_expit(eta[1]) + log_p0(eta[0], theta);
    double log_total = log_sum_exp(log_structural, log_count);
    double p0 = std::exp(log_p0(eta[0], theta));

    if (h == 1) {
      return pi * (1.0 - pi) * (1.0 - p0) / std::exp(log_total);
    }

    // Only the count term depends on the first predictor, so its share of the
    // zero probability weights the derivative of log P(count yields zero).
    double share = std::exp(log_count - log_total);
    double d_log_p0 = negbin ? -theta * mu / (theta + mu) : -mu;
    return share * d_log_p0;
  }

  // Complete-data expected information, which is positive by construction. The
  // observed information of a mixture can be negative, and only the proposal
  // depends on this, so a well-behaved approximation is preferable to an exact
  // quantity that sometimes has the wrong sign.
  double info_unit(int i, const double* eta, int h) const override {
    double pi = expit(eta[1]);

    if (h == 1) {
      return pi * (1.0 - pi);
    }

    double mu = std::exp(eta[0]);
    double count_info = negbin ? mu * theta / (mu + theta) : mu;
    return (1.0 - pi) * count_info;
  }

  void update_aux(const arma::mat& eta) override {
    if (!negbin || !update_theta) {
      return;
    }
    auto logf = [this, &eta](double log_theta) {
      double th = std::exp(log_theta);
      if (!(th > 0.0) || !std::isfinite(th)) {
        return R_NegInf;
      }
      double out = prior_shape * log_theta - prior_rate * th;
      for (int i = 0; i < N; i++) {
        out += w(i) * loglik_one(y(i), eta.colptr(i), th);
      }
      return out;
    };
    theta = std::exp(slice_sampler(std::log(theta), logf, 1.0, -20.0, 20.0));
    refresh_eta_free();
  }

  std::vector<std::string> aux_names() const override {
    if (!negbin) {
      return std::vector<std::string>();
    }
    return std::vector<std::string>{"theta"};
  }

  arma::vec aux_values() const override {
    if (!negbin) {
      return arma::vec();
    }
    return arma::vec{theta};
  }

  void set_aux(const arma::vec& values) override {
    if (negbin && values.n_elem > 0) {
      theta = values(0);
      refresh_eta_free();
    }
  }
};

// ---------------------------------------------------------------------------
// Ordered beta regression (Kubinec 2023). The response lies on the closed unit
// interval, with point masses at both ends. A single additive predictor drives
// both the probability of landing on an endpoint, through a pair of cutpoints
// as in an ordinal model, and the mean of the beta density in between. Because
// the predictor also enters the beta mean it is identified, so unlike the
// ordinal family neither cutpoint has to be pinned.
// ---------------------------------------------------------------------------

struct OrdBetaFamily : Family {
  double cut1;
  double cut2;
  double phi;
  double prior_shape;
  double prior_rate;
  bool update_phi;

  // log(y) and log(1 - y) appear in the beta density multiplied by shapes that
  // do move with eta, so they cannot go into the eta-free part, but they are
  // still fixed per observation and worth computing once.
  arma::vec log_y;
  arma::vec log1m_y;

  OrdBetaFamily(const arma::vec& y_, const arma::vec& w_, double cut1_,
                double cut2_, double phi_, double prior_shape_,
                double prior_rate_, bool update_phi_)
    : Family(y_, w_, 1), cut1(cut1_), cut2(cut2_), phi(phi_),
      prior_shape(prior_shape_), prior_rate(prior_rate_),
      update_phi(update_phi_) {
    log_y.set_size(N);
    log1m_y.set_size(N);
    for (int i = 0; i < N; i++) {
      bool interior = y_(i) > 0.0 && y_(i) < 1.0;
      log_y(i) = interior ? std::log(y_(i)) : 0.0;
      log1m_y(i) = interior ? std::log1p(-y_(i)) : 0.0;
    }
  }

  // lbeta(a, b) = lgamma(a) + lgamma(b) - lgamma(phi), and the last term is
  // free of eta because the two shapes always sum to phi. The unit constants of
  // the beta density go the same way.
  arma::vec compute_eta_free() const override {
    arma::vec out(N, arma::fill::zeros);
    double lg_phi = R::lgammafn(phi);
    for (int i = 0; i < N; i++) {
      if (y(i) > 0.0 && y(i) < 1.0) {
        out(i) = lg_phi - log_y(i) - log1m_y(i);
      }
    }
    return out;
  }

  // The hot path: the two shape-dependent log-gamma terms and the middle
  // interval, with everything else already accounted for above.
  double logdens_unit(int i, const double* eta) const override {
    double e = eta[0];
    double y_i = y(i);

    if (y_i <= 0.0) {
      return log1m_expit(e - cut1);
    }
    if (y_i >= 1.0) {
      return log_expit(e - cut2);
    }

    double mu = expit(e);
    double a = mu * phi;
    double b = phi - a;

    if (!(a > 0.0) || !(b > 0.0)) {
      return R_NegInf;
    }

    return log_diff_logistic(e - cut2, e - cut1) + a * log_y(i) +
      b * log1m_y(i) - R::lgammafn(a) - R::lgammafn(b);
  }

  // Complete log density, used by the nuisance-parameter updates, which vary
  // the cutpoints and the precision and so cannot use the cached eta-free part.
  static double loglik_one(double y_i, double eta, double c1, double c2,
                           double phi_) {
    if (y_i <= 0.0) {
      return log1m_expit(eta - c1);
    }
    if (y_i >= 1.0) {
      return log_expit(eta - c2);
    }
    // Probability of the continuous middle, times the beta density on it.
    double log_middle = log_diff_logistic(eta - c2, eta - c1);
    double mu = expit(eta);
    double a = mu * phi_;
    double b = phi_ - a;

    if (!(a > 0.0) || !(b > 0.0)) {
      return R_NegInf;
    }

    double log_beta = (a - 1.0) * std::log(y_i) +
      (b - 1.0) * std::log1p(-y_i) - Rf_lbeta(a, b);
    return log_middle + log_beta;
  }

  // The endpoint cases are ordinary logistic terms. For the interior, the
  // derivative splits into the probability of the middle interval and the beta
  // density. Differentiating the beta term, the two digamma contributions from
  // log Beta(a, b) keep only their difference, because the shapes move in
  // opposite directions: da/deta = -db/deta. The information uses the expected
  // beta information s^2 (psi'(a) + psi'(b)) rather than the observed version,
  // whose extra term has mean zero but can turn the total negative.
  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double e = eta[0];
    double y_i = y(i);

    if (y_i <= 0.0) {
      double a = expit(e - cut1);
      *d1 = -a;
      *d2 = a * (1.0 - a);
      return;
    }

    if (y_i >= 1.0) {
      double b = expit(e - cut2);
      *d1 = 1.0 - b;
      *d2 = b * (1.0 - b);
      return;
    }

    double upper = expit(e - cut1);
    double lower = expit(e - cut2);
    double span = upper - lower;

    if (!(span > 1e-300) || !std::isfinite(span)) {
      score_info_numeric(i, eta, h, d1, d2);
      return;
    }

    double f_up = upper * (1.0 - upper);
    double f_lo = lower * (1.0 - lower);
    double diff = f_up - f_lo;

    double d1_span = diff / span;
    double d2_span = ((f_lo * (1.0 - 2.0 * lower) -
                       f_up * (1.0 - 2.0 * upper)) * span + diff * diff) /
      (span * span);

    double mu = expit(e);
    double a = mu * phi;
    double b = phi - a;

    if (!(a > 0.0) || !(b > 0.0)) {
      score_info_numeric(i, eta, h, d1, d2);
      return;
    }

    double slope = phi * mu * (1.0 - mu);
    double centered = std::log(y_i) - std::log1p(-y_i) - R::digamma(a) +
      R::digamma(b);

    *d1 = d1_span + slope * centered;
    *d2 = d2_span + slope * slope * (R::trigamma(a) + R::trigamma(b));
  }

  void update_aux(const arma::mat& eta) override {
    const arma::rowvec e = eta.row(0);

    auto total = [this, &e](double c1, double c2, double p) {
      double out = 0.0;
      for (int i = 0; i < N; i++) {
        out += w(i) * loglik_one(y(i), e(i), c1, c2, p);
      }
      return out;
    };

    double c2_now = cut2;
    cut1 = slice_sampler(cut1, [&](double v) {
      if (!(v < c2_now)) {
        return R_NegInf;
      }
      return total(v, c2_now, phi);
    }, 0.5, -30.0, c2_now);

    double c1_now = cut1;
    cut2 = slice_sampler(cut2, [&](double v) {
      if (!(v > c1_now)) {
        return R_NegInf;
      }
      return total(c1_now, v, phi);
    }, 0.5, c1_now, 30.0);

    if (update_phi) {
      double c1 = cut1;
      double c2 = cut2;
      phi = std::exp(slice_sampler(std::log(phi), [&](double log_phi) {
        double p = std::exp(log_phi);
        if (!(p > 0.0) || !std::isfinite(p)) {
          return R_NegInf;
        }
        return prior_shape * log_phi - prior_rate * p + total(c1, c2, p);
      }, 0.5, -10.0, 15.0));
    }

    refresh_eta_free();
  }

  std::vector<std::string> aux_names() const override {
    return std::vector<std::string>{"cut1", "cut2", "phi"};
  }

  arma::vec aux_values() const override {
    return arma::vec{cut1, cut2, phi};
  }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem == 3) {
      cut1 = values(0);
      cut2 = values(1);
      phi = values(2);
      refresh_eta_free();
    }
  }
};

// ---------------------------------------------------------------------------
// Polya-Gamma augmented families.
//
// All three of these rest on the same identity. A likelihood proportional to
// exp(kappa * psi) / (1 + exp(psi))^b becomes, after introducing
// omega ~ PG(b, psi), proportional to exp(kappa * psi - omega * psi^2 / 2),
// which is Gaussian in psi. The only thing that differs between them is what
// psi, kappa and b are:
//
//   binomial logit    psi = eta,            kappa = s - w/2,     b = w
//   negative binomial psi = eta - log theta, kappa = (y - theta)/2, b = y + theta
//   multinomial       psi = eta_j - log C_j, kappa = y_j - n/2,   b = n
//
// where s is the success count, w the number of trials, and C_j the sum of
// exp(eta) over the categories other than j. The last of these is why the
// augmentation is refreshed per forest rather than per sweep: C_j depends on the
// other categories' predictors, which move during the sweep.
//
// Every one of these families is quadratic in the predictor, which is the whole
// point: the Laplace approximation stops being an approximation.
//
// The base class multiplies logdens_unit() by the prior weight, so these
// families carry unit weights and fold the real ones into kappa, omega and the
// reported likelihood. There is no way to split a Polya-Gamma draw into
// per-unit-weight pieces, since PG(w, c) is not w times PG(1, c).
// ---------------------------------------------------------------------------

namespace {

arma::vec unit_weights(arma::uword n) {
  return arma::vec(n, arma::fill::ones);
}

} // namespace

struct LogitAugmentedFamily : Concrete<LogitAugmentedFamily> {
  arma::vec trials;     // the prior weights, which are the binomial denominators
  arma::vec successes;
  arma::vec kappa;
  arma::vec omega;

  LogitAugmentedFamily(const arma::vec& y_, const arma::vec& w_)
    : Concrete<LogitAugmentedFamily>(y_, unit_weights(y_.n_elem), 1), trials(w_) {
    successes = w_ % y_;
    kappa = successes - 0.5 * w_;
    omega = arma::vec(N, arma::fill::ones);
  }

  static bool applies(const arma::vec& w) {
    return arma::all(w > 0.0);
  }

  TargetForm target_form(int h) const override {
    return TARGET_QUADRATIC;
  }

  double logdens_unit(int i, const double* eta) const override {
    double e = eta[0];
    return kappa(i) * e - 0.5 * omega(i) * e * e;
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    *d1 = kappa(i) - omega(i) * eta[0];
    *d2 = omega(i);
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return kappa(i) - omega(i) * eta[0];
  }

  double info_unit(int i, const double* eta, int h) const override {
    return omega(i);
  }

  void before_forest(int h, const arma::mat& eta) override {
    const arma::rowvec& e = eta.row(0);

    for (int i = 0; i < N; i++) {
      omega(i) = rpg(trials(i), e(i));
    }
  }

  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      out += successes(i) * log_expit(e(i)) +
        (trials(i) - successes(i)) * log1m_expit(e(i)) +
        R::lgammafn(trials(i) + 1.0) - R::lgammafn(successes(i) + 1.0) -
        R::lgammafn(trials(i) - successes(i) + 1.0);
    }

    return out;
  }
};

// The negative binomial as a Poisson-gamma mixture, which is the augmentation
// Hill et al. (2020, sec. 3.1.5) point at and a better bargain than the
// Polya-Gamma one this replaced.
//
// A negative binomial count with mean mu and size theta is a Poisson count whose
// rate is itself drawn from a gamma: y ~ Poisson(lambda) with
// lambda ~ Gamma(theta, rate theta / mu). Introduce that lambda and the
// predictor's whole contribution to the log joint is
//
//     -theta * eta - lambda * theta * exp(-eta),
//
// which is TARGET_EXP_DOWN -- the same shape as the gamma family, with lambda in
// place of the response and theta in place of the shape. So the leaf-level work
// collapses to three numbers from one pass, exactly as it does for a Poisson or
// a gamma response.
//
// The cost of the augmentation is one gamma draw per observation per sweep,
// against the two hundred a Polya-Gamma draw needed here: PG(y + theta, .) has a
// non-integer parameter, so it could not use Devroye's exact method and fell back
// on a truncated series, which cost more than the likelihood it was replacing.
// The conditional needed here is exact and takes one line.
//
// theta is drawn from its conditional with lambda integrated out -- that is, from
// the true negative binomial likelihood, exactly as the direct family does -- and
// lambda is then redrawn given the new theta. Updating a parameter from its
// collapsed conditional and then the latent variable it was collapsed over is a
// valid partially collapsed Gibbs step in that order (Van Dyk and Park 2008).
struct NegBinAugmentedFamily : Concrete<NegBinAugmentedFamily> {
  arma::vec count;
  arma::vec rate;     // lambda, the Poisson rate, redrawn every sweep
  double theta;
  double log_theta;
  double prior_shape;
  double prior_rate;
  bool update_theta;

  NegBinAugmentedFamily(const arma::vec& y_, const arma::vec& w_, double theta_,
                        double prior_shape_, double prior_rate_,
                        bool update_theta_)
    : Concrete<NegBinAugmentedFamily>(y_, unit_weights(y_.n_elem), 1), count(y_), theta(theta_),
      prior_shape(prior_shape_), prior_rate(prior_rate_),
      update_theta(update_theta_) {
    log_theta = std::log(theta);
    // A deterministic start, replaced by a proper draw before the first forest
    // moves. The counts themselves are the natural guess at their own rates.
    rate = arma::clamp(y_, 0.1, arma::datum::inf);
  }

  static bool applies(const arma::vec& w) {
    return arma::all(w == 1.0);
  }

  TargetForm target_form(int h) const override {
    return TARGET_EXP_DOWN;
  }

  double logdens_unit(int i, const double* eta) const override {
    return -theta * (eta[0] + rate(i) * std::exp(-eta[0]));
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double scaled = theta * rate(i) * std::exp(-eta[0]);
    *d1 = scaled - theta;
    *d2 = scaled;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return theta * (rate(i) * std::exp(-eta[0]) - 1.0);
  }

  double info_unit(int i, const double* eta, int h) const override {
    return theta * rate(i) * std::exp(-eta[0]);
  }

  // lambda | y, eta, theta is Gamma(y + theta, 1 + theta / mu).
  void before_forest(int h, const arma::mat& eta) override {
    const arma::rowvec& e = eta.row(0);

    for (int i = 0; i < N; i++) {
      double scale = 1.0 / (1.0 + theta * std::exp(-e(i)));
      rate(i) = Rf_rgamma(count(i) + theta, scale);
    }
  }

  void update_aux(const arma::mat& eta) override {
    if (!update_theta) {
      return;
    }

    const arma::rowvec& e = eta.row(0);
    auto logf = [this, &e](double candidate) {
      double th = std::exp(candidate);

      if (!(th > 0.0) || !std::isfinite(th)) {
        return R_NegInf;
      }

      double out = prior_shape * candidate - prior_rate * th;

      for (int i = 0; i < N; i++) {
        out += NegBinFamily::loglik_one(count(i), e(i), th);
      }

      return out;
    };

    theta = std::exp(slice_sampler(log_theta, logf, 1.0, -20.0, 20.0));
    log_theta = std::log(theta);
  }

  std::vector<std::string> aux_names() const override {
    return std::vector<std::string>{"theta"};
  }

  arma::vec aux_values() const override { return arma::vec{theta}; }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem > 0) {
      theta = values(0);
      log_theta = std::log(theta);
    }
  }

  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      out += NegBinFamily::loglik_one(count(i), e(i), theta);
    }

    return out;
  }
};

// Multinomial logistic, one category at a time. Conditional on the other
// categories the likelihood of category j is exactly binomial-logistic in
// eta_j - log C_j, so the same augmentation applies with one Polya-Gamma draw
// per observation per category.
struct MultinomAugmentedFamily : Concrete<MultinomAugmentedFamily> {
  arma::vec category;
  int num_cat;
  bool symmetric;
  arma::vec offset;   // log C_j for the forest currently being updated
  arma::vec kappa;
  arma::vec omega;
  int current;

  MultinomAugmentedFamily(const arma::vec& y_, const arma::vec& w_, int num_cat_,
                          bool symmetric_)
    : Concrete<MultinomAugmentedFamily>(y_, unit_weights(y_.n_elem), symmetric_ ? num_cat_ : num_cat_ - 1),
      category(y_), num_cat(num_cat_), symmetric(symmetric_), current(0) {
    offset = arma::vec(N, arma::fill::zeros);
    kappa = arma::vec(N, arma::fill::zeros);
    omega = arma::vec(N, arma::fill::ones);
  }

  static bool applies(const arma::vec& w) {
    return arma::all(w == 1.0);
  }

  TargetForm target_form(int h) const override {
    return TARGET_QUADRATIC;
  }

  // Which predictor belongs to the observed category, or -1 for the reference.
  int own(int i) const {
    int cat = static_cast<int>(category(i));
    return symmetric ? cat : cat - 1;
  }

  double logdens_unit(int i, const double* eta) const override {
    double psi = eta[current] - offset(i);
    return kappa(i) * psi - 0.5 * omega(i) * psi * psi;
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double psi = eta[h] - offset(i);
    *d1 = kappa(i) - omega(i) * psi;
    *d2 = omega(i);
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return kappa(i) - omega(i) * (eta[h] - offset(i));
  }

  double info_unit(int i, const double* eta, int h) const override {
    return omega(i);
  }

  void before_forest(int h, const arma::mat& eta) override {
    current = h;

    for (int i = 0; i < N; i++) {
      // log of the sum of exp(eta) over the other categories, including the
      // reference category's implicit zero when the coding has one.
      double m = symmetric ? R_NegInf : 0.0;

      for (int l = 0; l < H; l++) {
        if (l != h && eta(l, i) > m) {
          m = eta(l, i);
        }
      }

      double total = 0.0;

      if (!symmetric) {
        total += std::exp(-m);
      }

      for (int l = 0; l < H; l++) {
        if (l != h) {
          total += std::exp(eta(l, i) - m);
        }
      }

      offset(i) = m + std::log(total);
      kappa(i) = (own(i) == h ? 1.0 : 0.0) - 0.5;
      omega(i) = rpg(1.0, eta(h, i) - offset(i));
    }
  }

  double reported_loglik(const arma::mat& eta) const override {
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      double m = symmetric ? eta(0, i) : 0.0;

      for (int l = 0; l < H; l++) {
        if (eta(l, i) > m) {
          m = eta(l, i);
        }
      }

      double total = symmetric ? 0.0 : std::exp(-m);

      for (int l = 0; l < H; l++) {
        total += std::exp(eta(l, i) - m);
      }

      int j = own(i);
      out += (j >= 0 ? eta(j, i) : 0.0) - m - std::log(total);
    }

    return out;
  }
};

// ---------------------------------------------------------------------------
// Families whose likelihood or link is supplied from R.
//
// Both call back into the interpreter, where the cost of a call swamps the
// arithmetic inside it, so both take a whole leaf at a time: `wants_block()`
// is true and the block methods do the work. That turns one call per
// observation per Fisher-scoring step into one call per leaf per step.
// ---------------------------------------------------------------------------

namespace {

// Evaluate an R function on a block of predictors and check what comes back.
// The result has to be a plain numeric vector of the length that went in;
// anything else is a bug in the caller's function and is worth saying so about
// rather than reading past the end of a vector.
arma::vec call_r(const Rcpp::Function& f, const double* x, int n,
                 const char* what) {
  Rcpp::NumericVector arg(x, x + n);
  Rcpp::RObject value = f(arg);

  if (!Rf_isReal(value) && !Rf_isInteger(value)) {
    Rcpp::stop("the %s supplied to genbart() must return a numeric vector.",
               what);
  }

  arma::vec out = Rcpp::as<arma::vec>(value);

  if (static_cast<int>(out.n_elem) != n) {
    Rcpp::stop("the %s supplied to genbart() returned %d values for %d "
               "observations.", what, static_cast<int>(out.n_elem), n);
  }

  return out;
}

} // namespace

// A family with a link supplied from R, wrapped around one of the compiled
// families. The additive predictor is mapped to the scale the wrapped family
// works on -- the mean for a Gaussian response, the log mean for the counts,
// the log odds for a binomial -- by a function built on the R side from the
// caller's link, and the chain rule carries the derivatives back.
//
// Only families with a single additive predictor and a conventional link are
// wrapped, which is what makes the composition well defined: theta = t(eta)
// with t scalar and monotone.
//
// The information is scaled by t'(eta)^2 and the term in t''(eta) is dropped.
// That term is the score times t'', and the score has expectation zero, so what
// remains is exactly the expected Fisher information of the composite whenever
// the wrapped family reports the expected information. It is also guaranteed
// non-negative, where the full second derivative is not. Since these numbers
// only ever build a proposal, and the exact log density is what the acceptance
// ratio uses, dropping the term costs a little efficiency and nothing else.
struct LinkedFamily : Family {
  std::unique_ptr<Family> inner;
  Rcpp::Function theta;
  Rcpp::Function dtheta;
  bool has_dtheta;

  LinkedFamily(Family* inner_, const Rcpp::Function& theta_,
               const Rcpp::RObject& dtheta_)
    : Family(inner_->y, inner_->w, 1), inner(inner_), theta(theta_),
      dtheta(Rf_isFunction(dtheta_) ? Rcpp::Function(dtheta_) : theta_),
      has_dtheta(Rf_isFunction(dtheta_)) {}

  bool wants_block() const override { return true; }

  // Central difference when the caller's link brought no derivative with it.
  // Deterministic in eta, so the proposal stays reversible.
  arma::vec slope(const double* x, int n) const {
    if (has_dtheta) {
      return call_r(dtheta, x, n, "link derivative");
    }

    const double step = 1e-5;
    arma::vec up(n);
    arma::vec down(n);

    for (int k = 0; k < n; k++) {
      up(k) = x[k] + step;
      down(k) = x[k] - step;
    }

    arma::vec a = call_r(theta, up.memptr(), n, "link");
    arma::vec b = call_r(theta, down.memptr(), n, "link");
    return (a - b) / (2.0 * step);
  }

  void logdens_block(const int* idx, int n, const double* block,
                     double* out) const override {
    arma::vec th = call_r(theta, block, n, "link");
    inner->logdens_block(idx, n, th.memptr(), out);
  }

  void score_info_block(const int* idx, int n, const double* block, int,
                        double* d1, double* d2) const override {
    arma::vec th = call_r(theta, block, n, "link");
    inner->score_info_block(idx, n, th.memptr(), 0, d1, d2);
    arma::vec tp = slope(block, n);

    for (int k = 0; k < n; k++) {
      d1[k] *= tp(k);
      d2[k] *= tp(k) * tp(k);
    }
  }

  // The single-observation route, for the places that report a density rather
  // than drive the sampler. Correct but slow, and never on the hot path,
  // because wants_block() sends the sampler through the block methods. It
  // applies the same chain rule, so the two paths report the same numbers.
  double logdens_unit(int i, const double* eta) const override {
    arma::vec th = call_r(theta, eta, 1, "link");
    return inner->logdens_unit(i, th.memptr());
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    arma::vec th = call_r(theta, eta, 1, "link");
    inner->score_info_unit(i, th.memptr(), 0, d1, d2);
    arma::vec tp = slope(eta, 1);
    *d1 *= tp(0);
    *d2 *= tp(0) * tp(0);
  }

  double log_norm_const(int i) const override {
    return inner->log_norm_const(i);
  }

  arma::vec compute_eta_free() const override {
    return inner->eta_free_part();
  }

  double logdens_extra_total() const override {
    return inner->logdens_extra_total();
  }

  // The nuisance parameters belong to the wrapped family and are drawn on its
  // scale, so the predictors are transformed before they are handed over.
  void update_aux(const arma::mat& eta) override {
    arma::vec th = call_r(theta, eta.memptr(), static_cast<int>(eta.n_elem),
                          "link");
    inner->update_aux(arma::mat(th.memptr(), 1, eta.n_cols));
    refresh_eta_free();
  }

  std::vector<std::string> aux_names() const override {
    return inner->aux_names();
  }

  arma::vec aux_values() const override { return inner->aux_values(); }

  void set_aux(const arma::vec& values) override {
    inner->set_aux(values);
    refresh_eta_free();
  }
};

// A family whose log density is an R function. The caller supplies only the log
// density of one unit of prior weight, as a function of the response and the
// additive predictors; the derivatives are central differences of it unless a
// second function supplies them.
//
// The whole leaf goes to R in one call. The differences reuse the middle
// evaluation, so a Fisher-scoring step costs three calls rather than five, and
// the step is chosen for the second difference rather than the first: the
// gradient it produces is a little coarse, but it is a deterministic function
// of the state, which is what reversibility needs, and only the proposal
// depends on it.
struct RFamily : Family {
  Rcpp::Function dens;
  Rcpp::Function derivs;
  bool has_derivs;
  std::string label;

  RFamily(const arma::vec& y_, const arma::vec& w_, int H_,
          const Rcpp::Function& dens_, const Rcpp::RObject& derivs_,
          const std::string& label_)
    : Family(y_, w_, H_), dens(dens_),
      derivs(Rf_isFunction(derivs_) ? Rcpp::Function(derivs_) : dens_),
      has_derivs(Rf_isFunction(derivs_)), label(label_) {}

  bool wants_block() const override { return true; }

  // The response and the predictors of one block, in the shape the caller's
  // function expects: y a vector of length n, eta an n by H matrix.
  void unpack(const int* idx, int n, const double* block,
              Rcpp::NumericVector& y_out, Rcpp::NumericMatrix& eta_out) const {
    for (int k = 0; k < n; k++) {
      y_out[k] = y(idx[k]);

      for (int j = 0; j < H; j++) {
        eta_out(k, j) = block[static_cast<std::size_t>(k) * H + j];
      }
    }
  }

  arma::vec evaluate(const Rcpp::NumericVector& y_in,
                     const Rcpp::NumericMatrix& eta_in, int n) const {
    Rcpp::RObject value = dens(y_in, eta_in);

    if (!Rf_isReal(value) && !Rf_isInteger(value)) {
      Rcpp::stop("the log density supplied to custom_family() must return a "
                 "numeric vector.");
    }

    arma::vec out = Rcpp::as<arma::vec>(value);

    if (static_cast<int>(out.n_elem) != n) {
      Rcpp::stop("the log density supplied to custom_family() returned %d "
                 "values for %d observations.",
                 static_cast<int>(out.n_elem), n);
    }

    return out;
  }

  void logdens_block(const int* idx, int n, const double* block,
                     double* out) const override {
    Rcpp::NumericVector yv(n);
    Rcpp::NumericMatrix ev(n, H);
    unpack(idx, n, block, yv, ev);
    arma::vec value = evaluate(yv, ev, n);

    for (int k = 0; k < n; k++) {
      out[k] = w(idx[k]) * value(k);
    }
  }

  void score_info_block(const int* idx, int n, const double* block, int h,
                        double* d1, double* d2) const override {
    Rcpp::NumericVector yv(n);
    Rcpp::NumericMatrix ev(n, H);
    unpack(idx, n, block, yv, ev);

    if (has_derivs) {
      Rcpp::RObject value = derivs(yv, ev, h + 1);
      Rcpp::List parts(value);
      arma::vec score = Rcpp::as<arma::vec>(parts["score"]);
      arma::vec info = Rcpp::as<arma::vec>(parts["info"]);

      if (static_cast<int>(score.n_elem) != n ||
          static_cast<int>(info.n_elem) != n) {
        Rcpp::stop("the derivatives supplied to custom_family() must return "
                   "one score and one information value per observation.");
      }

      for (int k = 0; k < n; k++) {
        d1[k] = w(idx[k]) * score(k);
        d2[k] = clamp_info(w(idx[k]) * info(k));
      }

      return;
    }

    const double step = 1e-4;
    arma::vec mid = evaluate(yv, ev, n);

    for (int k = 0; k < n; k++) {
      ev(k, h) += step;
    }

    arma::vec up = evaluate(yv, ev, n);

    for (int k = 0; k < n; k++) {
      ev(k, h) -= 2.0 * step;
    }

    arma::vec down = evaluate(yv, ev, n);

    for (int k = 0; k < n; k++) {
      d1[k] = w(idx[k]) * 0.5 * (up(k) - down(k)) / step;
      d2[k] = clamp_info(-w(idx[k]) *
        (up(k) - 2.0 * mid(k) + down(k)) / (step * step));
    }
  }

  double logdens_unit(int i, const double* eta) const override {
    Rcpp::NumericVector yv(1);
    Rcpp::NumericMatrix ev(1, H);
    yv[0] = y(i);

    for (int j = 0; j < H; j++) {
      ev(0, j) = eta[j];
    }

    return evaluate(yv, ev, 1)(0);
  }

  // Supplied derivatives are used by the per-observation route as well, so that
  // the two routes report the same numbers rather than one of them quietly
  // falling back on differences.
  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    if (!has_derivs) {
      score_info_by_difference(i, eta, h, d1, d2);
      return;
    }

    int one = i;
    double a;
    double b;
    score_info_block(&one, 1, eta, h, &a, &b);
    // The block route applies the prior weight; the *_unit contract does not.
    *d1 = a / w(i);
    *d2 = b / w(i);
  }
};

// ---------------------------------------------------------------------------

namespace {

// The eta-free part is family-specific, so it cannot be filled from the base
// constructor; do it once here, after construction.
Family* finish(Family* family) {
  family->refresh_eta_free();
  return family;
}

} // namespace

namespace {

// Wrap the family in the caller's link, when one was supplied that the engine
// does not carry natively.
Family* with_link(Family* family, const List& opts) {
  if (!opts.containsElementNamed("link_theta")) {
    return family;
  }

  // The composition maps one predictor to one, which is the only case where it
  // is defined. The R side only offers the link for families with a single
  // predictor, so reaching this is a bug rather than a user error.
  if (family->H != 1) {
    Rcpp::stop("a link supplied from R needs a family with one additive "
               "predictor; this one has %d.", family->H);
  }

  RObject theta = as<RObject>(opts["link_theta"]);

  if (!Rf_isFunction(theta)) {
    return family;
  }

  RObject dtheta = R_NilValue;

  if (opts.containsElementNamed("link_dtheta")) {
    dtheta = as<RObject>(opts["link_dtheta"]);
  }

  return finish(new LinkedFamily(family, as<Function>(theta), dtheta));
}

} // namespace

Family* make_family(const std::string& name, const std::string& link,
                    const arma::vec& y, const arma::vec& w,
                    const List& opts) {

  if (name == "custom") {
    return finish(new RFamily(y, w, as<int>(opts["num_predictors"]),
                              as<Function>(opts["logdens"]),
                              as<RObject>(opts["derivatives"]),
                              as<std::string>(opts["name"])));
  }

  if (name == "gaussian") {
    return with_link(finish(new GaussianFamily(y, w,
                                              as<double>(opts["sigma_hat"]))),
                     opts);
  }

  if (name == "binomial") {
    BinomialFamily::Link l = BinomialFamily::LOGIT;
    if (link == "probit") {
      l = BinomialFamily::PROBIT;
    }
    else if (link == "cloglog") {
      l = BinomialFamily::CLOGLOG;
    }
    return with_link(finish(new BinomialFamily(y, w, l)), opts);
  }

  if (name == "poisson") {
    return with_link(finish(new PoissonFamily(y, w)), opts);
  }

  if (name == "negbin") {
    return with_link(finish(new NegBinFamily(y, w, as<double>(opts["theta"]),
                            as<double>(opts["theta_prior_shape"]),
                            as<double>(opts["theta_prior_rate"]),
                            as<bool>(opts["update_theta"]))), opts);
  }

  if (name == "Gamma") {
    return with_link(finish(new GammaFamily(y, w, as<double>(opts["shape"]),
                           as<double>(opts["shape_prior_shape"]),
                           as<double>(opts["shape_prior_rate"]),
                           as<bool>(opts["update_shape"]))), opts);
  }

  if (name == "ordinal") {
    int l = ORD_LOGIT;

    if (link == "probit") {
      l = ORD_PROBIT;
    }
    else if (link == "cloglog") {
      l = ORD_CLOGLOG;
    }

    return finish(new OrdinalFamily(y, w, l, as<int>(opts["num_cat"]),
                             as<arma::vec>(opts["cuts"]),
                             as<bool>(opts["update_cuts"])));
  }

  if (name == "multinomial") {
    return finish(new MultinomFamily(y, w, as<int>(opts["num_cat"]),
                                     as<bool>(opts["symmetric"])));
  }

  if (name == "aft") {
    AFTFamily::Dist d = AFTFamily::WEIBULL;
    if (link == "loglogistic") {
      d = AFTFamily::LOGLOGISTIC;
    }
    else if (link == "lognormal") {
      d = AFTFamily::LOGNORMAL;
    }
    return finish(new AFTFamily(y, w, d, as<arma::vec>(opts["event"]),
                         as<double>(opts["sigma_hat"]),
                         as<bool>(opts["update_sigma"])));
  }

  if (name == "location_scale") {
    return finish(new LocationScaleFamily(y, w));
  }

  if (name == "zip" || name == "zinb") {
    bool nb = name == "zinb";
    return finish(new ZeroInflatedFamily(y, w, nb,
                                  nb ? as<double>(opts["theta"]) : 1.0,
                                  nb ? as<double>(opts["theta_prior_shape"]) : 0.0,
                                  nb ? as<double>(opts["theta_prior_rate"]) : 0.0,
                                  nb ? as<bool>(opts["update_theta"]) : false));
  }

  if (name == "ordbeta") {
    return finish(new OrdBetaFamily(y, w, as<double>(opts["cut1"]),
                             as<double>(opts["cut2"]), as<double>(opts["phi"]),
                             as<double>(opts["phi_prior_shape"]),
                             as<double>(opts["phi_prior_rate"]),
                             as<bool>(opts["update_phi"])));
  }

  stop("Unsupported family '%s'.", name);
}

Family* augmented_family(const std::string& name, const std::string& link,
                         const arma::vec& y, const arma::vec& w,
                         const List& opts,
                         const std::vector<std::string>& enabled) {
  bool wanted = false;

  for (std::size_t k = 0; k < enabled.size(); k++) {
    if (enabled[k] == name) {
      wanted = true;
    }
  }

  if (!wanted) {
    return nullptr;
  }

  if (name == "binomial" && link == "probit" &&
      ProbitAugmentedFamily::applies(y, w)) {
    return finish(new ProbitAugmentedFamily(y, w));
  }

  if (name == "binomial" && link == "logit" &&
      LogitAugmentedFamily::applies(w)) {
    return finish(new LogitAugmentedFamily(y, w));
  }

  if (name == "negbin" && NegBinAugmentedFamily::applies(w)) {
    return finish(new NegBinAugmentedFamily(y, w, as<double>(opts["theta"]),
                                            as<double>(opts["theta_prior_shape"]),
                                            as<double>(opts["theta_prior_rate"]),
                                            as<bool>(opts["update_theta"])));
  }

  if (name == "ordinal" && link == "cloglog" &&
      OrdinalCloglogAugmentedFamily::applies(w)) {
    return finish(new OrdinalCloglogAugmentedFamily(
      y, w, as<int>(opts["num_cat"]), as<arma::vec>(opts["cuts"]),
      as<bool>(opts["update_cuts"])));
  }

  if (name == "ordinal" && link == "logit" &&
      OrdinalLogitAugmentedFamily::applies(w)) {
    return finish(new OrdinalLogitAugmentedFamily(
      y, w, as<int>(opts["num_cat"]), as<arma::vec>(opts["cuts"]),
      as<bool>(opts["update_cuts"])));
  }

  if (name == "ordinal" && link == "probit" &&
      OrdinalProbitAugmentedFamily::applies(w)) {
    return finish(new OrdinalProbitAugmentedFamily(
      y, w, as<int>(opts["num_cat"]), as<arma::vec>(opts["cuts"]),
      as<bool>(opts["update_cuts"])));
  }

  if (name == "multinomial" && MultinomAugmentedFamily::applies(w)) {
    return finish(new MultinomAugmentedFamily(y, w, as<int>(opts["num_cat"]),
                                               as<bool>(opts["symmetric"])));
  }

  return nullptr;
}

} // namespace genbart
