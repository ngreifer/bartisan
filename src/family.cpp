#include "family.h"
#include "polyagamma.h"
#include "slice.h"

using namespace Rcpp;

namespace bartisan {

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
      // Compared rather than passed to std::max(), which takes its arguments
      // by const reference: binding OMEGA_MIN to one is an ODR-use, and a
      // `static constexpr` member is only implicitly inline from C++17 on. Under
      // an older standard the symbol has no definition, which the optimizer
      // hides by folding the constant and a -O0 build does not.
      double drawn = rpg(2.0, std::fabs(r));
      omega(i) = drawn < OMEGA_MIN ? OMEGA_MIN : drawn;

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
    // `lt` rather than `log_theta`, which is the name of the cached member this
    // lambda would otherwise shadow: an edit inside here that meant the cache
    // would silently get the argument instead.
    auto logf = [this, &e](double lt) {
      double th = std::exp(lt);
      if (!(th > 0.0) || !std::isfinite(th)) {
        return R_NegInf;
      }
      double out = prior_shape * lt - prior_rate * th;
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

struct AFTFamily : Concrete<AFTFamily> {
  enum Dist { WEIBULL, LOGLOGISTIC, LOGNORMAL };
  Dist dist;
  arma::vec event;
  double sigma;
  double sigma_hat;
  bool update_sigma;

  AFTFamily(const arma::vec& y_, const arma::vec& w_, Dist dist_,
            const arma::vec& event_, double sigma_hat_, bool update_sigma_)
    : Concrete<AFTFamily>(y_, w_, 1), dist(dist_), event(event_),
      sigma(sigma_hat_), sigma_hat(sigma_hat_), update_sigma(update_sigma_) {}

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

  // The Weibull log-likelihood is *exactly* of the exponential form, censored
  // observations included:
  //
  //   delta (r - log sigma) - exp(r) = c - (delta / sigma) eta
  //                                      - exp(y / sigma) exp(-eta / sigma),
  //
  // with r = (y - eta) / sigma. Censoring only sets delta to zero, which drops
  // the linear term and leaves the shape intact, and the rate -1/sigma is the
  // same for every observation, which is what the form requires. So this needs
  // no data augmentation to escape the general path -- it only needed the rate
  // to be something other than +-1. The other two error distributions have no
  // such form and go through augmentation instead.
  TargetForm target_form(int h) const override {
    return dist == WEIBULL ? TARGET_EXP_DOWN : TARGET_GENERAL;
  }

  double exp_rate(int h) const override {
    return dist == WEIBULL ? -1.0 / sigma : 0.0;
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

  // The score and the information together, from one evaluation of the
  // transcendental rather than two. The base class's default forms them through
  // separate virtual calls to the two functions above, which for this family
  // means computing `r` twice and the exponential, the logistic or the inverse
  // Mills ratio twice. Every expression below is the one the function above it
  // forms, in the same order, so the two agree to the last bit -- which is what
  // the bit-identity harness in the test suite checks.
  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double r = (y(i) - eta[0]) / sigma;
    double delta = event(i);
    double s2 = sigma * sigma;

    switch (dist) {
    case WEIBULL: {
      double e = std::exp(r);
      *d1 = (e - delta) / sigma;
      *d2 = e / s2;
      break;
    }
    case LOGLOGISTIC: {
      double p = expit(r);
      *d1 = ((1.0 + delta) * p - delta) / sigma;
      *d2 = (1.0 + delta) * p * (1.0 - p) / s2;
      break;
    }
    default: {
      if (delta > 0.0) {
        *d1 = r / sigma;
        *d2 = 1.0 / s2;
      }
      else {
        double lambda = inv_mills(r);
        *d1 = lambda / sigma;
        *d2 = lambda * (lambda - r) / s2;
      }
      break;
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
// The log-normal and log-logistic accelerated failure time models, by imputing
// the censored log-failure-times. Both are quadratic once the imputation is in
// hand, and unlike the exponential form that route survives soft rules, so this
// is the version that helps at the default gate.
//
// Right-censoring is what makes the direct likelihood awkward: an observed
// failure contributes a density in (y - eta), but a censored one contributes a
// survival function, and the two have different shapes in eta. Imputing the
// failure time above its censoring time replaces the survival term with a
// density, and then every observation contributes the same shape.
//
// The scale is drawn from the *observed*-data likelihood, with the imputations
// integrated out, and only then are they redrawn -- the partially collapsed
// order of Van Dyk and Park (2008). Conditioning sigma on the current
// imputations instead would be valid but slower to mix, and the marginal draw
// costs nothing extra here because it is the same sum over observations either
// way.
// ---------------------------------------------------------------------------

struct LognormalAFTAugmentedFamily : Concrete<LognormalAFTAugmentedFamily> {
  arma::vec obs;      // log of the observed time, a failure or a censoring
  arma::vec event;    // 1 for an observed failure, 0 for right-censored
  arma::vec latent;   // the imputed log failure time; equal to obs when observed
  double sigma;
  double sigma_hat;
  bool update_sigma;

  LognormalAFTAugmentedFamily(const arma::vec& y_, const arma::vec& w_,
                              const arma::vec& event_, double sigma_hat_,
                              bool update_sigma_)
    : Concrete<LognormalAFTAugmentedFamily>(y_, w_, 1), obs(y_), event(event_),
      sigma(sigma_hat_), sigma_hat(sigma_hat_), update_sigma(update_sigma_) {
    // A deterministic start one scale unit past each censoring time, replaced by
    // a proper draw at the end of the first sweep.
    latent = obs;
    for (int i = 0; i < N; i++) {
      if (event(i) <= 0.0) {
        latent(i) += sigma;
      }
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

  TargetForm target_form(int h) const override { return TARGET_QUADRATIC; }

  double logdens_unit(int i, const double* eta) const override {
    double r = (latent(i) - eta[0]) / sigma;
    return -0.5 * r * r;
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double prec = 1.0 / (sigma * sigma);
    *d1 = (latent(i) - eta[0]) * prec;
    *d2 = prec;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return (latent(i) - eta[0]) / (sigma * sigma);
  }

  double info_unit(int i, const double* eta, int h) const override {
    return 1.0 / (sigma * sigma);
  }

  void update_aux(const arma::mat& eta) override {
    const arma::rowvec& e = eta.row(0);

    if (update_sigma) {
      auto logf = [this, &e](double log_sigma) {
        double s = std::exp(log_sigma);
        if (!(s > 0.0) || !std::isfinite(s)) {
          return R_NegInf;
        }
        double out = Rf_dcauchy(s, 0.0, sigma_hat, 1) + log_sigma;
        for (int i = 0; i < N; i++) {
          out += w(i) * AFTFamily::loglik_one(AFTFamily::LOGNORMAL, obs(i),
                                              e(i), event(i), s);
        }
        return out;
      };
      sigma = std::exp(slice_sampler(std::log(sigma), logf, 0.5, -20.0, 20.0));
    }

    for (int i = 0; i < N; i++) {
      if (event(i) > 0.0) {
        latent(i) = obs(i);
        continue;
      }
      // The failure happened after the censoring time, so the error is a
      // standard normal truncated below at the residual it would have had.
      double lo = (obs(i) - e(i)) / sigma;
      latent(i) = e(i) + sigma * truncated_normal_between(lo, R_PosInf);
    }

    y = latent;
  }

  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      out += w(i) * AFTFamily::loglik_one(AFTFamily::LOGNORMAL, obs(i), e(i),
                                          event(i), sigma);
    }

    return out;
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

// The log-logistic model needs a second layer: even with the failure time in
// hand the complete-data density is logistic, not normal. Writing it as a scale
// mixture of normals -- Polson, Scott and Windle (2013), Theorem 1 at a = 1 and
// b = 2, where the tilting term kappa = a - b / 2 vanishes -- gives each
// observation its own precision and makes the target quadratic. The same device
// the ordinal logit uses.
struct LoglogisticAFTAugmentedFamily
  : Concrete<LoglogisticAFTAugmentedFamily> {
  arma::vec obs;
  arma::vec event;
  arma::vec latent;
  arma::vec omega;    // the Polya-Gamma precision, one per observation
  double sigma;
  double sigma_hat;
  bool update_sigma;

  static constexpr double OMEGA_MIN = 1e-10;

  LoglogisticAFTAugmentedFamily(const arma::vec& y_, const arma::vec& w_,
                                const arma::vec& event_, double sigma_hat_,
                                bool update_sigma_)
    : Concrete<LoglogisticAFTAugmentedFamily>(y_, w_, 1), obs(y_),
      event(event_), sigma(sigma_hat_), sigma_hat(sigma_hat_),
      update_sigma(update_sigma_) {
    latent = obs;
    omega.set_size(N);
    for (int i = 0; i < N; i++) {
      if (event(i) <= 0.0) {
        latent(i) += sigma;
      }
      // The mean of PG(2, 0), replaced by a proper draw at the end of the first
      // sweep.
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

  TargetForm target_form(int h) const override { return TARGET_QUADRATIC; }

  double logdens_unit(int i, const double* eta) const override {
    double r = (latent(i) - eta[0]) / sigma;
    return -0.5 * omega(i) * r * r;
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double prec = omega(i) / (sigma * sigma);
    *d1 = (latent(i) - eta[0]) * prec;
    *d2 = prec;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return omega(i) * (latent(i) - eta[0]) / (sigma * sigma);
  }

  double info_unit(int i, const double* eta, int h) const override {
    return omega(i) / (sigma * sigma);
  }

  void update_aux(const arma::mat& eta) override {
    const arma::rowvec& e = eta.row(0);

    if (update_sigma) {
      auto logf = [this, &e](double log_sigma) {
        double s = std::exp(log_sigma);
        if (!(s > 0.0) || !std::isfinite(s)) {
          return R_NegInf;
        }
        double out = Rf_dcauchy(s, 0.0, sigma_hat, 1) + log_sigma;
        for (int i = 0; i < N; i++) {
          out += w(i) * AFTFamily::loglik_one(AFTFamily::LOGLOGISTIC, obs(i),
                                              e(i), event(i), s);
        }
        return out;
      };
      sigma = std::exp(slice_sampler(std::log(sigma), logf, 0.5, -20.0, 20.0));
    }

    for (int i = 0; i < N; i++) {
      if (event(i) > 0.0) {
        latent(i) = obs(i);
      }
      else {
        // Drawn from the logistic truncated below, by inverting its cumulative
        // distribution rather than conditioning on omega: the two are both
        // valid Gibbs steps, and this one leaves the imputation independent of
        // the precision it is about to be paired with. Carrying the *upper*
        // tail probability keeps the inversion accurate when the censoring time
        // is far above the predictor and the lower probability rounds to one.
        double q = expit(-(obs(i) - e(i)) / sigma) * (1.0 - unif_rand());
        if (!(q > 0.0)) {
          q = std::numeric_limits<double>::min();
        }
        latent(i) = e(i) + sigma * (std::log1p(-q) - std::log(q));
      }

      // Then the precision, given the residual now in place.
      double r = (latent(i) - e(i)) / sigma;
      // Compared rather than passed to std::max(), which takes its arguments
      // by const reference: binding OMEGA_MIN to one is an ODR-use, and a
      // `static constexpr` member is only implicitly inline from C++17 on. Under
      // an older standard the symbol has no definition, which the optimizer
      // hides by folding the constant and a -O0 build does not.
      double drawn = rpg(2.0, std::fabs(r));
      omega(i) = drawn < OMEGA_MIN ? OMEGA_MIN : drawn;
    }

    y = latent;
  }

  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      out += w(i) * AFTFamily::loglik_one(AFTFamily::LOGLOGISTIC, obs(i), e(i),
                                          event(i), sigma);
    }

    return out;
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
// Proportional hazards with a piecewise-constant baseline: the semiparametric
// survival model, with the log hazard ratio given a forest.
//
//   lambda(t | x) = lambda_0(t) exp(r(x)),   lambda_0(t) = lambda_b for t in bin b.
//
// Cox's *partial* likelihood couples observations through risk sets and so does
// not decompose into a sum over the observations reaching a leaf, which is what
// this sampler needs. The full likelihood of the piecewise-exponential model
// does decompose, and it approaches the partial likelihood as the bins shrink
// (Sinha, Ibrahim and Chen 2003), so proportional hazards is reachable after all
// -- and cheaply. Writing Lambda_0 for the cumulative baseline, the contribution
// of observation i is
//
//   delta_i (log lambda_{b_i} + eta) - Lambda_0(y_i) exp(eta),
//
// which is `a eta + b exp(eta)` with the rate exactly +1: the same exponential
// form the Poisson has, so a leaf update is one pass over the node rather than
// one per trial value. The baseline is a nuisance vector drawn from its exact
// gamma conditional, which is what makes it as cheap as it is flexible.
//
// The level of the predictor and the level of the baseline are identified only
// jointly, so the convention here is that the baseline carries it: the predictor
// starts at zero and is reported as a log hazard ratio against the fitted
// baseline. Basak, Linero, Maringe and Rubio (2024) is the reference for this
// construction, in the relative-survival setting.
// ---------------------------------------------------------------------------

struct PHFamily : Concrete<PHFamily> {
  arma::vec event;
  int num_bins;
  arma::vec edges;        // lower edge of each bin; edges(0) is 0
  arma::vec width;        // bin widths, the last one unused and set to zero
  arma::uvec bin_of;      // which bin each observation's time falls in
  arma::vec lambda;
  arma::vec cum_base;     // Lambda_0(y_i), rebuilt whenever lambda moves
  arma::vec events_in;    // A_b, free of everything drawn, so computed once
  double shape;           // a_lambda
  double rate;            // b_lambda, itself drawn
  bool update_lambda;

  // Scratch for the O(N + B) baseline update, allocated once.
  mutable arma::vec bin_total;
  mutable arma::vec bin_partial;

  // R draws gamma variates by shape and *scale*; every conditional here is
  // written by rate, which is the convention the model is stated in.
  static double gamma_by_rate(double shape_in, double rate_in) {
    if (!(rate_in > 0.0) || !std::isfinite(rate_in)) {
      return 0.0;
    }
    return Rf_rgamma(shape_in, 1.0 / rate_in);
  }

  PHFamily(const arma::vec& y_, const arma::vec& w_, const arma::vec& event_,
           const arma::vec& edges_, double shape_, double rate_,
           bool update_lambda_)
    : Concrete<PHFamily>(y_, w_, 1), event(event_), edges(edges_),
      shape(shape_), rate(rate_), update_lambda(update_lambda_) {
    num_bins = static_cast<int>(edges.n_elem);

    width.set_size(num_bins);
    for (int b = 0; b < num_bins - 1; b++) {
      width(b) = edges(b + 1) - edges(b);
    }
    // The last bin runs to infinity, but no observation lies *above* it, so its
    // width never multiplies anything.
    width(num_bins - 1) = 0.0;

    // Binary search rather than a scan: with one bin per event time a scan is
    // quadratic in the sample size.
    bin_of.set_size(N);
    for (int i = 0; i < N; i++) {
      const double* first = edges.memptr();
      const double* found = std::upper_bound(first, first + num_bins, y(i));
      bin_of(i) = static_cast<arma::uword>(
        std::max<std::ptrdiff_t>(found - first - 1, 0));
    }

    events_in.zeros(num_bins);
    for (int i = 0; i < N; i++) {
      events_in(bin_of(i)) += w(i) * event(i);
    }

    lambda.set_size(num_bins);
    lambda.fill(shape / std::max(rate, 1e-8));
    bin_total.set_size(num_bins);
    bin_partial.set_size(num_bins);
    cum_base.set_size(N);
    rebuild_baseline();
  }

  // Lambda_0(y) = sum of the whole bins below y, plus the part of y's own bin it
  // reaches. One O(B) scan and one O(N) pass.
  void rebuild_baseline() {
    arma::vec before(num_bins, arma::fill::zeros);
    for (int b = 1; b < num_bins; b++) {
      before(b) = before(b - 1) + lambda(b - 1) * width(b - 1);
    }
    for (int i = 0; i < N; i++) {
      arma::uword b = bin_of(i);
      cum_base(i) = before(b) + lambda(b) * (y(i) - edges(b));
    }
  }

  TargetForm target_form(int h) const override { return TARGET_EXP_UP; }

  // The `delta * log lambda` term moves with the baseline, not the predictor, so
  // it belongs to the part the sampler caches per sweep.
  arma::vec compute_eta_free() const override {
    arma::vec out(N);
    for (int i = 0; i < N; i++) {
      out(i) = event(i) > 0.0 ? std::log(lambda(bin_of(i))) : 0.0;
    }
    return out;
  }

  double logdens_unit(int i, const double* eta) const override {
    return event(i) * eta[0] - cum_base(i) * std::exp(eta[0]);
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double mu = cum_base(i) * std::exp(eta[0]);
    *d1 = event(i) - mu;
    *d2 = mu;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return event(i) - cum_base(i) * std::exp(eta[0]);
  }

  double info_unit(int i, const double* eta, int h) const override {
    return cum_base(i) * std::exp(eta[0]);
  }

  void update_aux(const arma::mat& eta) override {
    if (!update_lambda) {
      return;
    }

    const arma::rowvec& e = eta.row(0);

    // B_b = width_b * (total exposure of everyone who outlives the bin)
    //       + (the part-bins of everyone who leaves inside it).
    bin_total.zeros();
    bin_partial.zeros();
    for (int i = 0; i < N; i++) {
      double contribution = w(i) * std::exp(e(i));
      arma::uword b = bin_of(i);
      bin_total(b) += contribution;
      bin_partial(b) += contribution * (y(i) - edges(b));
    }

    double above = 0.0;
    for (int b = num_bins - 1; b >= 0; b--) {
      // `above` is the exposure of the observations whose bin is strictly higher.
      double exposure = width(b) * above + bin_partial(b);
      lambda(b) = gamma_by_rate(shape + events_in(b), rate + exposure);
      above += bin_total(b);
    }

    // A flat prior on the baseline's own rate leaves a gamma conditional, which
    // is what lets the bins borrow a level from each other rather than each
    // resting on its own handful of events.
    rate = gamma_by_rate(num_bins + 1.0, arma::accu(lambda));

    rebuild_baseline();
    refresh_eta_free();
  }

  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;
    for (int i = 0; i < N; i++) {
      out += w(i) * (event(i) * (std::log(lambda(bin_of(i))) + e(i)) -
                     cum_base(i) * std::exp(e(i)));
    }
    return out;
  }

  std::vector<std::string> aux_names() const override {
    std::vector<std::string> out;
    for (int b = 0; b < num_bins; b++) {
      out.push_back("lambda" + std::to_string(b + 1));
    }
    out.push_back("lambda_rate");
    return out;
  }

  arma::vec aux_values() const override {
    arma::vec out(num_bins + 1);
    out.head(num_bins) = lambda;
    out(num_bins) = rate;
    return out;
  }

  void set_aux(const arma::vec& values) override {
    if (static_cast<int>(values.n_elem) >= num_bins) {
      lambda = values.head(num_bins);
      if (static_cast<int>(values.n_elem) > num_bins) {
        rate = values(num_bins);
      }
      rebuild_baseline();
      refresh_eta_free();
    }
  }
};

// ---------------------------------------------------------------------------
// Gaussian location-scale regression: one forest for the mean and a second for
// the log standard deviation, so the variance is an unrestricted function of
// the predictors.
// ---------------------------------------------------------------------------

struct LocationScaleFamily : Concrete<LocationScaleFamily> {

  // Quadratic in the mean, and in the log standard deviation the *exponential*
  // form at rate -2: the log density is
  //
  //   const - eta1 - (y - eta0)^2 exp(-2 eta1) / 2,
  //
  // which is c + a eta1 + b exp(-2 eta1) with b free of eta1. So neither
  // predictor needs the general path, which is why this is asked per predictor
  // rather than per family. The rate is what made this reachable: the form used
  // to be hardcoded to exp(+-eta).
  TargetForm target_form(int h) const override {
    return h == 0 ? TARGET_QUADRATIC : TARGET_EXP_DOWN;
  }

  double exp_rate(int h) const override {
    return h == 0 ? 0.0 : -2.0;
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
    // The *observed* information for the log scale. The expected version is
    // constant at 2 and was used because it is quieter, but the exponential form
    // reads its coefficients off the curvature and so needs the curvature the
    // target actually has. Nothing is approximated either way: a Laplace fit is
    // a proposal, and here it stops being an approximation at all.
    double r = (y(i) - eta[0]) * std::exp(-eta[1]);
    return 2.0 * r * r;
  }
};

// ---------------------------------------------------------------------------
// Zero-inflated counts. Two additive predictors: the first is the log mean of
// the count component, the second the logit of the probability that an
// observation is a structural zero. Both are modeled nonparametrically, so the
// excess-zero mechanism is free to depend on the predictors.
// ---------------------------------------------------------------------------

// The zero-inflated log likelihood, as free functions so that the direct family
// and the augmented one below share exactly one definition of the model.
namespace {

// log P(the count component yields a zero).
double zi_log_p0(double eta_count, double theta, bool negbin) {
  double mu = std::exp(eta_count);

  if (!negbin) {
    return -mu;
  }

  return theta * (std::log(theta) - std::log(theta + mu));
}

double zi_count_logpmf(double y_i, double eta_count, double theta,
                       bool negbin) {
  double mu = std::exp(eta_count);

  if (!negbin) {
    return y_i * eta_count - mu - R::lgammafn(y_i + 1.0);
  }

  double log_theta_mu = std::log(theta + mu);

  return R::lgammafn(y_i + theta) - R::lgammafn(theta) -
    R::lgammafn(y_i + 1.0) + theta * (std::log(theta) - log_theta_mu) +
    y_i * (eta_count - log_theta_mu);
}

double zi_loglik_one(double y_i, const double* eta, double theta,
                     bool negbin) {
  if (y_i > 0.0) {
    return log1m_expit(eta[1]) + zi_count_logpmf(y_i, eta[0], theta, negbin);
  }

  // A zero arises either structurally or from the count component.
  return log_sum_exp(log_expit(eta[1]),
                     log1m_expit(eta[1]) + zi_log_p0(eta[0], theta, negbin));
}

} // namespace

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

  double log_p0(double eta_count, double th) const {
    return zi_log_p0(eta_count, th, negbin);
  }

  double loglik_one(double y_i, const double* eta, double th) const {
    return zi_loglik_one(y_i, eta, th, negbin);
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
// Beta regression for a response strictly inside the unit interval. The forest
// is on the logit of the mean and a precision is drawn, so the two shapes are
// mu * phi and (1 - mu) * phi and their sum is free of the predictor -- which is
// what makes the normalizing constant partly cacheable.
//
// This is the interior of `ordbeta()` with the endpoint machinery removed. Kept
// as its own family rather than as a special case of that one because the two
// answer different questions: this one says the response cannot reach 0 or 1,
// and ordered beta says it can and estimates how often it does. Fitting ordered
// beta to data with no boundary observations leaves its two cutpoints with
// nothing to identify them.
// ---------------------------------------------------------------------------

struct BetaFamily : Family {
  double phi;
  double prior_shape;
  double prior_rate;
  bool update_phi;

  // log(y) and log(1 - y) are multiplied by shapes that move with eta, so they
  // cannot go into the eta-free part, but they are fixed per observation.
  arma::vec log_y;
  arma::vec log1m_y;
  arma::vec logit_y;

  // The derivatives need two digamma and two trigamma evaluations per
  // observation per pass, which made this the slowest family in the package. But
  // both of the combinations they appear in,
  //
  //   psi(mu phi) - psi((1 - mu) phi)   and   psi'(mu phi) + psi'((1 - mu) phi),
  //
  // are functions of the *single scalar* mu, because the two shapes always sum to
  // phi and phi is fixed for the whole of a sweep -- it moves only in
  // `update_aux`. So they are tabulated once per sweep on a grid in the additive
  // predictor and interpolated, which turns four special-function calls into two
  // loads and a multiply.
  //
  // This is a proposal, not the target: `logdens_unit` below stays exact, and the
  // Metropolis step corrects. What the sampler requires of a Laplace fit is that
  // it be a *deterministic* function of the current state, so that the forward
  // and reverse moves rebuild the same proposal, and a fixed table is exactly
  // that. Interpolation error costs acceptance rate, not correctness.
  static const int TAB_N = 2049;
  static constexpr double TAB_L = 8.0;
  std::vector<double> tab_dpsi;
  std::vector<double> tab_tpsi;
  double tab_step;
  double tab_phi;

  BetaFamily(const arma::vec& y_, const arma::vec& w_, double phi_,
             double prior_shape_, double prior_rate_, bool update_phi_)
    : Family(y_, w_, 1), phi(phi_), prior_shape(prior_shape_),
      prior_rate(prior_rate_), update_phi(update_phi_) {
    log_y.set_size(N);
    log1m_y.set_size(N);
    logit_y.set_size(N);
    for (int i = 0; i < N; i++) {
      log_y(i) = std::log(y_(i));
      log1m_y(i) = std::log1p(-y_(i));
      logit_y(i) = log_y(i) - log1m_y(i);
    }
    tab_dpsi.resize(TAB_N);
    tab_tpsi.resize(TAB_N);
    tab_step = 2.0 * TAB_L / (TAB_N - 1);
    build_tables();
  }

  // Exact, for the grid itself and for predictors beyond it.
  void psi_exact(double e, double* dpsi, double* tpsi) const {
    double mu = expit(e);
    double a = mu * phi;
    double b = phi - a;
    *dpsi = R::digamma(a) - R::digamma(b);
    *tpsi = R::trigamma(a) + R::trigamma(b);
  }

  void build_tables() {
    for (int g = 0; g < TAB_N; g++) {
      psi_exact(-TAB_L + g * tab_step, &tab_dpsi[g], &tab_tpsi[g]);
    }
    tab_phi = phi;
  }

  void psi_lookup(double e, double* dpsi, double* tpsi) const {
    if (!(e > -TAB_L) || !(e < TAB_L)) {
      // The tails, where the grid would need to resolve a function growing like
      // exp(|e|) / phi. Rare enough to pay for exactly.
      psi_exact(e, dpsi, tpsi);
      return;
    }
    double t = (e + TAB_L) / tab_step;
    int g = static_cast<int>(t);
    // The guard above makes g <= TAB_N - 2 mathematically, but rounding in the
    // division could land it one past that, and the read below is of g + 1.
    if (g > TAB_N - 2) {
      g = TAB_N - 2;
    }
    double frac = t - g;
    *dpsi = tab_dpsi[g] + frac * (tab_dpsi[g + 1] - tab_dpsi[g]);
    *tpsi = tab_tpsi[g] + frac * (tab_tpsi[g + 1] - tab_tpsi[g]);
  }

  // lbeta(a, b) = lgamma(a) + lgamma(b) - lgamma(phi), and the last term is free
  // of eta because the two shapes always sum to phi. The unit constants of the
  // density go the same way.
  arma::vec compute_eta_free() const override {
    arma::vec out(N);
    double lg_phi = R::lgammafn(phi);
    for (int i = 0; i < N; i++) {
      out(i) = lg_phi - log_y(i) - log1m_y(i);
    }
    return out;
  }

  double logdens_unit(int i, const double* eta) const override {
    double mu = expit(eta[0]);
    double a = mu * phi;
    double b = phi - a;

    if (!(a > 0.0) || !(b > 0.0)) {
      return R_NegInf;
    }

    return a * log_y(i) + b * log1m_y(i) - R::lgammafn(a) - R::lgammafn(b);
  }

  // Complete log density, for the precision update, which varies phi and so
  // cannot use the cached eta-free part.
  static double loglik_one(double y_i, double eta, double phi_) {
    double mu = expit(eta);
    double a = mu * phi_;
    double b = phi_ - a;

    if (!(a > 0.0) || !(b > 0.0)) {
      return R_NegInf;
    }

    return (a - 1.0) * std::log(y_i) + (b - 1.0) * std::log1p(-y_i) -
      Rf_lbeta(a, b);
  }

  // The two digamma contributions from log Beta(a, b) keep only their
  // difference, because the shapes move in opposite directions: da/deta is
  // -db/deta. The information is the expected beta information
  // s^2 (psi'(a) + psi'(b)) rather than the observed one, whose extra term has
  // mean zero but can turn the total negative.
  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double e = eta[0];
    double mu = expit(e);
    double slope = phi * mu * (1.0 - mu);
    double dpsi;
    double tpsi;
    psi_lookup(e, &dpsi, &tpsi);

    *d1 = slope * (logit_y(i) - dpsi);
    *d2 = slope * slope * tpsi;
  }

  
  // Both of these default to a central difference in the base class. Delegating
  // to the analytic pair keeps the diagnostic entry point on the same derivatives
  // the sampler actually uses, which is what lets the derivative tests check
  // them against a difference rather than against another difference.
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

  void update_aux(const arma::mat& eta) override {
    if (!update_phi) {
      return;
    }

    const arma::rowvec& e = eta.row(0);

    phi = std::exp(slice_sampler(std::log(phi), [this, &e](double log_phi) {
      double p = std::exp(log_phi);
      if (!(p > 0.0) || !std::isfinite(p)) {
        return R_NegInf;
      }
      // Gamma prior on the precision, with the Jacobian for the log scale
      // folded into the shape term, as `ordbeta()` does.
      // lgamma(a + b) = lgamma(phi) is free of the observation, so it comes
      // out of the loop; that is one of the three log-gammas `lbeta` does.
      double lg_p = R::lgammafn(p);
      double out = prior_shape * log_phi - prior_rate * p;
      for (int i = 0; i < N; i++) {
        double mu = expit(e(i));
        double a = mu * p;
        double b = p - a;
        if (!(a > 0.0) || !(b > 0.0)) {
          return R_NegInf;
        }
        out += w(i) * ((a - 1.0) * log_y(i) + (b - 1.0) * log1m_y(i) + lg_p -
                       R::lgammafn(a) - R::lgammafn(b));
      }
      return out;
    }, 0.5, -10.0, 15.0));

    build_tables();
    refresh_eta_free();
  }

  std::vector<std::string> aux_names() const override {
    return std::vector<std::string>{"phi"};
  }

  arma::vec aux_values() const override { return arma::vec{phi}; }

  void set_aux(const arma::vec& values) override {
    if (values.n_elem > 0) {
      phi = values(0);
      if (phi != tab_phi) {
        build_tables();
      }
    }
  }
};

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
  arma::vec logit_y;

  // The same tabulation `BetaFamily` uses, for the same reason: the interior of
  // this likelihood *is* a beta density, so its derivatives carry the same two
  // digamma and two trigamma evaluations per observation per pass, and the same
  // two combinations are functions of mu alone once phi is fixed.
  static const int TAB_N = 2049;
  static constexpr double TAB_L = 8.0;
  std::vector<double> tab_dpsi;
  std::vector<double> tab_tpsi;
  double tab_step;
  double tab_phi;

  OrdBetaFamily(const arma::vec& y_, const arma::vec& w_, double cut1_,
                double cut2_, double phi_, double prior_shape_,
                double prior_rate_, bool update_phi_)
    : Family(y_, w_, 1), cut1(cut1_), cut2(cut2_), phi(phi_),
      prior_shape(prior_shape_), prior_rate(prior_rate_),
      update_phi(update_phi_) {
    log_y.set_size(N);
    log1m_y.set_size(N);
    logit_y.set_size(N);
    for (int i = 0; i < N; i++) {
      bool interior = y_(i) > 0.0 && y_(i) < 1.0;
      log_y(i) = interior ? std::log(y_(i)) : 0.0;
      log1m_y(i) = interior ? std::log1p(-y_(i)) : 0.0;
      logit_y(i) = log_y(i) - log1m_y(i);
    }
    tab_dpsi.resize(TAB_N);
    tab_tpsi.resize(TAB_N);
    tab_step = 2.0 * TAB_L / (TAB_N - 1);
    build_tables();
  }

  void psi_exact(double e, double* dpsi, double* tpsi) const {
    double mu = expit(e);
    double a = mu * phi;
    double b = phi - a;
    *dpsi = R::digamma(a) - R::digamma(b);
    *tpsi = R::trigamma(a) + R::trigamma(b);
  }

  void build_tables() {
    for (int g = 0; g < TAB_N; g++) {
      psi_exact(-TAB_L + g * tab_step, &tab_dpsi[g], &tab_tpsi[g]);
    }
    tab_phi = phi;
  }

  void psi_lookup(double e, double* dpsi, double* tpsi) const {
    if (!(e > -TAB_L) || !(e < TAB_L)) {
      psi_exact(e, dpsi, tpsi);
      return;
    }
    double t = (e + TAB_L) / tab_step;
    int g = static_cast<int>(t);
    // The guard above makes g <= TAB_N - 2 mathematically, but rounding in the
    // division could land it one past that, and the read below is of g + 1.
    if (g > TAB_N - 2) {
      g = TAB_N - 2;
    }
    double frac = t - g;
    *dpsi = tab_dpsi[g] + frac * (tab_dpsi[g + 1] - tab_dpsi[g]);
    *tpsi = tab_tpsi[g] + frac * (tab_tpsi[g + 1] - tab_tpsi[g]);
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
    double slope = phi * mu * (1.0 - mu);
    double dpsi;
    double tpsi;
    psi_lookup(e, &dpsi, &tpsi);

    *d1 = d1_span + slope * (logit_y(i) - dpsi);
    *d2 = d2_span + slope * slope * tpsi;
  }

  
  // Both of these default to a central difference in the base class. Delegating
  // to the analytic pair keeps the diagnostic entry point on the same derivatives
  // the sampler actually uses, which is what lets the derivative tests check
  // them against a difference rather than against another difference.
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

  void update_aux(const arma::mat& eta) override {
    const arma::rowvec e = eta.row(0);

    // Slice sampling only needs the log density up to an additive constant, and
    // this likelihood splits: the endpoint-and-middle part depends on the
    // cutpoints but not on phi, and the beta density on the interior depends on
    // phi but not on the cutpoints. So each update evaluates its own half and
    // drops the other as a constant, instead of rebuilding the whole likelihood
    // -- three log-gamma calls per observation -- on every slice evaluation.
    auto cut_total = [this, &e](double c1, double c2) {
      double out = 0.0;
      for (int i = 0; i < N; i++) {
        double y_i = y(i);
        if (y_i <= 0.0) {
          out += w(i) * log1m_expit(e(i) - c1);
        }
        else if (y_i >= 1.0) {
          out += w(i) * log_expit(e(i) - c2);
        }
        else {
          out += w(i) * log_diff_logistic(e(i) - c2, e(i) - c1);
        }
      }
      return out;
    };

    auto phi_total = [this, &e](double p) {
      // lgamma(a + b) = lgamma(phi) is free of the observation, so it comes out
      // of the loop; that is one of the three log-gammas `lbeta` would do.
      double lg_p = R::lgammafn(p);
      double out = 0.0;
      for (int i = 0; i < N; i++) {
        double y_i = y(i);
        if (!(y_i > 0.0 && y_i < 1.0)) {
          continue;
        }
        double mu = expit(e(i));
        double a = mu * p;
        double b = p - a;
        if (!(a > 0.0) || !(b > 0.0)) {
          return R_NegInf;
        }
        out += w(i) * ((a - 1.0) * log_y(i) + (b - 1.0) * log1m_y(i) + lg_p -
                       R::lgammafn(a) - R::lgammafn(b));
      }
      return out;
    };

    double c2_now = cut2;
    cut1 = slice_sampler(cut1, [&](double v) {
      if (!(v < c2_now)) {
        return R_NegInf;
      }
      return cut_total(v, c2_now);
    }, 0.5, -30.0, c2_now);

    double c1_now = cut1;
    cut2 = slice_sampler(cut2, [&](double v) {
      if (!(v > c1_now)) {
        return R_NegInf;
      }
      return cut_total(c1_now, v);
    }, 0.5, c1_now, 30.0);

    if (update_phi) {
      phi = std::exp(slice_sampler(std::log(phi), [&](double log_phi) {
        double p = std::exp(log_phi);
        if (!(p > 0.0) || !std::isfinite(p)) {
          return R_NegInf;
        }
        return prior_shape * log_phi - prior_rate * p + phi_total(p);
      }, 0.5, -10.0, 15.0));
    }

    build_tables();
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
      if (phi != tab_phi) {
        build_tables();
      }
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

// The zero-inflated families, rewritten with two latent variables so that both
// of their forests get an exploitable target.
//
// What blocks the direct family is the zero: `log[pi + (1 - pi) P_0(mu)]` is a
// log-sum-exp of the two components, so neither predictor has a shape and every
// trial value of a leaf parameter costs its own pass. Introducing the indicator
// the mixture is a mixture *over* -- z_i = 1 when observation i is a structural
// zero -- separates them. Conditional on z the two forests see:
//
//   count:     prod over {z = 0} of the count likelihood, a plain Poisson or
//              negative binomial, so the exponential form applies;
//   inflation: a Bernoulli logistic likelihood in z, so Polya-Gamma applies and
//              the target is quadratic.
//
// z is drawn from its exact conditional: it is zero whenever y > 0, and for
// y = 0 it is one with probability pi / (pi + (1 - pi) P_0), where P_0 is the
// true count probability of a zero.
//
// For the negative binomial a second augmentation goes on top of the first: the
// non-structural observations are given the Poisson-gamma rate of
// NegBinAugmentedFamily, which turns their target from a general one into the
// exponential form. z is drawn with that rate integrated out and the rate is
// redrawn immediately afterwards, which is a valid partially collapsed Gibbs
// step in that order (Van Dyk and Park 2008) and mixes better than conditioning
// z on a stale rate. theta is drawn the same way the direct family draws it,
// from the true zero-inflated likelihood.
//
// The rate of a structural zero is never used -- its observation's contribution
// to the count target is multiplied by (1 - z) -- so it is not drawn.
struct ZeroInflatedAugmentedFamily : Concrete<ZeroInflatedAugmentedFamily> {
  arma::vec count;
  arma::vec structural;   // z, one when the observation is a structural zero
  arma::vec rate;         // lambda, the Poisson rate, negative binomial only
  arma::vec kappa;        // Polya-Gamma constants for the inflation forest
  arma::vec omega;
  bool negbin;
  double theta;
  double log_theta;
  double prior_shape;
  double prior_rate;
  bool update_theta;

  ZeroInflatedAugmentedFamily(const arma::vec& y_, const arma::vec& w_,
                              bool negbin_, double theta_, double prior_shape_,
                              double prior_rate_, bool update_theta_)
    : Concrete<ZeroInflatedAugmentedFamily>(y_, unit_weights(y_.n_elem), 2),
      count(y_), negbin(negbin_), theta(theta_), prior_shape(prior_shape_),
      prior_rate(prior_rate_), update_theta(update_theta_) {
    log_theta = std::log(theta);
    structural = arma::vec(N, arma::fill::zeros);
    kappa = arma::vec(N, arma::fill::zeros);
    omega = arma::vec(N, arma::fill::ones);
    // Replaced by a proper draw before the first forest moves. The counts are
    // the natural guess at their own rates.
    rate = arma::clamp(y_, 0.1, arma::datum::inf);
  }

  static bool applies(const arma::vec& w) {
    return arma::all(w == 1.0);
  }

  // The count forest inherits the shape of whatever the count likelihood is
  // once the structural zeros are out of it: a Poisson rises with the
  // predictor, and the Poisson-gamma rewriting of a negative binomial falls
  // with it, exactly as for the unmixed families.
  TargetForm target_form(int h) const override {
    if (h == 1) {
      return TARGET_QUADRATIC;
    }
    return negbin ? TARGET_EXP_DOWN : TARGET_EXP_UP;
  }

  double count_kernel(int i, double eta_count) const {
    if (negbin) {
      return -theta * (eta_count + rate(i) * std::exp(-eta_count));
    }
    return count(i) * eta_count - std::exp(eta_count);
  }

  double logdens_unit(int i, const double* eta) const override {
    double live = 1.0 - structural(i);
    return live * count_kernel(i, eta[0]) +
      kappa(i) * eta[1] - 0.5 * omega(i) * eta[1] * eta[1];
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    if (h == 1) {
      *d1 = kappa(i) - omega(i) * eta[1];
      *d2 = omega(i);
      return;
    }

    double live = 1.0 - structural(i);

    if (negbin) {
      double scaled = theta * rate(i) * std::exp(-eta[0]);
      *d1 = live * (scaled - theta);
      *d2 = live * scaled;
      return;
    }

    double mu = std::exp(eta[0]);
    *d1 = live * (count(i) - mu);
    *d2 = live * mu;
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

  // z | y, eta, theta, with the Poisson-gamma rate integrated out. A positive
  // count cannot be structural; a zero is structural with the odds the two
  // components give it.
  void draw_structural(const arma::mat& eta) {
    for (int i = 0; i < N; i++) {
      if (count(i) > 0.0) {
        structural(i) = 0.0;
        continue;
      }

      double log_zero = negbin
        ? theta * (log_theta - std::log(theta + std::exp(eta(0, i))))
        : -std::exp(eta(0, i));
      double log_odds = eta(1, i) - log_zero;

      structural(i) = unif_rand() < expit(log_odds) ? 1.0 : 0.0;
    }
  }

  // lambda | y, z, eta, theta is Gamma(y + theta, 1 + theta / mu) for a
  // non-structural observation, and unused for a structural one.
  void draw_rate(const arma::mat& eta) {
    for (int i = 0; i < N; i++) {
      if (structural(i) > 0.0) {
        continue;
      }

      double scale = 1.0 / (1.0 + theta * std::exp(-eta(0, i)));
      rate(i) = Rf_rgamma(count(i) + theta, scale);
    }
  }

  void before_forest(int h, const arma::mat& eta) override {
    // Redrawn before each forest rather than once a sweep, so that the
    // inflation forest sees the indicator the count forest has just moved under
    // and the other way round.
    draw_structural(eta);

    if (h == 0) {
      if (negbin) {
        draw_rate(eta);
      }
      return;
    }

    for (int i = 0; i < N; i++) {
      kappa(i) = structural(i) - 0.5;
      omega(i) = rpg(1.0, eta(1, i));
    }
  }

  void update_aux(const arma::mat& eta) override {
    if (!negbin || !update_theta) {
      return;
    }

    // The collapsed conditional: the true zero-inflated negative binomial
    // likelihood, with both latent variables integrated out, which is the same
    // target the direct family slices on.
    auto logf = [this, &eta](double candidate) {
      double th = std::exp(candidate);

      if (!(th > 0.0) || !std::isfinite(th)) {
        return R_NegInf;
      }

      double out = prior_shape * candidate - prior_rate * th;

      for (int i = 0; i < N; i++) {
        out += zi_loglik_one(count(i), eta.colptr(i), th, true);
      }

      return out;
    };

    theta = std::exp(slice_sampler(log_theta, logf, 1.0, -20.0, 20.0));
    log_theta = std::log(theta);
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
      log_theta = std::log(theta);
    }
  }

  double reported_loglik(const arma::mat& eta) const override {
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      out += zi_loglik_one(count(i), eta.colptr(i), theta, negbin);
    }

    return out;
  }
};

// A draw from the inverse Wishart distribution with `df` degrees of freedom and
// scale matrix `scatter`, by Bartlett's decomposition of the Wishart it inverts.
// The dimension here is the number of response categories less one, so these are
// small matrices and the cubic work in them is not worth avoiding.
namespace {

arma::mat inverse_wishart(double df, const arma::mat& scatter) {
  int p = static_cast<int>(scatter.n_rows);
  arma::mat inverse_scale;

  if (!arma::inv_sympd(inverse_scale, scatter)) {
    return scatter / df;
  }

  arma::mat factor;

  if (!arma::chol(factor, inverse_scale, "lower")) {
    return scatter / df;
  }

  arma::mat bartlett(p, p, arma::fill::zeros);

  for (int j = 0; j < p; j++) {
    bartlett(j, j) = std::sqrt(Rf_rchisq(df - static_cast<double>(j)));

    for (int k = 0; k < j; k++) {
      bartlett(j, k) = norm_rand();
    }
  }

  arma::mat root = factor * bartlett;
  arma::mat wishart = root * root.t();
  arma::mat out;

  if (!arma::inv_sympd(out, wishart)) {
    return scatter / df;
  }

  return out;
}

} // namespace

// ---------------------------------------------------------------------------
// A Dirichlet process mixture of normals for the error distribution, which is
// DPMBART (George, Laud, Logan, McCulloch and Sparapani 2019).
//
// BART assumes the errors are i.i.d. normal, and that assumption does most of
// the work in its uncertainty quantification. Here it is dropped: each
// observation gets its own error mean and variance,
//
//     y_i = f(x_i) + e_i,   e_i ~ N(mu_i, sigma_i^2),
//     theta_i = (mu_i, sigma_i) ~ G,   G ~ DP(alpha G_0),
//
// and because the Dirichlet process is discrete the theta_i take far fewer
// distinct values than there are observations. The error distribution is then
// whatever mixture of normals the data ask for -- heavy tailed, skewed,
// bimodal -- rather than the one normal BART is committed to.
//
// **Conditional on the theta_i the target is still exactly quadratic in the
// predictor**, which is what keeps this cheap: the leaf draw is the closed form,
// and the only cost over a Gaussian fit is the mixture update. The paper reports
// the total roughly doubling, which is what is measured here too.
//
// The baseline G_0 is the conjugate normal-inverse-chi-square,
//
//     sigma^2 ~ nu lambda / chisq_nu,    mu | sigma ~ N(mu_0, sigma^2 / k_0),
//
// and it has to be conjugate: the Escobar and West (1995) draws that make the
// mixture update a few lines are exactly the closed forms conjugacy provides.
// That is why this family does not use the half-Cauchy scale prior the Gaussian
// family here uses -- there is no conjugate mixture update to be had from it.
// ---------------------------------------------------------------------------

struct DPMFamily : Concrete<DPMFamily> {
  // Per-observation error mean and variance, which are what the predictor's
  // target sees.
  arma::vec shift;
  arma::vec scale2;

  // The mixture, held as a list of atoms plus a label per observation. Atoms are
  // deleted by swapping the last one into the hole, which is why the labels have
  // to be repaired on deletion.
  std::vector<double> atom_mu;
  std::vector<double> atom_s2;
  std::vector<int> atom_count;
  std::vector<int> label;

  // Baseline parameters, all fixed: the paper chooses values rather than priors
  // for these, on the argument that BART plus a Dirichlet process is already
  // adaptable enough that keeping the rest simple is worth more than another
  // layer.
  double nu;
  double lambda;
  double mu_0;
  double k_0;

  double alpha;
  bool update_alpha;
  arma::vec alpha_grid;
  arma::vec alpha_logprior;

  // Scratch, so that a sweep allocates nothing.
  mutable std::vector<double> weight_buffer;

  DPMFamily(const arma::vec& y_, const arma::vec& w_, double nu_,
            double lambda_, double mu_0_, double k_0_, double alpha_,
            bool update_alpha_, const arma::vec& alpha_grid_,
            const arma::vec& alpha_logprior_)
    : Concrete<DPMFamily>(y_, w_, 1), nu(nu_), lambda(lambda_), mu_0(mu_0_),
      k_0(k_0_), alpha(alpha_), update_alpha(update_alpha_),
      alpha_grid(alpha_grid_), alpha_logprior(alpha_logprior_) {

    shift = arma::vec(N, arma::fill::zeros);
    scale2 = arma::vec(N).fill(lambda);
    label.assign(N, 0);

    // One atom to begin with, at the baseline's centre. The first mixture update
    // splits it as the data ask.
    atom_mu.push_back(mu_0);
    atom_s2.push_back(lambda);
    atom_count.push_back(N);
    shift.fill(mu_0);
    refresh_eta_free();
  }

  // Conditional on the mixture the log density is a Gaussian one in the
  // predictor, with an offset and a variance that differ by observation.
  double logdens_unit(int i, const double* eta) const override {
    double r = y(i) - shift(i) - eta[0];
    return -0.5 * r * r / scale2(i);
  }

  arma::vec compute_eta_free() const override {
    arma::vec out(N);

    for (int i = 0; i < N; i++) {
      out(i) = -0.5 * (LN_2PI + std::log(scale2(i)));
    }

    return out;
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    *d1 = (y(i) - shift(i) - eta[0]) / scale2(i);
    *d2 = 1.0 / scale2(i);
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    return (y(i) - shift(i) - eta[0]) / scale2(i);
  }

  double info_unit(int i, const double* eta, int h) const override {
    return 1.0 / scale2(i);
  }

  TargetForm target_form(int h) const override {
    return TARGET_QUADRATIC;
  }

  // The marginal density of one residual under the baseline, which is the
  // weight a new atom gets in the Polya urn. Normal-inverse-chi-square
  // integrates to a t.
  double log_marginal(double r) const {
    double s = std::sqrt(lambda * (1.0 + 1.0 / k_0));
    return R::dt((r - mu_0) / s, nu, 1) - std::log(s);
  }

  // A draw from the baseline's posterior given `count` residuals with the given
  // sum and sum of squares. Also used with count = 1 for a fresh atom.
  void draw_atom(double count, double total, double total_sq, double* out_mu,
                 double* out_s2) const {
    double k_n = k_0 + count;
    double nu_n = nu + count;
    double mean = count > 0.0 ? total / count : 0.0;
    double centered = total_sq - count * mean * mean;
    double pull = k_0 * count / k_n * (mean - mu_0) * (mean - mu_0);
    double scale_n = (nu * lambda + centered + pull) / nu_n;
    double mu_n = (k_0 * mu_0 + total) / k_n;

    double drawn_s2 = nu_n * scale_n / Rf_rchisq(nu_n);

    if (!(drawn_s2 > 0.0) || !std::isfinite(drawn_s2)) {
      drawn_s2 = lambda;
    }

    *out_s2 = drawn_s2;
    *out_mu = mu_n + std::sqrt(drawn_s2 / k_n) * norm_rand();
  }

  void drop_atom(int at) {
    int last = static_cast<int>(atom_mu.size()) - 1;

    if (at != last) {
      atom_mu[at] = atom_mu[last];
      atom_s2[at] = atom_s2[last];
      atom_count[at] = atom_count[last];

      for (int i = 0; i < N; i++) {
        if (label[i] == last) {
          label[i] = at;
        }
      }
    }

    atom_mu.pop_back();
    atom_s2.pop_back();
    atom_count.pop_back();
  }

  // Escobar and West's sampler, in the two steps the paper takes from Dey,
  // Muller and Sinha (1998, sec. 1.3.3): reassign every observation, then
  // redraw every occupied atom given the observations that landed in it.
  void update_mixture(const arma::mat& eta) {
    const arma::rowvec& e = eta.row(0);

    for (int i = 0; i < N; i++) {
      double r = y(i) - e(i);

      // Take the observation out of its atom, and remove the atom if it empties.
      int mine = label[i];
      atom_count[mine]--;

      if (atom_count[mine] == 0) {
        drop_atom(mine);
      }

      int atoms = static_cast<int>(atom_mu.size());
      weight_buffer.resize(atoms + 1);
      double largest = R_NegInf;

      for (int k = 0; k < atoms; k++) {
        double resid = r - atom_mu[k];
        weight_buffer[k] = std::log(static_cast<double>(atom_count[k])) -
          0.5 * (LN_2PI + std::log(atom_s2[k]) + resid * resid / atom_s2[k]);

        if (weight_buffer[k] > largest) {
          largest = weight_buffer[k];
        }
      }

      weight_buffer[atoms] = std::log(alpha) + log_marginal(r);

      if (weight_buffer[atoms] > largest) {
        largest = weight_buffer[atoms];
      }

      double total = 0.0;

      for (int k = 0; k <= atoms; k++) {
        weight_buffer[k] = std::exp(weight_buffer[k] - largest);
        total += weight_buffer[k];
      }

      double u = unif_rand() * total;
      int chosen = atoms;
      double running = 0.0;

      for (int k = 0; k <= atoms; k++) {
        running += weight_buffer[k];

        if (u <= running) {
          chosen = k;
          break;
        }
      }

      if (chosen == atoms) {
        double fresh_mu;
        double fresh_s2;
        draw_atom(1.0, r, r * r, &fresh_mu, &fresh_s2);
        atom_mu.push_back(fresh_mu);
        atom_s2.push_back(fresh_s2);
        atom_count.push_back(1);
        label[i] = atoms;
      }
      else {
        atom_count[chosen]++;
        label[i] = chosen;
      }
    }

    // Redraw each atom given everything assigned to it.
    int atoms = static_cast<int>(atom_mu.size());
    std::vector<double> total(atoms, 0.0);
    std::vector<double> total_sq(atoms, 0.0);
    std::vector<double> count(atoms, 0.0);

    for (int i = 0; i < N; i++) {
      double r = y(i) - e(i);
      int k = label[i];
      count[k] += 1.0;
      total[k] += r;
      total_sq[k] += r * r;
    }

    for (int k = 0; k < atoms; k++) {
      draw_atom(count[k], total[k], total_sq[k], &atom_mu[k], &atom_s2[k]);
    }

    for (int i = 0; i < N; i++) {
      shift(i) = atom_mu[label[i]];
      scale2(i) = atom_s2[label[i]];
    }
  }

  // alpha on a grid. The number of occupied atoms is all the data say about it:
  // P(I = k | alpha) is proportional to alpha^k Gamma(alpha) / Gamma(alpha + n)
  // times a Stirling number that does not involve alpha and so cancels.
  void update_concentration() {
    if (!update_alpha || alpha_grid.n_elem == 0) {
      return;
    }

    double atoms = static_cast<double>(atom_mu.size());
    double n = static_cast<double>(N);
    arma::vec weights(alpha_grid.n_elem);
    double largest = R_NegInf;

    for (arma::uword g = 0; g < alpha_grid.n_elem; g++) {
      double a = alpha_grid(g);
      weights(g) = atoms * std::log(a) + R::lgammafn(a) -
        R::lgammafn(a + n) + alpha_logprior(g);

      if (weights(g) > largest) {
        largest = weights(g);
      }
    }

    weights = arma::exp(weights - largest);
    double total = arma::accu(weights);
    double u = unif_rand() * total;
    double running = 0.0;

    for (arma::uword g = 0; g < alpha_grid.n_elem; g++) {
      running += weights(g);

      if (u <= running) {
        alpha = alpha_grid(g);
        return;
      }
    }

    alpha = alpha_grid(alpha_grid.n_elem - 1);
  }

  void update_aux(const arma::mat& eta) override {
    update_mixture(eta);
    update_concentration();
    refresh_eta_free();
  }

  // The predictive density of one error under the current draw: the occupied
  // atoms weighted by their sizes, plus the chance under the Dirichlet process
  // that the next observation opens an atom of its own, which is the baseline's
  // own marginal.
  double log_predictive(double r) const {
    double n = static_cast<double>(N);
    double total = alpha + n;
    double out = std::exp(std::log(alpha) - std::log(total) + log_marginal(r));

    for (std::size_t k = 0; k < atom_mu.size(); k++) {
      double resid = r - atom_mu[k];
      out += atom_count[k] / total *
        std::exp(-0.5 * (LN_2PI + std::log(atom_s2[k]) +
                         resid * resid / atom_s2[k]));
    }

    return std::log(out);
  }

  // Reported as the mixture's own log likelihood rather than the complete-data
  // one, so that it is the quantity a Gaussian fit's log likelihood is, and so
  // that summing `predict(type = "density")` over observations reproduces it.
  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      out += w(i) * log_predictive(y(i) - e(i));
    }

    return out;
  }

  // The mean and standard deviation of the fitted error distribution, which is
  // the mixture over the atoms weighted by how many observations sit in each.
  void error_moments(double* out_mean, double* out_sd) const {
    double n = static_cast<double>(N);
    double mean = 0.0;
    double second = 0.0;

    for (std::size_t k = 0; k < atom_mu.size(); k++) {
      double p = static_cast<double>(atom_count[k]) / n;
      mean += p * atom_mu[k];
      second += p * (atom_s2[k] + atom_mu[k] * atom_mu[k]);
    }

    *out_mean = mean;
    *out_sd = std::sqrt(std::max(second - mean * mean, 0.0));
  }

  // Nothing forces the mixture to be centered, so the sampler works in a chart
  // where the sum of the predictor and the error mean is the conditional mean
  // and neither piece is identified on its own. Reporting is done in the chart
  // where the mixture has mean zero, which puts the whole conditional mean on
  // the predictor -- so `eta` means for `dpm()` what it means for `gaussian()`.
  // The shift is applied to the recorded predictor, the recorded leaf values and
  // the recorded mixture together, so every density is untouched.
  arma::vec report_shift(const arma::mat& eta) const override {
    double mean;
    double sd;
    error_moments(&mean, &sd);

    // model.cpp records `eta - shift`, and the predictor has to move *up* by the
    // error mean, so the shift is its negative.
    return arma::vec{-mean};
  }

  // `center` is the raw mixture's mean, which is exactly the shift taken out
  // above. It is reported because the predictive density's new-component term
  // needs it -- the baseline is centered on the *raw* chart -- and because it
  // says how far from symmetric the fitted error came out. It is a bookkeeping
  // quantity, not an estimate of anything: the error mean in the reported chart
  // is zero by construction.
  std::vector<std::string> aux_names() const override {
    return std::vector<std::string>{"alpha", "clusters", "center",
                                    "error_sd"};
  }

  arma::vec aux_values() const override {
    double mean;
    double sd;
    error_moments(&mean, &sd);

    return arma::vec{alpha, static_cast<double>(atom_mu.size()), mean, sd};
  }

  // Only alpha is recoverable from these four numbers; the mixture itself is
  // reported separately, as a flat vector of atoms, because it has a different
  // number of components at every draw.
  void set_aux(const arma::vec& values) override {
    if (values.n_elem > 0) {
      alpha = values(0);
    }
  }

  // The atoms of the current draw, as (mean, standard deviation, weight)
  // triples. `model.cpp` appends these to one flat vector across draws, the same
  // way it stores the trees.
  arma::vec mixture_flat() const override {
    arma::vec out(3 * atom_mu.size());
    double n = static_cast<double>(N);

    // Centered, to match the chart the predictor is reported in.
    double mean;
    double sd;
    error_moments(&mean, &sd);

    for (std::size_t k = 0; k < atom_mu.size(); k++) {
      out(3 * k) = atom_mu[k] - mean;
      out(3 * k + 1) = std::sqrt(atom_s2[k]);
      out(3 * k + 2) = static_cast<double>(atom_count[k]) / n;
    }

    return out;
  }
};

// ---------------------------------------------------------------------------
// The accelerated failure time model with a Dirichlet process mixture for its
// errors: Henderson, Louis, Rosner and Varadhan (2020).
//
//   log T_i = m(x_i) + W_i,   W_i ~ a mean-constrained DP mixture of normals,
//
// with right-censored log-times imputed. This is `dpm()`'s error model joined to
// the survival families' censoring, and both halves already existed: everything
// about the mixture -- the Polya urn, the atom draws, the concentration, the
// centering that makes the predictor the conditional mean of log T, the reported
// error density -- is inherited unchanged from `DPMFamily`. What is added is the
// imputation, and an observed-data likelihood that credits a censored
// observation with the mixture's survival rather than its density.
//
// The paper's error model is a location mixture with one common scale; this one
// is a location-scale mixture, because that is what `DPMFamily` already is, so
// the error distribution here is the more flexible of the two.
// ---------------------------------------------------------------------------

struct DPMAFTFamily : DPMFamily {
  arma::vec obs;      // log of the observed time, an event or a censoring
  arma::vec event;

  DPMAFTFamily(const arma::vec& y_, const arma::vec& w_,
               const arma::vec& event_, double nu_, double lambda_,
               double mu_0_, double k_0_, double alpha_, bool update_alpha_,
               const arma::vec& alpha_grid_, const arma::vec& alpha_logprior_)
    : DPMFamily(y_, w_, nu_, lambda_, mu_0_, k_0_, alpha_, update_alpha_,
                alpha_grid_, alpha_logprior_),
      obs(y_), event(event_) {
    // A deterministic start one scale unit past each censoring time, replaced by
    // a proper draw at the end of the first sweep.
    for (int i = 0; i < N; i++) {
      if (event(i) <= 0.0) {
        y(i) = obs(i) + std::sqrt(lambda);
      }
    }
  }

  static bool applies(const arma::vec& w_) {
    for (arma::uword i = 0; i < w_.n_elem; i++) {
      if (w_(i) != 1.0) {
        return false;
      }
    }
    return true;
  }

  void update_aux(const arma::mat& eta) override {
    const arma::rowvec& e = eta.row(0);

    // Impute each censored log-time from the component it currently sits in,
    // truncated below at its censoring time. Conditioning on the label is what
    // makes this an ordinary Gibbs step; the label itself is redrawn immediately
    // afterwards by the mixture update, given the value drawn here.
    for (int i = 0; i < N; i++) {
      if (event(i) > 0.0) {
        y(i) = obs(i);
        continue;
      }

      int k = label[i];
      double centre = e(i) + atom_mu[k];
      double sd = std::sqrt(atom_s2[k]);
      double lo = (obs(i) - centre) / sd;
      y(i) = centre + sd * truncated_normal_between(lo, R_PosInf);
    }

    DPMFamily::update_aux(eta);
  }

  // The survival function of one error under the current draw, matching
  // `log_predictive()` term for term: the occupied atoms weighted by their sizes
  // plus the baseline's own marginal, which is a t.
  double log_survival(double r) const {
    double n = static_cast<double>(N);
    double total = alpha + n;
    double scale = std::sqrt(lambda * (1.0 + 1.0 / k_0));
    double out = alpha / total *
      R::pt((r - mu_0) / scale, nu, 0, 0);

    for (std::size_t k = 0; k < atom_mu.size(); k++) {
      out += atom_count[k] / total *
        R::pnorm5(r, atom_mu[k], std::sqrt(atom_s2[k]), 0, 0);
    }

    if (!(out > 0.0)) {
      return R_NegInf;
    }

    return std::log(out);
  }

  // The *observed-data* likelihood: a density for an event, a survival
  // probability for a censoring. Not the complete-data one the imputation works
  // with, so that this is comparable with what the other survival families
  // report.
  double reported_loglik(const arma::mat& eta) const override {
    const arma::rowvec& e = eta.row(0);
    double out = 0.0;

    for (int i = 0; i < N; i++) {
      double r = obs(i) - e(i);
      out += w(i) * (event(i) > 0.0 ? log_predictive(r) : log_survival(r));
    }

    return out;
  }
};

// ---------------------------------------------------------------------------
// Multinomial probit (Imai and van Dyk 2005), with the sum-of-trees mean of
// Kindo, Wang and Pena (2016) and the sampler of Xu et al. (2025).
//
// The outcome is the argmax of C + 1 latent utilities. Differencing against the
// reference category leaves C latent variables per observation,
// W_i ~ MVN(eta_i, Sigma), with
//
//     S_i = l  if max(W_i) = W_il >= 0,     S_i = 0  if max(W_i) < 0,
//
// and one forest per component of eta. What the probit link buys over the
// logistic one is Sigma: the categories may be correlated, which a multinomial
// logit cannot express at all.
//
// Only the scale of W is unidentified, so Sigma is normalized by the trace
// constraint trace(Sigma) = C (Burgette and Nordheim 2012), which is what makes
// a two-category fit identical to binary probit.
//
// **Conditional on W and Sigma the target is exactly quadratic in every
// component of eta**, which is the whole reason this is fast: the leaf draw is
// the closed form rather than a Laplace approximation, exactly as for the
// augmented ordinal probit. The score and information in component h are read
// off the precision matrix P = Sigma^{-1}:
//
//     d/deta_h = sum_k P_hk (W_k - eta_k),    -d^2/deta_h^2 = P_hh.
//
// The sampler is Algorithm [P2] of Xu et al. (2025): draw W by a Gibbs sweep of
// truncated normals, fit the forests to it, then draw the unnormalized
// covariance from its inverse Wishart conditional and rescale to the trace
// constraint. Their Algorithms [P1] and [P2] measured indistinguishable on every
// figure in that paper, and [P2] is the one with no expansion parameter at all.
// Both beat the earlier [KD] sampler, which fits the trees to the *unnormalized*
// utilities and so grows them much deeper -- average depths of 6 and 9 against
// about 2 -- because the quantity the stochastic search is chasing keeps being
// rescaled underneath it.
//
// There is no closed form for P(S_i = k | eta_i, Sigma): it is a C-dimensional
// Gaussian orthant probability. The reported log likelihood is therefore
// simulated, with a fixed set of standard normal draws held by the family so
// that it is a deterministic function of eta and Sigma and adds no Monte Carlo
// noise to the chain. Predictions simulate too, with fresh draws.
// ---------------------------------------------------------------------------

struct MultinomProbitFamily : Concrete<MultinomProbitFamily> {
  arma::vec category;      // observed category, 0 for the reference level
  int num_cat;
  arma::mat latent;        // W, C by N
  arma::mat sigma;         // normalized, trace(sigma) = C
  arma::mat prec;          // sigma inverse
  arma::mat chol_sigma;    // lower Cholesky factor of sigma
  arma::mat unit_draws;    // fixed standard normals, C by replicates
  double nu;               // inverse Wishart degrees of freedom
  bool update_sigma;
  int replicates;

  MultinomProbitFamily(const arma::vec& y_, const arma::vec& w_, int num_cat_,
                       double nu_, bool update_sigma_, int replicates_)
    : Concrete<MultinomProbitFamily>(y_, w_, num_cat_ - 1),
      category(y_), num_cat(num_cat_), nu(nu_),
      update_sigma(update_sigma_), replicates(replicates_) {
    int C = H;
    sigma = arma::eye<arma::mat>(C, C);
    prec = arma::eye<arma::mat>(C, C);
    chol_sigma = arma::eye<arma::mat>(C, C);

    // A start that satisfies the sign constraints, so that the first Gibbs
    // sweep over the latent variables has something coherent to condition on.
    latent = arma::mat(C, N);
    latent.fill(-1.0);

    for (int i = 0; i < N; i++) {
      int k = static_cast<int>(category(i));

      if (k > 0) {
        latent(k - 1, i) = 1.0;
      }
    }

    unit_draws = arma::mat(C, replicates);

    for (int r = 0; r < replicates; r++) {
      for (int l = 0; l < C; l++) {
        unit_draws(l, r) = norm_rand();
      }
    }
  }

  TargetForm target_form(int h) const override {
    return TARGET_QUADRATIC;
  }

  double logdens_unit(int i, const double* eta) const override {
    double out = 0.0;

    for (int l = 0; l < H; l++) {
      double r_l = latent(l, i) - eta[l];
      out -= 0.5 * prec(l, l) * r_l * r_l;

      for (int k = l + 1; k < H; k++) {
        out -= prec(l, k) * r_l * (latent(k, i) - eta[k]);
      }
    }

    return out;
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    double score = 0.0;

    for (int k = 0; k < H; k++) {
      score += prec(h, k) * (latent(k, i) - eta[k]);
    }

    *d1 = score;
    *d2 = prec(h, h);
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    double score = 0.0;

    for (int k = 0; k < H; k++) {
      score += prec(h, k) * (latent(k, i) - eta[k]);
    }

    return score;
  }

  double info_unit(int i, const double* eta, int h) const override {
    return prec(h, h);
  }

  // Step 1 of Algorithm [P2]: a Gibbs sweep over the latent variables, each
  // from its univariate conditional truncated to the region the observed
  // category defines. The conditional moments come from the precision matrix
  // rather than from a sweep operator, which makes each one O(C) instead of
  // O(C^3).
  void draw_latent(const arma::mat& eta) {
    for (int i = 0; i < N; i++) {
      int winner = static_cast<int>(category(i)) - 1;   // -1 for the reference

      for (int l = 0; l < H; l++) {
        double variance = 1.0 / (w(i) * prec(l, l));
        double centered = 0.0;

        for (int k = 0; k < H; k++) {
          if (k != l) {
            centered += prec(l, k) * (latent(k, i) - eta(k, i));
          }
        }

        double mean = eta(l, i) - centered / prec(l, l);
        double lower = R_NegInf;
        double upper = R_PosInf;

        if (winner < 0) {
          // Every utility is below the reference category's.
          upper = 0.0;
        }
        else if (l == winner) {
          // The winner is above zero and above every rival.
          lower = 0.0;

          for (int k = 0; k < H; k++) {
            if (k != l && latent(k, i) > lower) {
              lower = latent(k, i);
            }
          }
        }
        else {
          upper = latent(winner, i);
        }

        double scale = std::sqrt(variance);
        double draw = truncated_normal_between((lower - mean) / scale,
                                               (upper - mean) / scale);
        latent(l, i) = mean + scale * draw;
      }
    }
  }

  // Step 3: the unnormalized covariance from its inverse Wishart conditional,
  // then the rescaling that puts it back on the trace constraint. The scale
  // factor is applied to the latent variables as well, which is what makes this
  // a move in the normalized space rather than a reparameterization of it.
  void update_aux(const arma::mat& eta) override {
    if (!update_sigma || H < 2) {
      return;
    }

    int C = H;
    arma::mat scatter = arma::eye<arma::mat>(C, C);   // Psi, the identity
    double df = nu;

    for (int i = 0; i < N; i++) {
      arma::vec resid(C);

      for (int l = 0; l < C; l++) {
        resid(l) = latent(l, i) - eta(l, i);
      }

      scatter += w(i) * resid * resid.t();
      df += w(i);
    }

    arma::mat drawn = inverse_wishart(df, scatter);
    double factor = arma::trace(drawn) / static_cast<double>(C);

    if (!(factor > 0.0) || !std::isfinite(factor)) {
      return;
    }

    sigma = drawn / factor;
    double root = std::sqrt(factor);

    for (int i = 0; i < N; i++) {
      for (int l = 0; l < C; l++) {
        latent(l, i) = eta(l, i) + (latent(l, i) - eta(l, i)) / root;
      }
    }

    refresh_sigma();
  }

  void refresh_sigma() {
    if (!arma::inv_sympd(prec, sigma)) {
      prec = arma::inv(sigma);
    }
    if (!arma::chol(chol_sigma, sigma, "lower")) {
      chol_sigma = arma::eye<arma::mat>(H, H);
    }
  }

  void before_forest(int h, const arma::mat& eta) override {
    // Once a sweep, before the first forest moves, as in Algorithm [P2].
    if (h == 0) {
      draw_latent(eta);
    }
  }

  // With two categories the trace constraint pins the single variance at one, so
  // there is nothing to record.
  std::vector<std::string> aux_names() const override {
    std::vector<std::string> out;

    if (H < 2) {
      return out;
    }

    for (int l = 0; l < H; l++) {
      for (int k = 0; k <= l; k++) {
        out.push_back("sigma" + std::to_string(l + 1) + std::to_string(k + 1));
      }
    }

    return out;
  }

  arma::vec aux_values() const override {
    if (H < 2) {
      return arma::vec();
    }

    arma::vec out(H * (H + 1) / 2);
    int at = 0;

    for (int l = 0; l < H; l++) {
      for (int k = 0; k <= l; k++) {
        out(at++) = sigma(l, k);
      }
    }

    return out;
  }

  void set_aux(const arma::vec& values) override {
    if (H < 2 || static_cast<int>(values.n_elem) != H * (H + 1) / 2) {
      return;
    }

    int at = 0;

    for (int l = 0; l < H; l++) {
      for (int k = 0; k <= l; k++) {
        sigma(l, k) = values(at);
        sigma(k, l) = values(at);
        at++;
      }
    }

    refresh_sigma();
  }

  // The likelihood the caller asked for, which is a Gaussian orthant
  // probability and so is simulated. The draws are fixed at construction, so
  // this is a deterministic function of eta and Sigma: the chain sees no Monte
  // Carlo noise, only a bias that is the same at every iteration.
  double reported_loglik(const arma::mat& eta) const override {
    double out = 0.0;
    arma::vec value(H);

    for (int i = 0; i < N; i++) {
      int winner = static_cast<int>(category(i)) - 1;
      int hits = 0;

      for (int r = 0; r < replicates; r++) {
        for (int l = 0; l < H; l++) {
          double shift = 0.0;

          for (int k = 0; k <= l; k++) {
            shift += chol_sigma(l, k) * unit_draws(k, r);
          }

          value(l) = eta(l, i) + shift;
        }

        if (chosen(value) == winner) {
          hits++;
        }
      }

      double p = (hits > 0 ? static_cast<double>(hits) : 0.5) /
        static_cast<double>(replicates);
      out += w(i) * std::log(p);
    }

    return out;
  }

  // Which category a vector of latent utilities implies: the largest if it is
  // non-negative, and the reference otherwise.
  static int chosen(const arma::vec& value) {
    int at = -1;
    double best = 0.0;

    for (arma::uword l = 0; l < value.n_elem; l++) {
      if (value(l) >= best) {
        best = value(l);
        at = static_cast<int>(l);
      }
    }

    return at;
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
    Rcpp::stop("the %s supplied to bartisan() must return a numeric vector.",
               what);
  }

  arma::vec out = Rcpp::as<arma::vec>(value);

  if (static_cast<int>(out.n_elem) != n) {
    Rcpp::stop("the %s supplied to bartisan() returned %d values for %d "
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

  // A nuisance parameter is carried as an additive predictor whose forest is
  // pinned at depth zero -- one tree that can never split -- so it is a single
  // scalar drawn by the same Laplace-plus-Metropolis step as any leaf, with no
  // separate sampler. The engine does the pinning; all this class does is keep
  // the two kinds apart at the boundary with R. The caller's function sees
  // `eta` with `num_predictors` columns and `aux` as a plain numeric vector, and
  // never learns that the second kind is a forest.
  int num_predictors;
  int num_aux;
  std::vector<std::string> aux_labels;
  arma::vec aux_cache;

  RFamily(const arma::vec& y_, const arma::vec& w_, int num_predictors_,
          int num_aux_, const std::vector<std::string>& aux_labels_,
          const Rcpp::Function& dens_, const Rcpp::RObject& derivs_,
          const std::string& label_)
    : Family(y_, w_, num_predictors_ + num_aux_), dens(dens_),
      derivs(Rf_isFunction(derivs_) ? Rcpp::Function(derivs_) : dens_),
      has_derivs(Rf_isFunction(derivs_)), label(label_),
      num_predictors(num_predictors_), num_aux(num_aux_),
      aux_labels(aux_labels_),
      aux_cache(std::max(num_aux_, 1), arma::fill::zeros) {}

  bool wants_block() const override { return true; }

  // The response and the predictors of one block, in the shape the caller's
  // function expects: y a vector of length n, eta an n by num_predictors
  // matrix, and the nuisance parameters as a vector of length num_aux.
  //
  // The nuisance columns are constant down the block, because a pinned forest
  // has one leaf holding every observation and so shifts them all together, so
  // reading the first row is reading the parameter.
  void unpack(const int* idx, int n, const double* block, int from,
              Rcpp::NumericVector& y_out, Rcpp::NumericMatrix& eta_out,
              Rcpp::NumericVector& aux_out) const {
    for (int k = 0; k < n; k++) {
      y_out[k] = y(idx[from + k]);

      for (int j = 0; j < num_predictors; j++) {
        eta_out(k, j) = block[static_cast<std::size_t>(from + k) * H + j];
      }
    }

    for (int j = 0; j < num_aux; j++) {
      aux_out[j] =
        block[static_cast<std::size_t>(from) * H + num_predictors + j];
    }
  }

  // A nuisance parameter is one number, so the caller is handed a vector rather
  // than a column -- but the paired evaluations the sampler uses stack two values
  // of one component in a single block, so the column is constant only in runs.
  // These find the runs. Every call that is not a paired one has exactly one,
  // and a family with no nuisance parameters always has exactly one, so the
  // common path is unchanged.
  bool same_aux(const double* block, int a, int b) const {
    for (int j = 0; j < num_aux; j++) {
      if (block[static_cast<std::size_t>(a) * H + num_predictors + j] !=
          block[static_cast<std::size_t>(b) * H + num_predictors + j]) {
        return false;
      }
    }
    return true;
  }

  int aux_run_end(const double* block, int from, int n) const {
    if (num_aux == 0) {
      return n;
    }

    int stop = from + 1;
    while (stop < n && same_aux(block, from, stop)) {
      stop++;
    }
    return stop;
  }

  arma::vec evaluate(const Rcpp::NumericVector& y_in,
                     const Rcpp::NumericMatrix& eta_in,
                     const Rcpp::NumericVector& aux_in, int n) const {
    Rcpp::RObject value = num_aux > 0 ? dens(y_in, eta_in, aux_in)
                                      : dens(y_in, eta_in);

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
    int from = 0;

    while (from < n) {
      int stop = aux_run_end(block, from, n);
      int m = stop - from;
      Rcpp::NumericVector yv(m);
      Rcpp::NumericMatrix ev(m, num_predictors);
      Rcpp::NumericVector av(std::max(num_aux, 1));
      unpack(idx, m, block, from, yv, ev, av);
      arma::vec value = evaluate(yv, ev, av, m);

      for (int k = 0; k < m; k++) {
        out[from + k] = w(idx[from + k]) * value(k);
      }

      from = stop;
    }
  }

  void score_info_block(const int* idx, int n, const double* block, int h,
                        double* d1, double* d2) const override {
    int from = 0;

    while (from < n) {
      int stop = aux_run_end(block, from, n);
      score_info_run(idx, from, stop - from, block, h, d1, d2);
      from = stop;
    }
  }

  void score_info_run(const int* idx, int from, int n, const double* block,
                      int h, double* d1_all, double* d2_all) const {
    double* d1 = d1_all + from;
    double* d2 = d2_all + from;
    Rcpp::NumericVector yv(n);
    Rcpp::NumericMatrix ev(n, num_predictors);
    Rcpp::NumericVector av(std::max(num_aux, 1));
    unpack(idx, n, block, from, yv, ev, av);

    // Supplied derivatives cover the additive predictors only. A nuisance
    // parameter is always differenced, which is affordable in a way that
    // differencing a predictor is not: its forest has one leaf, so it costs
    // three calls per sweep rather than three per leaf visit, and asking a
    // caller to hand-derive it would be a poor trade.
    if (has_derivs && h < num_predictors) {
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
        d1[k] = w(idx[from + k]) * score(k);
        d2[k] = clamp_info(w(idx[from + k]) * info(k));
      }

      return;
    }

    const double step = 1e-4;
    arma::vec mid = evaluate(yv, ev, av, n);

    // Differencing moves whichever of the two the index names.
    bool is_aux = h >= num_predictors;
    int col = is_aux ? h - num_predictors : h;

    if (is_aux) {
      av[col] += step;
    }
    else {
      for (int k = 0; k < n; k++) {
        ev(k, col) += step;
      }
    }

    arma::vec up = evaluate(yv, ev, av, n);

    if (is_aux) {
      av[col] -= 2.0 * step;
    }
    else {
      for (int k = 0; k < n; k++) {
        ev(k, col) -= 2.0 * step;
      }
    }

    arma::vec down = evaluate(yv, ev, av, n);

    for (int k = 0; k < n; k++) {
      d1[k] = w(idx[from + k]) * 0.5 * (up(k) - down(k)) / step;
      d2[k] = clamp_info(-w(idx[from + k]) *
        (up(k) - 2.0 * mid(k) + down(k)) / (step * step));
    }
  }

  double logdens_unit(int i, const double* eta) const override {
    Rcpp::NumericVector yv(1);
    Rcpp::NumericMatrix ev(1, num_predictors);
    Rcpp::NumericVector av(std::max(num_aux, 1));
    yv[0] = y(i);

    for (int j = 0; j < num_predictors; j++) {
      ev(0, j) = eta[j];
    }
    for (int j = 0; j < num_aux; j++) {
      av[j] = eta[num_predictors + j];
    }

    return evaluate(yv, ev, av, 1)(0);
  }

  // The nuisance parameters are reported as parameters rather than as the
  // predictors they are carried as, so they land in `fit$aux` under their own
  // names and `fit$eta` holds only what the caller asked for.
  int num_pinned() const override { return num_aux; }

  std::vector<std::string> aux_names() const override { return aux_labels; }

  arma::vec aux_values() const override {
    return num_aux > 0 ? arma::vec(aux_cache.head(num_aux))
                       : arma::vec(0, arma::fill::zeros);
  }

  void set_aux(const arma::vec& values) override {
    for (int j = 0; j < num_aux && j < static_cast<int>(values.n_elem); j++) {
      aux_cache(j) = values(j);
    }
  }

  // Called once at the end of every sweep, which is where the drawn values are
  // read off the predictors they live in and cached for reporting.
  void update_aux(const arma::mat& eta) override {
    for (int j = 0; j < num_aux; j++) {
      aux_cache(j) = eta(num_predictors + j, 0);
    }
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

// A varying-coefficient model, as a decorator.
//
//   g(mu_i) = eta_i0 + sum_j basis(i, j) * eta_i(j+1)
//
// The wrapped family never learns this is happening: it is handed one predictor,
// mu, exactly as it would be in a single-forest model, and hands back a score
// and an information with respect to it. The chain rule is the whole of what
// this class adds -- d mu / d eta_h is 1 for the control function and
// basis(i, h - 1) for a coefficient -- so *every* family gets varying
// coefficients with no per-family code, which is the reason to write it this
// way rather than as a family of its own.
//
// This is Deshpande, Bai, Balocchi, Starling and Weiss (2026), and the case of
// one binary column is the Bayesian causal forest of Hahn, Murray and Carvalho
// (2020).
struct VaryingCoefficientFamily : Family {
  std::unique_ptr<Family> inner;
  arma::mat basis;

  // Data-adaptive coding, Hahn, Murray and Carvalho (2020) sec. 5.3. Rather than
  // fix what the covariate's levels are coded as, draw them: the model becomes
  // mu = eta_0 + b_{z_i} * eta_1 with b_k ~ N(0, 1/2) independently, so every
  // pairwise contrast b_j - b_k is marginally N(0, 1) and no level is a
  // reference. `coding(i, j)` is observation i's level in term j, or -1 when
  // that term's coding is fixed; `coding_levels(j)` is how many levels it has.
  arma::imat coding;
  arma::ivec coding_levels;
  std::vector<arma::vec> b;
  std::vector<std::string> b_labels;

  // One entry per forest: which of the wrapped family's additive predictors it
  // feeds, and which basis column it is multiplied by (-1 for a control
  // function, and for a pinned nuisance forest passing straight through).
  arma::ivec param;
  arma::ivec column;

  VaryingCoefficientFamily(Family* inner_, const arma::mat& basis_,
                           const arma::imat& coding_,
                           const arma::ivec& coding_levels_,
                           const std::vector<std::string>& b_labels_,
                           const arma::ivec& param_,
                           const arma::ivec& column_)
    : Family(inner_->y, inner_->w, static_cast<int>(param_.n_elem)),
      inner(inner_), basis(basis_), coding(coding_),
      coding_levels(coding_levels_), b_labels(b_labels_),
      param(param_), column(column_) {

    if (coding_levels.n_elem == 0) {
      return;
    }

    // Initialized at evenly spaced quantiles of the prior rather than drawn, so
    // the levels start apart: with every b_k equal the coefficient forest is
    // multiplied by a constant and has nothing to identify it.
    for (arma::uword j = 0; j < coding_levels.n_elem; j++) {
      int levels = coding_levels(j);
      arma::vec start(std::max(levels, 1), arma::fill::zeros);

      for (int k = 0; k < levels; k++) {
        start(k) = R::qnorm((k + 0.5) / levels, 0.0, 1.0, 1, 0) / std::sqrt(2.0);
      }

      b.push_back(start);

      if (levels > 0) {
        refresh_coding_column(static_cast<int>(j));
      }
    }
  }

  bool has_coding() const {
    return arma::any(coding_levels > 0);
  }

  // Exact only where the wrapped family's leaf target is quadratic; see
  // `Family::coding_is_exact()`. Asked of the predictor the drawn coding
  // actually feeds, since a family can be quadratic in one and not another.
  // With no drawn coding there is nothing to be exact about, so every fixed
  // centring passes.
  bool coding_is_exact() const override {
    return coding_not_exact() < 0;
  }

  // Which forest's drawn coding is not exact, or -1 if every one is. The index
  // rather than a flag, so the refusal can name the coefficient at fault: with
  // several additive predictors the family may be quadratic in one and not
  // another, and `location_scale()` is exactly that -- a drawn coding is fine on
  // the mean and not on the log standard deviation.
  int coding_not_exact() const override {
    for (arma::uword j = 0; j < coding_levels.n_elem; j++) {
      if (coding_levels(j) > 0 &&
          inner->target_form(param_of_column(static_cast<int>(j))) !=
            TARGET_QUADRATIC) {
        return forest_of_column(static_cast<int>(j));
      }
    }

    return -1;
  }

  // The forest a basis column belongs to, and through it the predictor it feeds.
  // Each column is multiplied by exactly one forest, so the search is a lookup.
  int param_of_column(int j) const {
    for (arma::uword h = 0; h < column.n_elem; h++) {
      if (column(h) == j) {
        return param(h);
      }
    }

    return 0;
  }

  int forest_of_column(int j) const {
    for (arma::uword h = 0; h < column.n_elem; h++) {
      if (column(h) == j) {
        return static_cast<int>(h);
      }
    }

    return 0;
  }

  void refresh_coding_column(int j) {
    for (int i = 0; i < N; i++) {
      int level = coding(i, j);
      basis(i, j) = level < 0 ? 0.0 : b[j](level);
    }
  }

  bool wants_block() const override { return inner->wants_block(); }

  // The wrapped family's trailing nuisance forests stay pinned. They are not
  // additive predictors and carry no coefficient, so the wrapper leaves them
  // exactly where they were -- at the end -- and has to say so, or the engine
  // treats each as an ordinary forest and the nuisance parameter is never drawn.
  int num_pinned() const override { return inner->num_pinned(); }

  // d mu_p / d eta_h, where p is the predictor forest h feeds. A control function
  // is not multiplied by anything and a coefficient is multiplied by its own
  // column; every other predictor does not involve eta_h at all, which is what
  // makes the chain rule below one factor rather than a sum.
  double slope(int i, int h) const {
    int j = column(h);
    return j < 0 ? 1.0 : basis(i, j);
  }

  // The wrapped family's predictors at observation i, from this family's. Each
  // one is its control function plus every coefficient of it times its column.
  void combine_all(int i, const double* eta, double* mu) const {
    for (int p = 0; p < inner->H; p++) {
      mu[p] = 0.0;
    }

    for (arma::uword h = 0; h < param.n_elem; h++) {
      mu[param(h)] += slope(i, static_cast<int>(h)) * eta[h];
    }
  }

  double logdens_unit(int i, const double* eta) const override {
    std::vector<double> mu(inner->H);
    combine_all(i, eta, mu.data());
    return inner->logdens_unit(i, mu.data());
  }

  void score_info_unit(int i, const double* eta, int h, double* d1,
                       double* d2) const override {
    std::vector<double> mu(inner->H);
    combine_all(i, eta, mu.data());
    inner->score_info_unit(i, mu.data(), param(h), d1, d2);

    double s = slope(i, h);
    *d1 *= s;
    *d2 *= s * s;
  }

  double dlogdens_unit(int i, const double* eta, int h) const override {
    std::vector<double> mu(inner->H);
    combine_all(i, eta, mu.data());
    return inner->dlogdens_unit(i, mu.data(), param(h)) * slope(i, h);
  }

  double info_unit(int i, const double* eta, int h) const override {
    std::vector<double> mu(inner->H);
    combine_all(i, eta, mu.data());
    double s = slope(i, h);
    return inner->info_unit(i, mu.data(), param(h)) * s * s;
  }

  // The block form gets the same treatment. `block` is H doubles per
  // observation, laid out the way the sampler stores eta; the inner family wants
  // its own H per observation, in the same layout.
  void combine_block(const int* idx, int n, const double* block,
                     std::vector<double>* mu) const {
    mu->assign(static_cast<std::size_t>(n) * inner->H, 0.0);

    for (int k = 0; k < n; k++) {
      combine_all(idx[k], block + static_cast<std::size_t>(k) * H,
                  mu->data() + static_cast<std::size_t>(k) * inner->H);
    }
  }

  void logdens_block(const int* idx, int n, const double* block,
                     double* out) const override {
    std::vector<double> mu;
    combine_block(idx, n, block, &mu);

    inner->logdens_block(idx, n, mu.data(), out);
  }

  void score_info_block(const int* idx, int n, const double* block, int h,
                        double* d1, double* d2) const override {
    std::vector<double> mu;
    combine_block(idx, n, block, &mu);

    inner->score_info_block(idx, n, mu.data(), param(h), d1, d2);

    for (int k = 0; k < n; k++) {
      double s = slope(idx[k], h);
      d1[k] *= s;
      d2[k] *= s * s;
    }
  }

  // The map is linear in each eta_h, so a target that is quadratic in the
  // predictor this forest feeds is quadratic in eta_h and the closed-form leaf
  // draw survives. The exponential form does not: `exp_rate()` is one scalar per
  // predictor and a varying coefficient makes the rate its own basis column,
  // which varies by observation. Those families fall back to the general path.
  TargetForm target_form(int h) const override {
    if (inner->target_form(param(h)) == TARGET_QUADRATIC) {
      return TARGET_QUADRATIC;
    }

    return TARGET_GENERAL;
  }

  double log_norm_const(int i) const override {
    return inner->log_norm_const(i);
  }

  double logdens_extra_total() const override {
    return inner->logdens_extra_total();
  }

  arma::vec compute_eta_free() const override {
    return inner->eta_free_part();
  }

  // The nuisance parameters belong to the wrapped family and are drawn on its
  // scale, so the predictors are combined before they are handed over. The
  // coding coefficients are drawn here too, since this is the once-a-sweep hook
  // and they are a nuisance parameter of exactly the same kind.
  void update_aux(const arma::mat& eta) override {
    draw_coding(eta);

    arma::mat mu(inner->H, eta.n_cols);

    for (arma::uword i = 0; i < eta.n_cols; i++) {
      combine_all(static_cast<int>(i), eta.colptr(i), mu.colptr(i));
    }

    inner->update_aux(mu);
    refresh_eta_free();
  }

  // One conjugate normal draw per level.
  //
  // Only the observations at level k carry b_k, so the design is block diagonal
  // and the levels are independent given everything else -- no matrix to build
  // or invert. Writing the target as a weighted least squares problem in mu,
  // which is what a quadratic target is, the level's precision is its prior
  // precision of 2 plus the sum of d2 * eta_1^2 over its own rows, and its mean
  // follows from the score at the current value.
  //
  // Terms are drawn one at a time with mu recomputed between them, since mu
  // depends on all of them.
  void draw_coding(const arma::mat& eta) {
    if (!has_coding()) {
      return;
    }

    for (arma::uword j = 0; j < coding_levels.n_elem; j++) {
      int levels = coding_levels(j);

      if (levels <= 0) {
        continue;
      }

      int h = forest_of_column(static_cast<int>(j));
      int p = param(h);
      arma::vec info(levels, arma::fill::zeros);
      arma::vec score(levels, arma::fill::zeros);
      std::vector<double> mu(inner->H);

      for (int i = 0; i < N; i++) {
        int level = coding(i, static_cast<int>(j));

        if (level < 0) {
          continue;
        }

        combine_all(i, eta.colptr(i), mu.data());
        double d1;
        double d2;
        inner->score_info(i, mu.data(), p, &d1, &d2);

        double slope_i = eta(h, i);
        info(level) += d2 * slope_i * slope_i;
        score(level) += d1 * slope_i;
      }

      for (int k = 0; k < levels; k++) {
        // The prior is N(0, 1/2), so its precision is 2 and every pairwise
        // contrast is marginally N(0, 1).
        double prec = info(k) + 2.0;
        double mean = (b[j](k) * info(k) + score(k)) / prec;
        b[j](k) = mean + norm_rand() / std::sqrt(prec);
      }

      refresh_coding_column(static_cast<int>(j));
    }
  }

  std::vector<std::string> aux_names() const override {
    std::vector<std::string> out = inner->aux_names();
    out.insert(out.end(), b_labels.begin(), b_labels.end());
    return out;
  }

  arma::vec aux_values() const override {
    arma::vec inner_values = inner->aux_values();

    if (b_labels.empty()) {
      return inner_values;
    }

    arma::vec out(inner_values.n_elem + b_labels.size());
    out.head(inner_values.n_elem) = inner_values;

    arma::uword at = inner_values.n_elem;

    for (arma::uword j = 0; j < b.size(); j++) {
      for (arma::uword k = 0; k < b[j].n_elem; k++) {
        if (coding_levels(j) > 0) {
          out(at++) = b[j](k);
        }
      }
    }

    return out;
  }

  void set_aux(const arma::vec& values) override {
    inner->set_aux(values);
    refresh_eta_free();
  }

  arma::vec aux_values_shifted(const arma::vec& shift) const override {
    if (b_labels.empty()) {
      return inner->aux_values_shifted(shift);
    }

    return aux_values();
  }

  arma::vec mixture_flat() const override { return inner->mixture_flat(); }
};

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

namespace {

Family* make_base_family(const std::string& name, const std::string& link,
                         const arma::vec& y, const arma::vec& w,
                         const List& opts) {

  if (name == "custom") {
    // The nuisance-parameter fields are optional, so that opts assembled by hand
    // -- as the density and derivative entry points are called in tests -- still
    // describe a family with no nuisance parameters.
    int n_aux = opts.containsElementNamed("num_aux")
      ? as<int>(opts["num_aux"]) : 0;
    std::vector<std::string> aux_labels;
    if (opts.containsElementNamed("aux_names")) {
      aux_labels = as<std::vector<std::string> >(opts["aux_names"]);
    }

    return finish(new RFamily(y, w, as<int>(opts["num_predictors"]), n_aux,
                              aux_labels,
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

  if (name == "dpm_aft") {
    return finish(new DPMAFTFamily(y, w, as<arma::vec>(opts["event"]),
                                   as<double>(opts["nu"]),
                                   as<double>(opts["lambda"]),
                                   as<double>(opts["mu_0"]),
                                   as<double>(opts["k_0"]),
                                   as<double>(opts["alpha"]),
                                   as<bool>(opts["update_alpha"]),
                                   as<arma::vec>(opts["alpha_grid"]),
                                   as<arma::vec>(opts["alpha_logprior"])));
  }

  if (name == "ph") {
    return finish(new PHFamily(y, w, as<arma::vec>(opts["event"]),
                               as<arma::vec>(opts["edges"]),
                               as<double>(opts["lambda_shape"]),
                               as<double>(opts["lambda_rate"]),
                               as<bool>(opts["update_lambda"])));
  }

  if (name == "location_scale") {
    return finish(new LocationScaleFamily(y, w));
  }

  if (name == "dpm") {
    return finish(new DPMFamily(y, w, as<double>(opts["nu"]),
                                as<double>(opts["lambda"]),
                                as<double>(opts["mu_0"]),
                                as<double>(opts["k_0"]),
                                as<double>(opts["alpha"]),
                                as<bool>(opts["update_alpha"]),
                                as<arma::vec>(opts["alpha_grid"]),
                                as<arma::vec>(opts["alpha_logprior"])));
  }

  if (name == "mnp") {
    return finish(new MultinomProbitFamily(y, w, as<int>(opts["num_cat"]),
                                           as<double>(opts["nu"]),
                                           as<bool>(opts["update_sigma"]),
                                           as<int>(opts["replicates"])));
  }

  if (name == "zip" || name == "zinb") {
    bool nb = name == "zinb";
    return finish(new ZeroInflatedFamily(y, w, nb,
                                  nb ? as<double>(opts["theta"]) : 1.0,
                                  nb ? as<double>(opts["theta_prior_shape"]) : 0.0,
                                  nb ? as<double>(opts["theta_prior_rate"]) : 0.0,
                                  nb ? as<bool>(opts["update_theta"]) : false));
  }

  if (name == "beta") {
    return finish(new BetaFamily(y, w, as<double>(opts["phi"]),
                          as<double>(opts["phi_prior_shape"]),
                          as<double>(opts["phi_prior_rate"]),
                          as<bool>(opts["update_phi"])));
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

} // namespace

// The family the sampler sees. `vc_basis` empty is the ordinary case and the
// base family is returned untouched, so a model with no varying coefficient
// reaches the engine exactly as it did before.
// Everything a varying-coefficient model needs beyond the basis, kept together
// so the two places that build one stay in step.
struct VaryingSpec {
  arma::mat basis;
  arma::imat coding;
  arma::ivec coding_levels;
  std::vector<std::string> b_labels;

  // One entry per forest the sampler builds. `param` is the wrapped family's
  // additive predictor it feeds, zero-based; `column` is the basis column it is
  // multiplied by, or -1 for a control function -- and for a custom family's
  // pinned nuisance forests, which pass straight through.
  arma::ivec param;
  arma::ivec column;

  bool adaptive() const {
    return coding_levels.n_elem > 0 && arma::any(coding_levels > 0);
  }
};

VaryingSpec varying_spec(const List& opts, const arma::mat& basis) {
  VaryingSpec out;
  out.basis = basis;
  out.coding = arma::imat(basis.n_rows, basis.n_cols);
  out.coding.fill(-1);
  out.coding_levels = arma::ivec(basis.n_cols, arma::fill::zeros);

  if (opts.containsElementNamed("vc_coding") && !Rf_isNull(opts["vc_coding"])) {
    out.coding = as<arma::imat>(opts["vc_coding"]);
    out.coding_levels = as<arma::ivec>(opts["vc_coding_levels"]);
    out.b_labels = as<std::vector<std::string> >(opts["vc_coding_names"]);
  }

  // Absent only from a path that predates several additive predictors, where the
  // map is the one this reconstructs: a control function followed by every
  // coefficient, all feeding the family's single predictor.
  if (opts.containsElementNamed("vc_param") && !Rf_isNull(opts["vc_param"])) {
    out.param = as<arma::ivec>(opts["vc_param"]);
    out.column = as<arma::ivec>(opts["vc_column"]);
  }
  else {
    out.param = arma::ivec(basis.n_cols + 1, arma::fill::zeros);
    out.column = arma::regspace<arma::ivec>(-1, static_cast<int>(basis.n_cols) - 1);
  }

  return out;
}

Family* wrap_varying(Family* base, const VaryingSpec& spec) {
  // Whether a drawn coding is exact here is settled by the family that reaches
  // the sampler, not by this one: augmentation can replace a non-quadratic
  // target with a quadratic one, which is how a logit or probit binomial
  // qualifies. `coding_is_exact()` is asked once in `model.cpp`, after the
  // augmentation choice is made.
  return finish(new VaryingCoefficientFamily(base, spec.basis, spec.coding,
                                             spec.coding_levels,
                                             spec.b_labels,
                                             spec.param, spec.column));
}

Family* make_family(const std::string& name, const std::string& link,
                    const arma::vec& y, const arma::vec& w,
                    const List& opts, const arma::mat& vc_basis) {

  Family* base = make_base_family(name, link, y, w, opts);

  if (vc_basis.n_cols == 0) {
    return base;
  }

  if (static_cast<int>(vc_basis.n_rows) != base->N) {
    std::unique_ptr<Family> guard(base);
    stop("the varying-coefficient basis has %d rows and the response has %d.",
         static_cast<int>(vc_basis.n_rows), base->N);
  }

  return wrap_varying(base, varying_spec(opts, vc_basis));
}

namespace {

Family* augmented_base(const std::string& name, const std::string& link,
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

  if (name == "aft" && link == "lognormal" &&
      LognormalAFTAugmentedFamily::applies(w)) {
    return finish(new LognormalAFTAugmentedFamily(
      y, w, as<arma::vec>(opts["event"]), as<double>(opts["sigma_hat"]),
      as<bool>(opts["update_sigma"])));
  }

  if (name == "aft" && link == "loglogistic" &&
      LoglogisticAFTAugmentedFamily::applies(w)) {
    return finish(new LoglogisticAFTAugmentedFamily(
      y, w, as<arma::vec>(opts["event"]), as<double>(opts["sigma_hat"]),
      as<bool>(opts["update_sigma"])));
  }

  if (name == "multinomial" && MultinomAugmentedFamily::applies(w)) {
    return finish(new MultinomAugmentedFamily(y, w, as<int>(opts["num_cat"]),
                                               as<bool>(opts["symmetric"])));
  }

  if ((name == "zip" || name == "zinb") &&
      ZeroInflatedAugmentedFamily::applies(w)) {
    bool nb = name == "zinb";
    return finish(new ZeroInflatedAugmentedFamily(
      y, w, nb,
      nb ? as<double>(opts["theta"]) : 1.0,
      nb ? as<double>(opts["theta_prior_shape"]) : 0.0,
      nb ? as<double>(opts["theta_prior_rate"]) : 0.0,
      nb ? as<bool>(opts["update_theta"]) : false));
  }

  return nullptr;
}

} // namespace

// The rewriting builds a fresh family from the name and the options, so a
// varying-coefficient wrapper put on by `make_family()` would be discarded here.
// Wrapping the rewritten one is what keeps the two paths equivalent -- without
// it, an augmented family reaches the sampler claiming one additive predictor
// while the rest of the fit expects 1 + J of them.
Family* augmented_family(const std::string& name, const std::string& link,
                         const arma::vec& y, const arma::vec& w,
                         const List& opts,
                         const std::vector<std::string>& enabled,
                         const arma::mat& vc_basis) {
  Family* base = augmented_base(name, link, y, w, opts, enabled);

  if (base == nullptr || vc_basis.n_cols == 0) {
    return base;
  }

  return wrap_varying(base, varying_spec(opts, vc_basis));
}

} // namespace bartisan
