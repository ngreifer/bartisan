#ifndef GENBART_UTILS_H
#define GENBART_UTILS_H

#include <RcppArmadillo.h>
#include <limits>
#include <cmath>

namespace genbart {

static const double LN_2PI = 1.8378770664093454835606594728112;

// A standard normal restricted to (lo, hi), by inverting the cumulative
// distribution on the log scale. Either bound may be infinite.
//
// The interval is mirrored into the left tail first, so that the two
// probabilities being interpolated between are the small ones. That matters:
// with many cutpoints an ordinal category can sit far enough out that both
// endpoints round to one on the natural scale, and interpolating there collapses
// every draw onto an endpoint.
inline double truncated_normal_between(double lo, double hi) {
  if (!(hi > lo)) {
    return lo;
  }

  if (lo + hi > 0.0) {
    return -truncated_normal_between(-hi, -lo);
  }

  // Lower-tail log probabilities, with lp_lo <= lp_hi.
  double lp_hi = R::pnorm5(hi, 0.0, 1.0, 1, 1);
  double lp_lo = std::isinf(lo) ? R_NegInf : R::pnorm5(lo, 0.0, 1.0, 1, 1);

  double u = unif_rand();

  if (u <= 0.0) {
    u = std::numeric_limits<double>::min();
  }

  // log(p_lo + u * (p_hi - p_lo)), factored through the larger of the two so
  // that the exponential is of a non-positive number.
  double d = lp_lo - lp_hi;
  double inner = u + (1.0 - u) * (std::isinf(d) ? 0.0 : std::exp(d));

  if (!(inner > 0.0)) {
    inner = std::numeric_limits<double>::min();
  }

  double z = R::qnorm5(lp_hi + std::log(inner), 0.0, 1.0, 1, 1);

  // The inversion can land a hair outside; the target requires it inside.
  if (z < lo) {
    return lo;
  }
  if (z > hi) {
    return hi;
  }

  return z;
}

// Numerically stable logistic transform and its logs. expit(x) overflows for
// large negative x if written as 1 / (1 + exp(-x)), so branch on the sign.
inline double expit(double x) {
  if (x >= 0.0) {
    return 1.0 / (1.0 + std::exp(-x));
  }
  double e = std::exp(x);
  return e / (1.0 + e);
}

inline double logit(double p) {
  return std::log(p) - std::log1p(-p);
}

// log(expit(x)) == -log1p(exp(-x))
inline double log_expit(double x) {
  if (x >= 0.0) {
    return -std::log1p(std::exp(-x));
  }
  return x - std::log1p(std::exp(x));
}

// log(1 - expit(x)) == log(expit(-x))
inline double log1m_expit(double x) {
  return log_expit(-x);
}

inline double log_sum_exp(double a, double b) {
  if (a == R_NegInf) {
    return b;
  }
  if (b == R_NegInf) {
    return a;
  }
  double m = std::max(a, b);
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

double log_sum_exp(const arma::vec& x);

// Sample an index from a probability vector, and from a column of a sparse
// matrix of probabilities.
int sample_class(const arma::vec& probs);
int sample_class(int n);
int sample_class_col(const arma::sp_mat& probs, int col);

// log of a Gamma(shape, 1) draw, computed directly for small shapes where the
// gamma draw itself underflows to zero. Method of Liu, Martin and Syring.
double rlgam(double shape);

// Metropolis update of a normal precision under a half-Cauchy prior on the
// standard deviation, using an independence proposal from the conditionally
// conjugate gamma. Used for the residual scale and the leaf scale.
double cauchy_jacobian(double tau, double sigma_hat);
double half_cauchy_update_precision_mh(double sse, double n, double prec_old,
                                       double sigma_scale);
double half_cauchy_update_precision_mh(const arma::vec& r, double prec_old,
                                       double sigma_scale);

// log(1 - exp(-d)) for d > 0, switching form at log(2) to keep both branches
// away from cancellation.
inline double log1mexp(double d) {
  if (!(d > 0.0)) {
    return R_NegInf;
  }
  if (d < M_LN2) {
    return std::log(-std::expm1(-d));
  }
  return std::log1p(-std::exp(-d));
}

// Inverse of the trigamma function by Newton's method, after Smyth's
// implementation in limma. Used to turn a target log-scale spread into the
// shape of a gamma prior.
double trigamma_inverse(double x);

} // namespace genbart

#endif
