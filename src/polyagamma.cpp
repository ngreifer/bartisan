#include "polyagamma.h"

#include <RcppArmadillo.h>
#include <cmath>

namespace bartisan {

namespace {

// Where the density of J*(1, 0) switches between its two series expansions.
// Devroye's choice; the two forms agree there to many digits.
const double TRUNC = 0.64;

// Terms kept from the sum-of-gammas representation before the tail is replaced.
// Twenty is enough because the tail is replaced by a gamma matched to both its
// mean and its variance rather than by its mean alone; see rpg_series() for the
// size of what that leaves out.
const int SERIES_TERMS = 20;

// Above this the sum of Devroye draws costs more than it is worth, and the
// series representation is used instead.
const double SUM_MAX = 20.0;

// The nth term of the alternating series for the density of J*(1, 0), in
// whichever of its two forms converges at x.
double series_term(int n, double x) {
  double m = n + 0.5;

  if (x > TRUNC) {
    return M_PI * m * std::exp(-0.5 * m * m * M_PI * M_PI * x);
  }

  return std::pow(2.0 / (M_PI * x), 1.5) * M_PI * m * std::exp(-2.0 * m * m / x);
}

// The distribution function of an inverse Gaussian with mean mu and unit shape.
// An infinite mean is the limit that the c = 0 case needs.
double pinvgauss(double t, double mu) {
  double r = std::sqrt(1.0 / t);

  if (!std::isfinite(mu)) {
    return 2.0 * R::pnorm5(-r, 0.0, 1.0, 1, 0);
  }

  return R::pnorm5(r * (t / mu - 1.0), 0.0, 1.0, 1, 0) +
    std::exp(2.0 / mu) * R::pnorm5(-r * (t / mu + 1.0), 0.0, 1.0, 1, 0);
}

// An inverse Gaussian with mean 1/z and unit shape, truncated to (0, t).
double truncated_invgauss(double z, double t) {
  double mu = z > 0.0 ? 1.0 / z : R_PosInf;

  // When the mean sits above the truncation point, rejection from the
  // untruncated law almost always fails, so the reciprocal is sampled through a
  // truncated chi-square instead.
  if (!(mu <= t)) {
    double x = t;
    double alpha = 0.0;

    while (unif_rand() > alpha) {
      double e1 = exp_rand();
      double e2 = exp_rand();

      while (e1 * e1 > 2.0 * e2 / t) {
        e1 = exp_rand();
        e2 = exp_rand();
      }

      double denom = 1.0 + t * e1;
      x = t / (denom * denom);
      alpha = std::exp(-0.5 * z * z * x);
    }

    return x;
  }

  double x = t + 1.0;

  while (x > t) {
    double y = norm_rand();
    y *= y;
    x = mu + 0.5 * mu * mu * y -
      0.5 * mu * std::sqrt(4.0 * mu * y + mu * mu * y * y);

    if (unif_rand() > mu / (mu + x)) {
      x = mu * mu / x;
    }
  }

  return x;
}

// PG(1, c), by Devroye's method. Exact: the proposal is a two-piece mixture and
// the alternating series decides acceptance without ever truncating.
double rpg1(double c) {
  double z = 0.5 * std::fabs(c);
  double k = 0.125 * M_PI * M_PI + 0.5 * z * z;
  double mu = z > 0.0 ? 1.0 / z : R_PosInf;
  double p = 0.5 * M_PI * std::exp(-k * TRUNC) / k;
  double q = 2.0 * std::exp(-z) * pinvgauss(TRUNC, mu);
  double weight = p / (p + q);

  for (int attempt = 0; attempt < 1000; attempt++) {
    double x = unif_rand() < weight
      ? TRUNC + exp_rand() / k
      : truncated_invgauss(z, TRUNC);

    double s = series_term(0, x);
    double y = unif_rand() * s;

    for (int n = 1; n < 200; n++) {
      if (n % 2 == 1) {
        s -= series_term(n, x);

        if (y <= s) {
          return 0.25 * x;
        }
      }
      else {
        s += series_term(n, x);

        if (y > s) {
          break;
        }
      }
    }
  }

  // Unreachable in practice; the acceptance probability of the outer loop is
  // bounded well away from zero.
  return 0.25 * TRUNC;
}

// sum_{k >= 1} 1 / ((k - 1/2)^2 + a^2), which is (pi / (2a)) tanh(pi a).
double reciprocal_sum(double a) {
  if (a < 1e-8) {
    return 0.5 * M_PI * M_PI;
  }

  return 0.5 * M_PI * std::tanh(M_PI * a) / a;
}

// sum_{k >= 1} 1 / ((k - 1/2)^2 + a^2)^2, which is minus the derivative of the
// above with respect to a^2. The limit at a = 0 is pi^4 / 6.
double reciprocal_square_sum(double a) {
  if (a < 1e-4) {
    return M_PI * M_PI * M_PI * M_PI / 6.0;
  }

  double u = a * a;
  double t = std::tanh(M_PI * a);
  double sech2 = 1.0 - t * t;
  return 0.25 * M_PI * (t / (u * a) - M_PI * sech2 / u);
}

// PG(b, c) from the representation
//
//   PG(b, c) = (1 / (2 pi^2)) sum_{k >= 1} g_k / ((k - 1/2)^2 + c^2/(4 pi^2)),
//
// with g_k independent Gamma(b, 1). Twenty terms are kept and everything beyond
// them is replaced by a single gamma matched to the discarded tail's mean and
// variance, both of which are available in closed form from the two sums above.
//
// What that leaves out is only the *shape* of the tail, and the tail is a sum of
// many tiny independent pieces, so it is very nearly deterministic anyway. At
// twenty terms the tail carries about 1% of the mean and 0.3% of the standard
// deviation of the whole, and matching two moments removes both -- what remains
// is a third-moment discrepancy in a component that small. It is not exactly
// zero, which is why an integer b takes Devroye's exact route instead.
double rpg_series(double b, double c) {
  double a = std::fabs(c) / (2.0 * M_PI);
  double aa = a * a;
  double total = 0.0;
  double kept = 0.0;
  double kept_square = 0.0;

  for (int k = 1; k <= SERIES_TERMS; k++) {
    double d = (k - 0.5) * (k - 0.5) + aa;
    kept += 1.0 / d;
    kept_square += 1.0 / (d * d);
    total += Rf_rgamma(b, 1.0) / d;
  }

  double tail_mean = std::max(reciprocal_sum(a) - kept, 0.0);
  double tail_var = std::max(reciprocal_square_sum(a) - kept_square, 0.0);

  if (tail_mean > 0.0 && tail_var > 0.0) {
    // Gamma(shape, scale) matched to mean b * tail_mean and variance
    // b * tail_var, in the same units as the kept terms.
    double shape = b * tail_mean * tail_mean / tail_var;
    double scale = tail_var / tail_mean;
    total += Rf_rgamma(shape, scale);
  }
  else {
    total += b * tail_mean;
  }

  return total / (2.0 * M_PI * M_PI);
}

} // namespace

double rpg(double b, double c) {
  if (!(b > 0.0) || !std::isfinite(b)) {
    return 0.0;
  }

  double rounded = std::floor(b + 0.5);
  bool whole = std::fabs(b - rounded) < 1e-9;

  if (whole && rounded <= SUM_MAX) {
    double out = 0.0;
    int reps = static_cast<int>(rounded);

    for (int k = 0; k < reps; k++) {
      out += rpg1(c);
    }

    return out;
  }

  return rpg_series(b, c);
}

} // namespace bartisan
