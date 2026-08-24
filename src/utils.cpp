#include "utils.h"

namespace genbart {

double log_sum_exp(const arma::vec& x) {
  double m = x.max();
  if (!std::isfinite(m)) {
    return m;
  }
  return m + std::log(arma::sum(arma::exp(x - m)));
}

int sample_class(const arma::vec& probs) {
  double u = unif_rand();
  double cumulative = 0.0;
  int k_max = static_cast<int>(probs.n_elem);
  for (int k = 0; k < k_max; k++) {
    cumulative += probs(k);
    if (u < cumulative) {
      return k;
    }
  }
  return k_max - 1;
}

int sample_class(int n) {
  int out = static_cast<int>(std::floor(unif_rand() * n));
  return out >= n ? n - 1 : out;
}

int sample_class_col(const arma::sp_mat& probs, int col) {
  double u = unif_rand();
  double cumulative = 0.0;
  arma::sp_mat::const_col_iterator it = probs.begin_col(col);
  arma::sp_mat::const_col_iterator it_end = probs.end_col(col);
  int last = 0;
  for (; it != it_end; ++it) {
    cumulative += (*it);
    last = static_cast<int>(it.row());
    if (u < cumulative) {
      return last;
    }
  }
  return last;
}

double rlgam(double shape) {
  if (shape >= 0.1) {
    return std::log(Rf_rgamma(shape, 1.0));
  }

  double a = shape;
  double L = 1.0 / a - 1.0;
  double w = std::exp(-1.0) * a / (1.0 - a);
  double ww = 1.0 / (1.0 + w);
  double z = 0.0;
  do {
    double u = unif_rand();
    if (u <= ww) {
      z = -std::log(u / ww);
    }
    else {
      z = std::log(unif_rand()) / L;
    }
    double eta = z >= 0 ? -z : std::log(w) + std::log(L) + L * z;
    double h = -z - std::exp(-z / a);
    if (h - eta > std::log(unif_rand())) {
      break;
    }
  } while (true);

  return -z / a;
}

double cauchy_jacobian(double tau, double sigma_hat) {
  double sigma = std::pow(tau, -0.5);
  double out = Rf_dcauchy(sigma, 0.0, sigma_hat, 1);
  out += -M_LN2 - 1.5 * std::log(tau);
  return out;
}

double half_cauchy_update_precision_mh(const arma::vec& r, double prec_old,
                                       double sigma_scale) {
  return half_cauchy_update_precision_mh(arma::dot(r, r),
                                         static_cast<double>(r.n_elem),
                                         prec_old, sigma_scale);
}

double half_cauchy_update_precision_mh(double sse, double n, double prec_old,
                                       double sigma_scale) {
  if (!(sse > 0.0) || !(n > 0.0)) {
    return prec_old;
  }
  double shape = 0.5 * n + 1.0;
  double scale = 2.0 / sse;
  double tau_prop = Rf_rgamma(shape, scale);
  if (!(tau_prop > 0.0) || !std::isfinite(tau_prop)) {
    return prec_old;
  }

  double log_ratio = cauchy_jacobian(tau_prop, sigma_scale) -
    cauchy_jacobian(prec_old, sigma_scale);

  return std::log(unif_rand()) < log_ratio ? tau_prop : prec_old;
}

double trigamma_inverse(double x) {
  if (x > 1e7) {
    return 1.0 / std::sqrt(x);
  }
  if (x < 1e-6) {
    return 1.0 / x;
  }

  double y = 0.5 + 1.0 / x;
  for (int i = 0; i < 50; i++) {
    double tri = R::trigamma(y);
    double dif = tri * (1.0 - tri / x) / R::tetragamma(y);
    y += dif;
    if (-dif / y < 1e-8) {
      break;
    }
  }
  return y;
}

} // namespace genbart

//' Whether the installed shared library was compiled with optimization
//'
//' Compilers define `__OPTIMIZE__` when they are optimizing, so this is exact
//' rather than a guess. It exists because an unoptimized build of this package
//' is between five and twenty times slower, and nothing else about it looks
//' wrong -- which makes it very easy to spend a long time drawing conclusions
//' from the wrong numbers.
//'
//' @return `TRUE` if the library was optimized.
//' @keywords internal
// [[Rcpp::export(.genbart_optimized)]]
bool genbart_optimized() {
#ifdef __OPTIMIZE__
  return true;
#else
  return false;
#endif
}
