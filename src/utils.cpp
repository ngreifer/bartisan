#include "utils.h"

namespace bartisan {

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

} // namespace bartisan

//' Whether the installed shared library was compiled with optimization
//'
//' Compilers define `__OPTIMIZE__` when they are optimizing, so this is exact
//' rather than a guess. It exists because an unoptimized build of this package
//' is between five and twenty times slower, and nothing else about it looks
//' wrong, which makes it very easy to spend a long time drawing conclusions
//' from the wrong numbers.
//'
//' @returns `TRUE` if the library was optimized.
//' @keywords internal
// [[Rcpp::export(.bartisan_optimized)]]
bool bartisan_optimized() {
#ifdef __OPTIMIZE__
  return true;
#else
  return false;
#endif
}

//' Smoothed empirical distribution function, and the sums a bandwidth selector needs
//'
//' The empirical distribution function convolved with a kernel, evaluated on a
//' grid, together with two sums over the same kernel weights that the bandwidth
//' selectors of Bergmann and Zaehle (2026) are built from.
//'
//' This is in compiled code because the obvious way to write it is quadratic:
//' one pass over every observation for every grid point is `n * m` kernel
//' evaluations, which is tens of millions per predictor on a large fit, and a
//' bandwidth selector repeats the whole thing at every candidate bandwidth.
//'
//' Almost all of that work is avoidable, because almost every weight is a 0 or a
//' 1 rather than something in between. An observation further below a grid point
//' than the kernel reaches has passed it entirely and contributes a weight of 1;
//' one further above has not been reached and contributes 0. Only those within
//' reach need evaluating. With the data and the grid both sorted, the two ends of
//' that window move forward monotonically as the grid advances, so each is found
//' by a pointer that never goes back, and the sweep costs `O(n + m)` plus one
//' kernel evaluation per observation actually inside a window.
//'
//' The Epanechnikov kernel makes this exact rather than approximate, its support
//' being `[-1, 1]`, so the window is `h` either side and nothing outside it is
//' discarded. The Gaussian reaches everywhere and is cut off at five bandwidths,
//' beyond which the weight being rounded to 0 or 1 is under `3e-7`.
//'
//' @param x `numeric`; the observations, **sorted ascending**.
//' @param grid `numeric`; where to evaluate, **sorted ascending**.
//' @param h `numeric`; the bandwidth, strictly positive.
//' @param kernel `integer`; 0 for Epanechnikov, 1 for Gaussian.
//' @returns A matrix with one row per grid point and three columns: the smoothed
//'   distribution function, the sum of squared kernel weights, and the sum of
//'   squared differences between each kernel weight and the step it smooths.
//'   The last two are the terms that distinguish the two selectors.
//' @keywords internal
// [[Rcpp::export(.bartisan_smooth_cdf)]]
Rcpp::NumericMatrix bartisan_smooth_cdf(const std::vector<double>& x,
                                        const std::vector<double>& grid,
                                        double h, int kernel) {
  int n = static_cast<int>(x.size());
  int m = static_cast<int>(grid.size());

  if (n < 1) {
    Rcpp::stop("`x` must have at least one value.");
  }

  if (!(h > 0.0) || !std::isfinite(h)) {
    Rcpp::stop("`h` must be a finite positive bandwidth.");
  }

  double reach = (kernel == 0) ? h : 5.0 * h;

  Rcpp::NumericMatrix out(m, 3);

  // `lo` counts the observations the kernel has entirely passed, each worth a
  // weight of exactly 1; `hi` is one past the last one it has reached at all.
  // `step` counts those at or below the grid point, which is the unsmoothed
  // empirical distribution function. None of the three moves backward, because
  // the grid is sorted.
  int lo = 0;
  int hi = 0;
  int step = 0;

  for (int k = 0; k < m; k++) {
    double v = grid[static_cast<std::size_t>(k)];

    while (lo < n && x[static_cast<std::size_t>(lo)] < v - reach) {
      lo++;
    }

    while (hi < n && x[static_cast<std::size_t>(hi)] <= v + reach) {
      hi++;
    }

    while (step < n && x[static_cast<std::size_t>(step)] <= v) {
      step++;
    }

    // The observations already passed contribute 1 to the weight, 1 to its
    // square, and 0 to the difference from the step, which is also 1 there.
    double s1 = static_cast<double>(lo);
    double s2 = static_cast<double>(lo);
    double s3 = 0.0;

    for (int i = lo; i < hi; i++) {
      double xi = x[static_cast<std::size_t>(i)];
      double u = (v - xi) / h;
      double w;

      if (kernel == 0) {
        // The integral of 0.75 * (1 - u^2) from -1, which is 0 at u = -1 and 1
        // at u = 1, so the two branches above join it continuously.
        w = 0.75 * (u - u * u * u / 3.0) + 0.5;
      }
      else {
        w = R::pnorm(u, 0.0, 1.0, 1, 0);
      }

      double indicator = (xi <= v) ? 1.0 : 0.0;

      s1 += w;
      s2 += w * w;
      s3 += (w - indicator) * (w - indicator);
    }

    // Beyond the window on the low side the weight is 1 and the step is 1, so
    // the difference is 0 and nothing is missed. On the high side both are 0.
    out(k, 0) = s1 / static_cast<double>(n);
    out(k, 1) = s2;
    out(k, 2) = s3;
  }

  return out;
}
