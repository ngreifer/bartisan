#include <RcppArmadillo.h>
#include "polyagamma.h"
#include "utils.h"

//' Polya-Gamma draws, for checking the sampler
//' @param n number of draws.
//' @param b,c parameters of the distribution.
//' @return A numeric vector of draws.
//' @keywords internal
// [[Rcpp::export(.bartisan_rpg)]]
Rcpp::NumericVector bartisan_rpg(int n, double b, double c) {
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = bartisan::rpg(b, c);
  }
  return out;
}

// Exposed for the same reason as the Polya-Gamma sampler above: the truncated
// normal draw is a few lines of log-scale arithmetic that the ordinal
// augmentation leans on for every observation of every sweep, and a quiet bias
// in the far tails would be invisible in a fit.
// [[Rcpp::export(.bartisan_rtruncnorm)]]
Rcpp::NumericVector bartisan_rtruncnorm(int n, double lo, double hi) {
  Rcpp::NumericVector out(n);
  Rcpp::RNGScope scope;

  for (int i = 0; i < n; i++) {
    out[i] = bartisan::truncated_normal_between(lo, hi);
  }

  return out;
}
