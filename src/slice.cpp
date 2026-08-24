#include "slice.h"

namespace genbart {

double slice_sampler(double x0, const std::function<double(double)>& logf,
                     double w, double lower, double upper, int max_steps) {

  double log_fx0 = logf(x0);
  if (!std::isfinite(log_fx0)) {
    return x0;
  }

  double log_y = log_fx0 - exp_rand();

  double u = w * unif_rand();
  double L = x0 - u;
  double R = x0 + (w - u);

  for (int step = 0; step < max_steps; step++) {
    if (L <= lower) {
      break;
    }
    double fl = logf(L);
    if (!std::isfinite(fl) || fl <= log_y) {
      break;
    }
    L -= w;
  }

  for (int step = 0; step < max_steps; step++) {
    if (R >= upper) {
      break;
    }
    double fr = logf(R);
    if (!std::isfinite(fr) || fr <= log_y) {
      break;
    }
    R += w;
  }

  if (L < lower) {
    L = lower;
  }
  if (R > upper) {
    R = upper;
  }

  // Shrink towards x0 on rejection. Bail out to the current value rather than
  // looping forever if the bracket collapses without an acceptance, which can
  // happen when the log density is not quite unimodal because of round-off.
  for (int step = 0; step < max_steps; step++) {
    double x1 = (R - L) * unif_rand() + L;
    double log_fx1 = logf(x1);

    if (std::isfinite(log_fx1) && log_fx1 >= log_y) {
      return x1;
    }

    if (x1 > x0) {
      R = x1;
    }
    else {
      L = x1;
    }

    if (!(R - L > 0.0)) {
      break;
    }
  }

  return x0;
}

} // namespace genbart
