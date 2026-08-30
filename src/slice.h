#ifndef GENBART_SLICE_H
#define GENBART_SLICE_H

#include <RcppArmadillo.h>
#include <functional>

namespace bartisan {

// Univariate stepping-out slice sampler (Neal 2003), taking the log density as
// a std::function so that every nuisance-parameter update can be written as a
// lambda at the call site rather than as its own functor class.
//
// Both loops are capped. Linero's original expands the interval in an
// unbounded `while (true)`, which spins forever if the log density returns a
// non-finite value, and that is reachable whenever a nuisance parameter wanders
// into a region where the likelihood underflows.
double slice_sampler(double x0, const std::function<double(double)>& logf,
                     double w, double lower, double upper,
                     int max_steps = 100);

} // namespace bartisan

#endif
