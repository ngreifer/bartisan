#ifndef GENBART_POLYAGAMMA_H
#define GENBART_POLYAGAMMA_H

namespace genbart {

// A draw from the Polya-Gamma distribution PG(b, c), for b > 0.
//
// The distribution earns its place through one identity (Polson, Scott and
// Windle 2013). For any likelihood of the form
//
//     exp(kappa * psi) / (1 + exp(psi))^b,
//
// introducing omega ~ PG(b, psi) leaves the conditional density of psi
// proportional to exp(kappa * psi - omega * psi^2 / 2): a Gaussian. Three of
// this package's families are of that form -- binomial with a logit link, the
// negative binomial in its success-probability parameterization, and one
// category of a multinomial conditional on the others -- so for all three the
// augmentation turns a target that is not quadratic in the additive predictor
// into one that is. Everything the sampler does with a quadratic target is then
// exact rather than approximate: Fisher scoring lands in one step and the leaf
// refresh becomes a Gibbs step.
//
// Integer b uses Devroye's alternating-series method, which is exact. A
// non-integer b -- which is what the negative binomial needs, since its b is
// y + theta -- uses the infinite-sum-of-gammas representation truncated and
// completed by its exact mean; see the note in polyagamma.cpp for the size of
// the error that leaves.
double rpg(double b, double c);

} // namespace genbart

#endif
