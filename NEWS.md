# *bartisan* (development version)

* Added `dpm_ls()`, a heteroskedastic `dpm()`: a second forest models the log of the error's standard deviation, so the shape of the error distribution is estimated as `dpm()` estimates it and its scale varies with the predictors.

* Added `dpm_aft_ls()` and `lognormal_aft_ls()`, versions of `dpm_aft()` and `lognormal_aft()` with a second forest on the log of the scale of the log survival time, so that survival curves can cross. `dpm_aft_ls()` is the NFT BART model of Sparapani et al. (2023).

* `posterior_predict()`, and with it `pp_check()` and `simulate()`, now draws event times for `dpm_aft()` fits, which previously had no posterior predictive sampler, and for `dpm_aft_ls()` fits.

* The error distribution of a `dpm()`, `dpm_ls()`, `dpm_aft()`, or `dpm_aft_ls()` fit now has its mixture weights drawn at each iteration given the sizes of its components, rather than fixed at the share of observations in each. With well separated components, whose observations rarely change component, the fixed shares left the uncertainty in the weights out of every density, survival curve, and interval computed from the draws, and the intervals were too narrow. Results differ from those of earlier versions with the same seed.

* The baseline that the mixture components of `dpm()`, `dpm_ls()`, `dpm_aft()`, and `dpm_aft_ls()` are drawn from now has its variance scale and its mean-scale ratio drawn under gamma priors centered at the values they were previously fixed at, which the new argument `update_baseline = FALSE` restores. Fixed at values calibrated from a linear fit, they made an error with narrow, well separated modes come out too wide and too close to its center, and intervals for the error density and for survival curves covered well below their nominal level. Results differ from those of earlier versions with the same seed.

* Added `add_chains()` to run more chains of a fitted model and add their draws to the fit. The added chains are the ones a larger `chains` would have run, so the result is identical to a fit that requested the total number of chains from the start with the same seed. `diagnose()` now suggests it where more chains would help.

* A fit with a single chain now draws from the same kind of random number stream as each chain of a multi-chain fit, so its chain is the first chain of a fit with more chains from the same seed. Single-chain results will therefore differ from those of version 0.1.0 with the same seed. Multi-chain `bartisan()` fits run directly after `set.seed()` are unaffected; `bcf()` fits are not when the propensity score is estimated, since that model runs a single chain by default.

* Added interfaces to three prediction workflows: `SL.bartisan()` as a *SuperLearner* library member, a `"bartisan"` engine for `parsnip::bart()` that can be used throughout *tidymodels*, including in a `"censored regression"` mode for the survival families that works with *censored* and the survival metrics of *yardstick*, and `mice.impute.bartisan()`, which lets `mice::mice()` impute with `method = "bartisan"` by drawing from the posterior predictive distribution. See `?SL.bartisan`, `` ?`bartisan-parsnip` ``, and `?mice.impute.bartisan`.

* Fixed a bug in which consecutive multi-chain calls to `bartisan()` or `bcf()` without a call to `set.seed()` in between produced identical draws.

* Fixed a bug in which the intercept-only fit that the additive predictor starts from ignored the offset, so the forests had to carry the difference between the two against a prior centered at zero. In a Poisson rate model with a log exposure averaging 5, the drawn leaf scale nearly tripled and intervals were about 17% wider than they should have been, and how much wider depended on the units the exposure was measured in. The starting value now includes the offset, so adding a constant to the offset leaves the fit unchanged. Results of fits with an offset differ from those of earlier versions.

* Added `anchor` to `bartisan_control()`, which sets the value on the link scale that each additive predictor is shrunk toward, in place of the intercept-only fit. A fit can then be shrunk toward a known base rate rather than the sample's, or toward a probability of .5, as `dbarts::bart()` does by default for a binary outcome. An offset is added to an anchor rather than absorbed by it.

* Added `hurdle_poisson()`, `hurdle_negbin()`, `hurdle_Gamma()`, and `hurdle_lognormal()`, hurdle (two-part) models in which one forest models the probability of a zero and another the positive values, which for the two count families follow a distribution truncated at zero. See `?bartisan-families` and `vignette("families")`.

# *bartisan* 0.1.0

* First version!
