# *bartisan* (development version)

* Added `dpm_ls()`, a heteroskedastic `dpm()`: a second forest models the log of the error's standard deviation, so the shape of the error distribution is estimated as `dpm()` estimates it and its scale varies with the predictors.

* Added `dpm_aft_ls()` and `lognormal_aft_ls()`, versions of `dpm_aft()` and `lognormal_aft()` with a second forest on the log of the scale of the log survival time, so that survival curves can cross. `dpm_aft_ls()` is the NFT BART model of Sparapani et al. (2023).

* `posterior_predict()`, and with it `pp_check()` and `simulate()`, now draws event times for `dpm_aft()` fits, which previously had no posterior predictive sampler, and for `dpm_aft_ls()` fits.

* Added `add_chains()` to run more chains of a fitted model and add their draws to the fit. The added chains are the ones a larger `chains` would have run, so the result is identical to a fit that requested the total number of chains from the start with the same seed. `diagnose()` now suggests it where more chains would help.

* A fit with a single chain now draws from the same kind of random number stream as each chain of a multi-chain fit, so its chain is the first chain of a fit with more chains from the same seed. Single-chain results will therefore differ from those of version 0.1.0 with the same seed. Multi-chain `bartisan()` fits run directly after `set.seed()` are unaffected; `bcf()` fits are not when the propensity score is estimated, since that model runs a single chain by default.

* Added interfaces to three prediction workflows: `SL.bartisan()` as a *SuperLearner* library member, a `"bartisan"` engine for `parsnip::bart()` that can be used throughout *tidymodels*, including in a `"censored regression"` mode for the survival families that works with *censored* and the survival metrics of *yardstick*, and `mice.impute.bartisan()`, which lets `mice::mice()` impute with `method = "bartisan"` by drawing from the posterior predictive distribution. See `?SL.bartisan`, `` ?`bartisan-parsnip` ``, and `?mice.impute.bartisan`.

* Fixed a bug in which consecutive multi-chain calls to `bartisan()` or `bcf()` without a call to `set.seed()` in between produced identical draws.

# *bartisan* 0.1.0

* First version!
