# *bartisan* (development version)

* Added `add_chains()` to run more chains of a fitted model and add their draws to the fit. The added chains are the ones a larger `chains` would have run, so the result is identical to a fit that requested the total number of chains from the start with the same seed. `diagnose()` now suggests it where more chains would help.

* A fit with a single chain now draws from the same kind of random number stream as each chain of a multi-chain fit, so its chain is the first chain of a fit with more chains from the same seed. Single-chain results will therefore differ from those of version 0.1.0 with the same seed. Multi-chain `bartisan()` fits run directly after `set.seed()` are unaffected; `bcf()` fits are not when the propensity score is estimated, since that model runs a single chain by default.

* Fixed a bug in which consecutive multi-chain calls to `bartisan()` or `bcf()` without a call to `set.seed()` in between produced identical draws.

# *bartisan* 0.1.0

* First version!
