## Resubmission

This is a resubmission of a new package, responding to the review of version 0.1.0.

> Please only write package names, software names and API (application programming interface) names in single quotes in title and description. [...] -> please omit them around SoftBart

Done. SoftBart is the name of a model here, and it is now written without quotes.

> Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace (including the package directory and getwd()).

The only code that wrote files was `vignettes/saved-results.R`, a helper the vignettes source so that the output of their slow model fits can be shown without refitting them on CRAN. When `NOT_CRAN` was "true", it ran the fits and wrote their output into `vignettes/results/`. That code is no longer part of the package. The helper now only reads the output saved in `vignettes/results/`, or runs the code without saving anything, and writes no files in either case. No function writes files, and the one test that writes a file does so with `tempfile()` and deletes it afterward.

## Spelling

The incoming check may flag three words in the Description as possibly misspelled. Linero is the surname of the author of the two cited papers, SoftBart is the name of the model introduced in Linero and Yang (2018), and Tweedie is the name of a family of distributions.
