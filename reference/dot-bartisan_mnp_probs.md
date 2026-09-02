# Category probabilities of a multinomial probit fit, by simulation

The probability that the argmax of a correlated Gaussian vector falls in
each category is an orthant probability with no closed form, so it is
simulated. Fresh draws are taken on every call, which makes the estimate
unbiased; the Monte Carlo error is then averaged down by the posterior
draws, so a modest number of replicates per draw is enough for a
posterior mean.

## Usage

``` r
.bartisan_mnp_probs(eta_draws, sigma, replicates)
```

## Arguments

- eta_draws:

  a list of one draws-by-observations matrix per latent variable.

- sigma:

  a matrix of draws by the lower triangle of the covariance matrix,
  column-major within a row, as `aux` stores it.

- replicates:

  simulation replicates per draw and observation.

## Value

An array of draws by observations by categories.
