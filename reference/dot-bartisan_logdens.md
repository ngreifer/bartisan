# Conditional log density of the outcome at stored posterior draws

Evaluates the family's log density for each observation at each
posterior draw of the additive predictors, using the nuisance parameters
drawn at the same iteration. This is the likelihood contribution of an
observation, so it is a density for a continuous response, a probability
for a discrete one, and a survival probability for a censored survival
time.

## Usage

``` r
.bartisan_logdens(y, weights, eta_draws, family_name, link, family_opts, aux)
```

## Arguments

- y:

  the outcome, coerced as the family expects.

- weights:

  prior weights.

- eta_draws:

  a list of `H` matrices of draws by observations.

- family_name, link, family_opts:

  the family specification.

- aux:

  a matrix of draws by nuisance parameters, with zero columns when the
  family has none.

## Value

A matrix of draws by observations.
