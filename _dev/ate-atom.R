# The atom at zero in the posterior of the ATE.
#
# Observation from the causal vignette: whenever the treatment's `prop_used` is
# below the confidence level, one credible bound comes out exactly zero. The RHC
# fit has `prop_used` 0.949 for `rhc` and reports an ATE of 0.0594 [0, 0.111].
#
# The arithmetic is exact, not approximate. Under the DART sparsity prior a draw
# in which the treatment is in no tree makes the contrast exactly zero, so the
# posterior is a mixture: mass 1 - prop_used at the point zero, and the rest
# spread over nonzero values. If the nonzero part is entirely positive, the
# 2.5% quantile is zero as soon as 1 - prop_used > 0.025, that is as soon as
# prop_used < 0.975. At 0.949 the lower bound must be zero, whatever the data
# say about the size of the effect.
#
# So the question is what to do about it. Three candidates, measured below:
#   sparsity = TRUE     the default, and where the atom comes from
#   sparsity = FALSE    classic BART, every predictor weighted alike
#   split_prior         fixed weights, with extra weight on the treatment

library(bartisan)
library(marginaleffects)

data(rhc)

model <- death ~ rhc + age + sex + race + edu + aps + meanbp + resp +
  hema + pafi + paco2 + crea + surv2m + card

fit_one <- function(label, ...) {
  set.seed(2026)
  fit <- bartisan(model, data = rhc, family = binomial(), chains = 4,
                  control = bartisan_control(...))

  vi <- variable_importance(fit)
  used <- vi[["prop_used"]][vi[["variable"]] == "rhc"]

  ate <- avg_comparisons(fit, variables = "rhc")

  # The atom itself, read off the draws rather than inferred: the share of
  # posterior draws in which the average contrast is exactly zero.
  draws <- avg_comparisons(fit, variables = "rhc")
  post <- attr(draws, "posterior_draws")
  atom <- if (is.null(post)) NA_real_ else mean(post == 0)

  cat(sprintf("%-26s prop_used %.3f  atom %.3f  ATE %.4f [%.4f, %.4f]%s\n",
              label, used, atom, ate$estimate, ate$conf.low, ate$conf.high,
              if (isTRUE(all.equal(ate$conf.low, 0))) "  <- bound exactly 0" else ""))
  invisible(fit)
}

fit_one("sparsity = TRUE (default)")
fit_one("sparsity = FALSE", sparsity = FALSE)
fit_one("split_prior rhc = 1", split_prior = c(rhc = 1))
fit_one("split_prior rhc = 5", split_prior = c(rhc = 5))
fit_one("split_prior rhc = 20", split_prior = c(rhc = 20))
