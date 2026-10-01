# Which rows of diagnose() on a dpm() fit measure something a user reads?
#
# What is tested. On the simple Gaussian example recorded in `_dev/TASKS.md`
# (n = 2000, y = sin(4 x1) + x2 + N(0, 0.3^2), three uniform predictors, four
# chains, defaults), `diagnose()` flagged `aux.center` at R-hat 1.63 with a bulk
# ESS of 6.6, `aux.clusters` at 1.34 and `aux.alpha` at 1.26, while the fitted
# values averaged over observations had R-hat 1.00 and ESS 3597. `center` is
# the raw mixture mean, the coordinate along the forest-level ridge that the
# likelihood does not identify, so it is slow by construction and nothing
# reported depends on it. `clusters` and `alpha` describe the mixture's
# representation; the quantity a user reads is the error density that
# representation integrates to, which `error_density()` and
# `predict(type = "density")` report. The claim: the three mixture rows mix
# slowly while every quantity a user reads from the fit -- the predictor, the
# error standard deviation, the error density on a grid, and the log
# likelihood -- mixes acceptably.
#
# What is measured. One fit at the defaults (200 + 800, four chains) and, when
# DPM_DRAWS is set, a longer one. Rows:
#
#   from diagnose(fit)    aux.center, aux.clusters, aux.alpha, aux.error_sd,
#                         loglik, eta (average) and eta (worst 5%)
#   added here            the error density at five points, -2, -1, 0, 1 and 2
#                         residual standard deviations, computed per draw from
#                         the stored mixture with bartisan:::dpm_predictive(),
#                         with split R-hat and bulk ESS from the posterior package
#
# Read R-hat and ESS of the density rows against the mixture rows.
#
# What each outcome means. If the density rows are near R-hat 1.00 with ESS in
# the hundreds while center, clusters and alpha stay far above, the mixture's
# slow coordinates do not reach any reported quantity and diagnose() should
# grade them apart from the reported rows, the way it grades `splits.*`. If the
# density rows are as slow as the mixture rows, the slow mixing reaches
# `predict(type = "density")` and `error_density()`, the rows stay reported,
# and the advice should name those two as what more draws are for. A longer run
# that brings the mixture rows down says the problem is mixing rather than
# chains settling in different places.
#
# Run with: Rscript _dev/dpm-diagnose-rows.R
#           DPM_DRAWS=3200 Rscript _dev/dpm-diagnose-rows.R
# Writes:   _dev/dpm-diagnose-rows-<draws>.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

DRAWS <- as.integer(Sys.getenv("DPM_DRAWS", "800"))
CHAINS <- 4L
OUT <- sprintf("_dev/dpm-diagnose-rows-%d.rds", DRAWS)

set.seed(1001)
n <- 2000L
d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
d$y <- sin(4 * d$x1) + d$x2 + rnorm(n, sd = 0.3)

pr <- prog_init(total = 2L, title = sprintf("dpm diagnose rows (%d draws)", DRAWS),
                unit = "step", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

fit <- bartisan(y ~ x1 + x2 + x3, data = d, family = dpm(), chains = CHAINS,
                num_burn = 200L, num_draws = DRAWS)
prog_tick(pr, label = "fit")

tab <- as.data.frame(diagnose(fit)[["table"]])
keep <- grepl("^aux\\.|^loglik|^eta", tab[["quantity"]])
from_fit <- tab[keep, c("quantity", "rhat", "ess_bulk", "ess_tail")]

# The error density at a few residual values, per draw, from the stored mixture.
# The residual sd is 0.3 by construction, so the grid is -2 to 2 of it.
at <- c(-0.6, -0.3, 0, 0.3, 0.6)
S <- nrow(fit[["aux"]])
dens <- t(vapply(seq_len(S), function(s) bartisan:::dpm_predictive(fit, s, at),
                 numeric(length(at))))

by_chain <- function(x) matrix(x, ncol = CHAINS)
density_rows <- data.frame(
  quantity = sprintf("density at %+.1f", at),
  rhat = apply(dens, 2L, function(x) posterior::rhat(by_chain(x))),
  ess_bulk = apply(dens, 2L, function(x) posterior::ess_bulk(by_chain(x))),
  ess_tail = apply(dens, 2L, function(x) posterior::ess_tail(by_chain(x))))

out <- rbind(from_fit, density_rows)
rownames(out) <- NULL
prog_tick(pr, label = "rows")

saveRDS(list(rows = out, draws = DRAWS, chains = CHAINS, complete = TRUE), OUT)

on.exit()
prog_end(pr, "done")

print(out, digits = 3)
