# Differential test against stochtree, a second independent implementation.
#
# What is tested. `_dev/AUDITING.md` sec. 7 asks for another implementation's
# answer as an oracle, and `_dev/differential-dbarts.R` does it for *dbarts*.
# *stochtree* fits the same Gaussian and probit BART models with the same
# Chipman, George and McCulloch tree prior, from a different code base again.
# The claim is the same one that held against *dbarts*: on a Gaussian and a
# probit Friedman problem the two posteriors for E[y | x] at held-out points
# agree to within what their differing scale priors and cutpoint schemes allow.
#
# What is measured. The Friedman function on five uniform predictors, n = 500
# to fit and 500 held out, noise standard deviation 1 for the Gaussian case and
# a probit response from the standardized function for the binary one. Both
# fits use 200 trees, a split prior of 0.95 and 2, 500 warmup and 1000 kept
# draws, hard rules, and bartisan uses `x_transform = "range"` so that its
# cutpoints sit on the raw scale. *stochtree* is given `num_gfr = 0` so that it
# runs its plain MCMC rather than the grow-from-root warm start, which bartisan
# has no counterpart for, and its probit model is requested through the current
# `outcome_model` argument. The deprecated `probit_outcome_model` flag it still
# accepts gives a badly compressed fit, the latent predictor spanning -0.62 to
# 1.80 on this design against -4.43 to 5.81, so a comparison made through it
# reads as a disagreement between the packages and is not one. From each, the posterior draws of E[y | x] at the
# held-out points: *stochtree* through `predict()$y_hat`, which it returns as
# observations by draws, and bartisan through
# `predict(type = "response", draws = TRUE)`. The columns are those of
# `_dev/differential-dbarts.R`: the correlation of the two posterior means over
# the test points, the root mean squared difference of those means over the
# mean posterior standard deviation, the ratio of mean posterior standard
# deviations, each fit's error against the truth, each fit's pointwise 95%
# coverage, and the Gaussian residual standard deviation each reports.
#
# What each outcome would mean. Correlation above 0.98, a difference under half
# a posterior standard deviation and a standard deviation ratio within 0.75 to
# 1.33 is agreement, as it was against *dbarts*, and two independent oracles
# then agree with this package on the models all three fit. A departure that
# *dbarts* did not show points at whichever of the two differs from the pair,
# and the first thing to read is the prior each puts on the residual scale and
# the leaf scale.
#
# Run with: Rscript _dev/differential-stochtree.R
# Writes:   _dev/differential-stochtree.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)
suppressMessages(library(stochtree))

set.seed(2026)
n <- 500L
p <- 5L
friedman <- function(X) {
  10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5]
}
X_train <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
X_test <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
f_train <- friedman(X_train)
f_test <- friedman(X_test)

control <- bartisan_control(num_trees = 200L, num_burn = 500L, num_draws = 1000L,
                            gate = "hard", k = 2, gamma = 0.95, beta = 2,
                            x_transform = "range", verbose = FALSE)

compare <- function(draws_b, draws_s, truth) {
  mean_b <- colMeans(draws_b)
  mean_s <- colMeans(draws_s)
  sd_b <- apply(draws_b, 2L, sd)
  sd_s <- apply(draws_s, 2L, sd)
  cover <- function(draws) {
    lo <- apply(draws, 2L, quantile, 0.025)
    hi <- apply(draws, 2L, quantile, 0.975)
    mean(lo <= truth & truth <= hi)
  }
  data.frame(cor = cor(mean_b, mean_s),
             rel_rmse = sqrt(mean((mean_b - mean_s)^2)) / mean(c(sd_b, sd_s)),
             sd_ratio = mean(sd_b) / mean(sd_s),
             rmse_b = sqrt(mean((mean_b - truth)^2)),
             rmse_s = sqrt(mean((mean_s - truth)^2)),
             cover_b = cover(draws_b), cover_s = cover(draws_s))
}

pr <- prog_init(total = 2L, title = "Differential test against stochtree",
                unit = "case", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

# Gaussian.
y <- f_train + rnorm(n, 0, 1)
d_train <- data.frame(y = y, X_train)
d_test <- data.frame(X_test)
set.seed(1)
fit_b <- bartisan(y ~ ., data = d_train, family = gaussian(), control = control)
draws_b <- predict(fit_b, newdata = d_test, type = "response", draws = TRUE)
set.seed(1)
fit_s <- bart(X_train = X_train, y_train = y, X_test = X_test,
              num_gfr = 0L, num_burnin = 500L, num_mcmc = 1000L,
              mean_forest_params = list(num_trees = 200L, alpha = 0.95, beta = 2))
# Returned as observations by draws.
draws_s <- t(predict(fit_s, X_test)$y_hat)
gaussian_row <- cbind(case = "gaussian", compare(draws_b, draws_s, f_test),
                      sigma_b = mean(fit_b$aux[, "sigma"]),
                      sigma_s = mean(sqrt(fit_s$sigma2_global_samples)))
prog_tick(pr, label = "gaussian")

# Probit, from the standardized function.
lat <- 1.2 * (f_train - mean(f_train)) / sd(f_train)
lat_test <- 1.2 * (f_test - mean(f_train)) / sd(f_train)
yb <- rbinom(n, 1L, pnorm(lat))
db_train <- data.frame(yb = yb, X_train)
set.seed(2)
fit_bp <- bartisan(yb ~ ., data = db_train, family = binomial("probit"),
                   control = control)
draws_bp <- predict(fit_bp, newdata = d_test, type = "response", draws = TRUE)
set.seed(2)
fit_sp <- bart(X_train = X_train, y_train = yb, X_test = X_test,
               num_gfr = 0L, num_burnin = 500L, num_mcmc = 1000L,
               general_params = list(
                 outcome_model = OutcomeModel(outcome = "binary",
                                              link = "probit")),
               mean_forest_params = list(num_trees = 200L, alpha = 0.95, beta = 2))
# On the probit scale, so it is put on the probability scale to match.
draws_sp <- pnorm(t(predict(fit_sp, X_test)$y_hat))
probit_row <- cbind(case = "probit", compare(draws_bp, draws_sp, pnorm(lat_test)),
                    sigma_b = NA_real_, sigma_s = NA_real_)
prog_tick(pr, label = "probit")

res <- rbind(gaussian_row, probit_row)
saveRDS(list(rows = res, complete = TRUE), "_dev/differential-stochtree.rds")

on.exit()
prog_end(pr, "done")
print(res, digits = 3, row.names = FALSE)
