# Differential test against dbarts, an independent implementation of BART.
#
# What is tested. `_dev/AUDITING.md` sec. 7: where two implementations fit the
# same model, the other's posterior is an oracle. dbarts fits Gaussian and
# probit BART with the Chipman, George and McCulloch tree prior, and bartisan
# fits both with the same prior under hard rules. The claim: on a Gaussian and
# a probit Friedman problem the two posteriors for E[y | x] at held-out points
# agree, to within what their differing scale priors and cutpoint schemes
# allow: bartisan puts a half-Cauchy on the Gaussian residual scale where
# dbarts puts a scaled inverse chi-square, and both set the leaf scale from the
# data.
#
# What is measured. The Friedman function on five uniform predictors, n = 500
# to fit and 500 held out, noise sd 1 for the Gaussian case and a probit
# response from the standardized function for the binary one. Both fits use
# 200 trees, k = 2, base 0.95, power 2, 500 warmup and 1000 kept draws, hard
# rules; bartisan with x_transform = "range" so that its cutpoints sit on the
# raw scale as dbarts' do. From each, the posterior draws of E[y | x] at the
# held-out points: dbarts through predict(type = "ev"), bartisan through
# predict(type = "response", draws = TRUE). Per case:
#
#   cor         correlation of the two posterior means over the test points
#   rel_rmse    RMSE of the difference of the means, over the mean posterior sd
#   sd_ratio    bartisan's mean posterior sd over dbarts'
#   rmse_b, rmse_d   each fit's RMSE against the truth
#   cover_b, cover_d each fit's 95% pointwise coverage of the truth
#   sigma_b, sigma_d the Gaussian residual sd each fit reports
#
# What each outcome means. cor above 0.98, rel_rmse below 0.5 and sd_ratio
# within 0.75 to 1.33 is agreement, and the two priors on sigma account for
# what is left. A shifted or shrunken mean, or an sd ratio outside that band,
# is a difference in the sampler and not in the prior, to be chased.
#
# Run with: Rscript _dev/differential-dbarts.R
# Writes:   _dev/differential-dbarts.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)
suppressMessages(library(dbarts))

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

compare <- function(draws_b, draws_d, truth) {
  mean_b <- colMeans(draws_b)
  mean_d <- colMeans(draws_d)
  sd_b <- apply(draws_b, 2L, sd)
  sd_d <- apply(draws_d, 2L, sd)
  cover <- function(draws) {
    lo <- apply(draws, 2L, quantile, 0.025)
    hi <- apply(draws, 2L, quantile, 0.975)
    mean(lo <= truth & truth <= hi)
  }
  data.frame(cor = cor(mean_b, mean_d),
             rel_rmse = sqrt(mean((mean_b - mean_d)^2)) / mean(c(sd_b, sd_d)),
             sd_ratio = mean(sd_b) / mean(sd_d),
             rmse_b = sqrt(mean((mean_b - truth)^2)),
             rmse_d = sqrt(mean((mean_d - truth)^2)),
             cover_b = cover(draws_b), cover_d = cover(draws_d))
}

pr <- prog_init(total = 2L, title = "Differential test against dbarts", unit = "case",
                kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

# Gaussian.
y <- f_train + rnorm(n, 0, 1)
d_train <- data.frame(y = y, X_train)
d_test <- data.frame(X_test)
set.seed(1)
fit_b <- bartisan(y ~ ., data = d_train, family = gaussian(), control = control)
draws_b <- predict(fit_b, newdata = d_test, type = "response", draws = TRUE)
set.seed(1)
fit_d <- bart(X_train, y, ntree = 200L, k = 2, power = 2, base = 0.95,
              nskip = 500L, ndpost = 1000L, keeptrees = TRUE, verbose = FALSE)
draws_d <- predict(fit_d, newdata = X_test, type = "ev")
gaussian_row <- cbind(case = "gaussian", compare(draws_b, draws_d, f_test),
                      sigma_b = mean(fit_b$aux[, "sigma"]), sigma_d = mean(fit_d$sigma))
prog_tick(pr, label = "gaussian")

# Probit, from the standardized function so that the probabilities span the
# unit interval without piling at its ends.
lat <- 1.2 * (f_train - mean(f_train)) / sd(f_train)
lat_test <- 1.2 * (f_test - mean(f_train)) / sd(f_train)
yb <- rbinom(n, 1L, pnorm(lat))
db_train <- data.frame(yb = yb, X_train)
set.seed(2)
fit_bp <- bartisan(yb ~ ., data = db_train, family = binomial("probit"), control = control)
draws_bp <- predict(fit_bp, newdata = d_test, type = "response", draws = TRUE)
set.seed(2)
fit_dp <- bart(X_train, yb, ntree = 200L, k = 2, power = 2, base = 0.95,
               nskip = 500L, ndpost = 1000L, keeptrees = TRUE, verbose = FALSE)
draws_dp <- predict(fit_dp, newdata = X_test, type = "ev")
probit_row <- cbind(case = "probit", compare(draws_bp, draws_dp, pnorm(lat_test)),
                    sigma_b = NA_real_, sigma_d = NA_real_)
prog_tick(pr, label = "probit")

res <- rbind(gaussian_row, probit_row)
saveRDS(list(rows = res, complete = TRUE), "_dev/differential-dbarts.rds")

on.exit()
prog_end(pr, "done")
print(res, digits = 3, row.names = FALSE)
