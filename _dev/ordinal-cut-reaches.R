# Does the slow scalar behind `aux.cut1` reach anything a user reports?
#
# WHAT IS BEING TESTED
#
# `_dev/ordinal-cut-chart.R` establishes that `aux.cut1` is not a cutpoint in
# the reported chart but the level of the additive predictor against the pinned
# first threshold, and that it is the one slowly mixing scalar in an ordinal
# fit. The reported predictor has that level subtracted out of it by
# construction, so neither `eta.eta` row can show it: they are about the shape
# of the fitted function, not its height.
#
# That leaves a question the cutpoint comparison cannot answer. If the level
# cancels out of everything a user reports, `aux.cut1` is a row about an
# internal bookkeeping choice and the remedy is to say so plainly. If it does
# not cancel, the row is the only warning a user gets about a quantity they will
# actually compute, and it is misnamed rather than surplus.
#
# The claim under test: **a category probability at the low end of the scale
# inherits `aux.cut1`'s effective sample size rather than the predictor rows'.**
# It could be false. Category probabilities are functions of `c_k - eta_i`,
# which is identified and chart-free, so it is entirely possible that every one
# of them mixes at the predictor's rate and the level is visible nowhere else.
#
# WHAT THE OUTCOME MEASUREMENT IS
#
# One ordinal fit at the main run's cheapest cell: n = 500, 20 categories, the
# same Friedman latent function and seed, hard rules, eight chains of 1000 draws
# after 200 warmup. For 50 held-out observations, the posterior draws of
# `P(Y = k | x)` at k = 1, at k = 10 and at k = 20, each column passed through
# `bartisan:::diagnosis_stats()`, the statistic `diagnose()` reports.
#
# **The column to read is the median `ess_bulk` across observations, at k = 1
# against k = 10 and k = 20**, and both against the fit's own `aux.cut1` and
# `eta.eta (worst 5% of observations)` rows. The cumulative probability
# `P(Y <= 1 | x)` is recorded beside the category probability because it is the
# quantity the first threshold alone determines, so it is the sharper version of
# the same comparison and the one to read if the two disagree.
#
# WHAT EACH OUTCOME WOULD MEAN
#
#   The lowest category's probability mixes near `aux.cut1`'s rate while the
#   middle one mixes near the predictor's. The slow scalar reaches a quantity a
#   user computes, the row belongs in the table, and `?diagnose` has to say that
#   the first reported cutpoint is what governs the low-category probabilities
#   rather than leaving a reader to infer it from a row named `cut1`.
#
#   Every category's probability mixes at the predictor's rate. The level
#   cancels in everything reported, `aux.cut1` is an artifact of the chart with
#   no consequence for a user, and the documentation should say that rather than
#   raising the row's profile.
#
#   The lowest and the highest both mix slowly and the middle does not. The
#   binding quantity is not the level but the tails of the latent distribution,
#   which is a different finding and would want the cutpoint comparison read
#   again before anything is written down.

library(bartisan)

N <- 500L
N_TEST <- 50L
P <- 10L
K_ORD <- 20L
CHAINS <- 8L

FAMILIES <- c("gaussian", "probit", "ordinal")

set.seed(52000 + 211 * 1 + match("ordinal", FAMILIES))
total <- N + N_TEST
X <- matrix(stats::runif(total * P), total, P,
            dimnames = list(NULL, paste0("x", seq_len(P))))
f <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
  10 * X[, 4] + 5 * X[, 5]
eta <- (f - mean(f)) / stats::sd(f)

latent <- eta + stats::rlogis(total)
cuts <- stats::quantile(latent, seq_len(K_ORD - 1L) / K_ORD)
y <- factor(findInterval(latent, cuts), levels = 0:(K_ORD - 1L),
            ordered = TRUE)

dat <- cbind(data.frame(y = y), as.data.frame(X))
train <- dat[seq_len(N), ]
test <- dat[N + seq_len(N_TEST), ]

set.seed(4041)
fit <- bartisan(y ~ ., data = train, family = ordinal(), chains = CHAINS,
                num_burn = 200L, num_draws = 1000L, num_thin = 1L,
                gate = "hard")

probs <- predict(fit, newdata = test, type = "prob", draws = TRUE)

ess_of <- function(m) {
  vapply(seq_len(ncol(m)), function(j) {
    bartisan:::diagnosis_stats(bartisan:::as_chains(m[, j], CHAINS))[["ess_bulk"]]
  }, numeric(1L))
}

cat("median bulk effective sample size over", N_TEST, "held-out observations\n\n")

for (k in c(1L, 10L, 20L)) {
  e <- ess_of(probs[, , k])
  cat(sprintf("P(Y = %2d | x)   %6.0f   (range %.0f to %.0f)\n",
              k, stats::median(e), min(e), max(e)))
}

cum1 <- probs[, , 1L]
e <- ess_of(cum1)
cat(sprintf("\nP(Y <= 1 | x)   %6.0f   (the first threshold alone)\n",
            stats::median(e)))

cat("\nthe fit's own rows:\n")
tb <- diagnose(fit)[["table"]]
keep <- tb$quantity %in% c("aux.cut1", "aux.cut10", "aux.cut19", "loglik") |
  grepl("^eta", tb$quantity)
print(tb[keep, c("quantity", "rhat", "ess_bulk")], row.names = FALSE)

saveRDS(list(probs_ess = lapply(c(1L, 10L, 20L), function(k) ess_of(probs[, , k])),
             table = tb), "_dev/ordinal-cut-reaches.rds")
