# Is there anything for a shared forest to win?
#
# The design is the one in Section 4 of Linero, Sinha and Lipsitz (2020), built
# from the paper's description rather than from their code: a binary Z and a
# continuous Y driven by the same Friedman function, conditionally independent
# given x, with Y informative and Z weak. What is scored is how well the model
# recovers Pr(Z = 1 | x), by the cross-entropy between the true and estimated
# probabilities on a fresh test set.
#
# Three arms, and the third is the point of running this before writing any
# engine code:
#
#   separate  a probit fit to Z alone with all P predictors, which is what the
#             package does today and what the paper compares against
#   oracle    the same fit given only the five predictors that matter, which is
#             the ceiling a perfect transfer of variable selection would reach
#   pooled    separate forests drawing their splitting proportions from one
#             pooled Dirichlet, added later if the gap above is worth closing
#
# The oracle is not a method, it is the size of the prize. If it is close to
# `separate` there is nothing to share.
#
# Usage: Rscript _dev/shared-forests-sim.R <reps> <arms>
options(parallelly.availableCores.fallback = 4, parallelly.maxWorkers.localhost = Inf)
suppressMessages({
  library(bartisan); library(future)
})
plan(multisession, workers = 4)

ROOT <- "/Users/NoahGreifer/Dropbox/Research/R/bartisan"
args <- commandArgs(TRUE)
REPS <- if (length(args) > 0) as.integer(args[[1]]) else 20L
ARMS <- if (length(args) > 1) strsplit(args[[2]], ",")[[1]] else
  c("separate", "oracle")

friedman <- function(x) {
  10 * sin(pi * x[, 1L] * x[, 2L]) + 20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] + 5 * x[, 5L]
}

# Standardized, so that Phi(signal * h_std) spans the unit interval the way the
# paper's sigma_theta = 4 does, and so that the signal argument means the same
# thing at every P.
h_std <- function(x) {
  h <- friedman(x)
  (h - 14.4) / 4.8
}

simulate_one <- function(n, p, signal, sd_y, seed) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  colnames(x) <- paste0("x", seq_len(p))
  hs <- h_std(x)
  list(x = x,
       z = stats::rbinom(n, 1L, stats::pnorm(signal * hs)),
       y = friedman(x) + stats::rnorm(n, sd = sd_y),
       pi = stats::pnorm(signal * hs))
}

# The cross-entropy the paper scores, which is the Kullback-Leibler divergence
# from the truth to the estimate, averaged over a test set.
cross_entropy <- function(truth, est) {
  est <- pmin(pmax(est, 1e-6), 1 - 1e-6)
  truth <- pmin(pmax(truth, 1e-6), 1 - 1e-6)
  mean(truth * log(truth / est) + (1 - truth) * log((1 - truth) / (1 - est)))
}

fit_arm <- function(arm, train, test, p) {
  keep <- if (identical(arm, "oracle")) 1:5 else seq_len(p)
  d <- data.frame(z = train[["z"]], train[["x"]][, keep, drop = FALSE])
  nd <- data.frame(test[["x"]][, keep, drop = FALSE])

  fit <- bartisan(z ~ ., data = d, family = stats::binomial("probit"),
                  num_burn = 500, num_draws = 1000, verbose = FALSE)

  as.vector(predict(fit, newdata = nd, type = "response"))
}

cells <- expand.grid(p = c(5L, 20L, 50L, 100L, 250L), signal = 1,
                     KEEP.OUT.ATTRS = FALSE)

out <- list()
for (i in seq_len(nrow(cells))) {
  p <- cells[["p"]][i]
  signal <- cells[["signal"]][i]

  for (r in seq_len(REPS)) {
    train <- simulate_one(250L, p, signal, 1, seed = 1000L * i + r)
    test <- simulate_one(1000L, p, signal, 1, seed = 500000L + 1000L * i + r)

    for (arm in ARMS) {
      t0 <- Sys.time()
      est <- fit_arm(arm, train, test, p)
      out[[length(out) + 1L]] <- data.frame(
        p = p, signal = signal, rep = r, arm = arm,
        loss = cross_entropy(test[["pi"]], est),
        secs = as.numeric(Sys.time() - t0, units = "secs"))
    }
  }
  cat(sprintf("p = %3d done (%d reps x %d arms)\n", p, REPS, length(ARMS)))
  utils::flush.console()
}
res <- do.call(rbind, out)
saveRDS(res, file.path(ROOT, "_dev/shared-forests-sim.rds"))

agg <- aggregate(loss ~ p + arm, data = res, FUN = mean)
wide <- reshape(agg, idvar = "p", timevar = "arm", direction = "wide")
names(wide) <- sub("^loss\\.", "", names(wide))
cat(sprintf("\ncross-entropy loss for Pr(Z = 1 | x), mean over %d reps, n = 250\n\n",
            REPS))
print(wide, row.names = FALSE, digits = 3)

if (all(c("separate", "oracle") %in% ARMS)) {
  cat("\nhow much room a perfect transfer of variable selection would buy:\n")
  for (pp in unique(wide[["p"]])) {
    row <- wide[wide[["p"]] == pp, ]
    cat(sprintf("  p = %3d  separate %.4f  oracle %.4f  ratio %.2fx\n",
                pp, row[["separate"]], row[["oracle"]],
                row[["separate"]] / row[["oracle"]]))
  }
}
cat("\nDONE\n")
