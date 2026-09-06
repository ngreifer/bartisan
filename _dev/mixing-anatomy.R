# What `diagnose()` is actually complaining about on a fit that will not settle,
# and whether the coordinate it complains about is one anybody reports.
#
# The question: the lalonde fit fails its diagnostics at 20000 burn and 10000
# draws over six chains. Is the sampler failing, or is the diagnostic keyed to a
# coordinate that has no bearing on what gets read off the fit? The four
# statistics are computed the way `diagnose()` computes its own rows, but for a
# ladder of quantities: one observation's additive predictor, an average over
# observations, and the treatment effect.
#
# Usage: Rscript _dev/mixing-anatomy.R <draws> [family]
options(parallelly.availableCores.fallback = 4, parallelly.maxWorkers.localhost = Inf)
suppressMessages({
  library(bartisan); library(future)
})
plan(multisession, workers = 4)

ROOT   <- "/Users/NoahGreifer/Dropbox/Research/R/bartisan"
args   <- commandArgs(TRUE)
DRAWS  <- if (length(args) > 0) as.integer(args[[1]]) else 4000L
FAM    <- if (length(args) > 1) args[[2]] else "dpm"
CHAINS <- 4L

data("lalonde", package = "cobalt")
family <- switch(FAM, dpm = dpm(), gaussian = gaussian(), stop("family?"))

set.seed(20260905)
t0 <- Sys.time()
fit <- bartisan(re78 ~ ., data = lalonde, family = family, chains = CHAINS,
                num_burn = DRAWS, num_draws = DRAWS, verbose = FALSE)
secs <- as.numeric(Sys.time() - t0, units = "secs")

d <- diagnose(fit)
cat(sprintf("\n=== lalonde, family = %s, %d chains x %d draws (burn %d), %.0fs\n\n",
            FAM, CHAINS, DRAWS, DRAWS, secs))
print(d)

stats <- function(x) bartisan:::diagnosis_stats(bartisan:::as_chains(x, CHAINS))

# The treatment effect, accumulated over blocks of units so that two full
# posterior prediction matrices never exist at once.
ate <- numeric(DRAWS * CHAINS)
blocks <- split(seq_len(nrow(lalonde)), ceiling(seq_len(nrow(lalonde)) / 100L))
for (b in blocks) {
  d0 <- d1 <- lalonde[b, , drop = FALSE]
  d0[["treat"]] <- 0
  d1[["treat"]] <- 1
  ate <- ate + rowSums(predict(fit, newdata = d1, draws = TRUE) -
                       predict(fit, newdata = d0, draws = TRUE))
}
ate <- ate / nrow(lalonde)

# Both scales, because a family whose own parameters wander can leave the
# additive predictor wandering with them while the response is steady. The
# response side is done a block of units at a time so that a full prediction
# matrix never exists: at 30000 draws over four chains that would be 590 MB.
eta <- fit[["eta"]][[1L]]
by_column <- function(m) vapply(seq_len(ncol(m)), function(j) stats(m[, j]), numeric(4L))
pc_eta <- by_column(eta)

# One prediction per block, used twice: for the per-unit statistics and for the
# running total that gives the average fitted value. Predicting a second time to
# get the average would double the most expensive step here.
by_block <- lapply(blocks, function(b) {
  m <- predict(fit, newdata = lalonde[b, , drop = FALSE], draws = TRUE)
  list(pc = by_column(m), total = rowSums(m))
})
pc_resp <- do.call(cbind, lapply(by_block, `[[`, "pc"))
avg_fitted <- Reduce(`+`, lapply(by_block, `[[`, "total")) / nrow(lalonde)
rm(by_block)

# What the chains' disagreement is worth against the posterior's own width, for
# one unit at a time. R-hat is a ratio and says nothing about the scale of the
# answer; this does.
spread <- function(m) {
  vapply(seq_len(ncol(m)), function(j) {
    ch <- bartisan:::as_chains(m[, j], CHAINS)
    stats::sd(colMeans(ch)) / stats::sd(m[, j])
  }, numeric(1L))
}
sp_eta <- spread(eta)

worst5 <- function(pc) c(stats::quantile(pc[1L, ], 0.95, names = FALSE),
                         stats::quantile(pc[2L, ], 0.95, names = FALSE),
                         stats::quantile(pc[3L, ], 0.05, names = FALSE),
                         stats::quantile(pc[4L, ], 0.05, names = FALSE))

rows <- rbind(
  "one unit's eta, median unit"      = apply(pc_eta,  1L, stats::median),
  "one unit's eta, worst 5%"         = worst5(pc_eta),
  "one unit's fitted mean, median"   = apply(pc_resp, 1L, stats::median),
  "one unit's fitted mean, worst 5%" = worst5(pc_resp),
  "mean of eta over units"           = stats(rowMeans(eta)),
  "mean fitted value over units"     = stats(avg_fitted),
  "ATE"                              = stats(ate))
cat("\n--- the same four statistics, by what is being summarized\n")
print(round(as.data.frame(rows), 3L))

cat(sprintf("\nunits whose eta has R-hat > 1.01:          %3.0f%%\n",
            100 * mean(pc_eta[1L, ] > 1.01, na.rm = TRUE)))
cat(sprintf("units whose fitted mean has R-hat > 1.01:  %3.0f%%\n",
            100 * mean(pc_resp[1L, ] > 1.01, na.rm = TRUE)))
cat(sprintf("\nspread of a unit's chain means, as a fraction of that unit's posterior sd:\n"))
cat(sprintf("   median unit %.3f, worst 5%% %.3f\n",
            stats::median(sp_eta), stats::quantile(sp_eta, 0.95, names = FALSE)))

by_chain <- colMeans(bartisan:::as_chains(ate, CHAINS))
cat(sprintf("\nATE  posterior mean %.1f, posterior sd %.1f\n", mean(ate), stats::sd(ate)))
cat(sprintf("     chain means: %s\n", paste(sprintf("%.1f", by_chain), collapse = "  ")))
cat(sprintf("     spread of chain means / posterior sd = %.3f\n",
            stats::sd(by_chain) / stats::sd(ate)))
cat(sprintf("     95%% interval: [%.0f, %.0f]\n",
            stats::quantile(ate, 0.025), stats::quantile(ate, 0.975)))

saveRDS(list(draws = DRAWS, family = FAM, secs = secs, table = d[["table"]],
             rows = rows, ate = ate, by_chain = by_chain,
             pc_eta = pc_eta, pc_resp = pc_resp, sp_eta = sp_eta,
             avg_eta = rowMeans(eta), loglik = fit[["loglik"]]),
        file.path(ROOT, sprintf("_dev/mixing-anatomy-%s-%d.rds", FAM, DRAWS)))
cat("\nDONE\n")
