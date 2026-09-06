# Are chains and draws interchangeable?
#
# Two questions that get different answers, which is the point. Effective sample
# size is a statement about how much information the draws carry, and draws from
# a second chain count the same as draws from the end of the first, so at a fixed
# total the arrangement should not matter. R-hat is a statement about whether the
# chains agree, and it is not symmetric at all: more draws per chain drives it
# toward 1, because each chain sees more of the posterior, while more chains give
# it more ways to find a disagreement. Held total fixed, the two effects run in
# opposite directions.
#
# Part 1 holds the total number of kept draws fixed and moves them between
# chains. Part 2 holds the chain count fixed and lengthens the chains.
#
# Usage: Rscript _dev/chains-vs-draws.R [seeds]
options(parallelly.availableCores.fallback = 4, parallelly.maxWorkers.localhost = Inf)
suppressMessages({
  library(bartisan); library(future)
})
plan(multisession, workers = 4)

ROOT  <- "/Users/NoahGreifer/Dropbox/Research/R/bartisan"
args  <- commandArgs(TRUE)
SEEDS <- if (length(args) > 0) as.integer(args[[1]]) else 3L

data("lalonde", package = "cobalt")

one_cell <- function(chains, draws, burn, seed) {
  set.seed(seed)
  t0 <- Sys.time()
  fit <- bartisan(re78 ~ ., data = lalonde, family = dpm(), chains = chains,
                  num_burn = burn, num_draws = draws, verbose = FALSE)
  secs <- as.numeric(Sys.time() - t0, units = "secs")

  stats <- function(x) bartisan:::diagnosis_stats(bartisan:::as_chains(x, chains))
  eta <- fit[["eta"]][[1L]]
  pc  <- vapply(seq_len(ncol(eta)), function(j) stats(eta[, j]), numeric(4L))
  lg  <- stats(fit[["loglik"]])
  av  <- stats(rowMeans(eta))

  data.frame(chains = chains, draws = draws, burn = burn, seed = seed,
             total = chains * draws, secs = secs,
             # what `diagnose()` keys its warnings to
             eta_rhat_worst5 = stats::quantile(pc[1L, ], 0.95, names = FALSE),
             eta_rhat_median = stats::median(pc[1L, ]),
             eta_share_bad   = mean(pc[1L, ] > 1.01, na.rm = TRUE),
             eta_ess_worst5  = stats::quantile(pc[3L, ], 0.05, names = FALSE),
             eta_ess_median  = stats::median(pc[3L, ]),
             # a quantity somebody reports
             avg_rhat = av[["rhat"]], avg_ess = av[["ess_bulk"]],
             loglik_rhat = lg[["rhat"]], loglik_ess = lg[["ess_bulk"]])
}

cells <- rbind(
  # Part 1: 16000 kept draws, arranged four ways.
  data.frame(part = 1L, chains = c(2L, 4L, 8L, 16L),
             draws = c(8000L, 4000L, 2000L, 1000L), burn = 2000L),
  # Part 2: four chains, lengthened.
  data.frame(part = 2L, chains = 4L,
             draws = c(1000L, 2000L, 4000L, 8000L, 16000L),
             burn = c(1000L, 2000L, 4000L, 8000L, 16000L)))

out <- list()
for (s in seq_len(SEEDS)) {
  for (i in seq_len(nrow(cells))) {
    r <- one_cell(cells[["chains"]][i], cells[["draws"]][i], cells[["burn"]][i],
                  1000L + s)
    r[["part"]] <- cells[["part"]][i]
    out[[length(out) + 1L]] <- r
    cat(sprintf("part %d  %2d chains x %5d draws  seed %d  %5.0fs\n",
                r[["part"]], r[["chains"]], r[["draws"]], r[["seed"]], r[["secs"]]))
    utils::flush.console()
  }
}
res <- do.call(rbind, out)
saveRDS(res, file.path(ROOT, "_dev/chains-vs-draws.rds"))

agg <- function(d) {
  a <- aggregate(cbind(secs, eta_rhat_worst5, eta_rhat_median, eta_share_bad,
                       eta_ess_worst5, eta_ess_median, avg_rhat, avg_ess,
                       loglik_rhat, loglik_ess) ~ chains + draws + total,
                 data = d, FUN = mean)
  a[order(a[["chains"]]), ]
}
show <- function(a, title) {
  cat(sprintf("\n%s\n", title))
  cat(sprintf("%7s %7s %7s %6s | %9s %9s %8s | %9s %8s | %8s\n",
              "chains", "draws", "total", "secs", "eta rhat", "eta rhat",
              "% bad", "eta ESS", "eta ESS", "avg ESS"))
  cat(sprintf("%7s %7s %7s %6s | %9s %9s %8s | %9s %8s | %8s\n",
              "", "", "", "", "worst5", "median", "", "worst5", "median", ""))
  for (i in seq_len(nrow(a)))
    cat(sprintf("%7d %7d %7d %6.0f | %9.3f %9.3f %7.0f%% | %9.0f %8.0f | %8.0f\n",
                a[["chains"]][i], a[["draws"]][i], a[["total"]][i], a[["secs"]][i],
                a[["eta_rhat_worst5"]][i], a[["eta_rhat_median"]][i],
                100 * a[["eta_share_bad"]][i], a[["eta_ess_worst5"]][i],
                a[["eta_ess_median"]][i], a[["avg_ess"]][i]))
}
show(agg(res[res[["part"]] == 1L, ]), "Part 1: 16000 kept draws, moved between chains")
show(agg(res[res[["part"]] == 2L, ]), "Part 2: four chains, lengthened")
cat("\nDONE\n")
