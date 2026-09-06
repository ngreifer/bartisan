# Where in a fit the mixing trouble actually is.
#
# `diagnose()` summarizes the additive predictor over its worst 5% of
# observations, and on a hard dataset that row fails at any number of draws. The
# question this answers is whether the row is describing the fit or only its
# hardest coordinate: the same statistics are computed for one observation at a
# time and for the average over observations, across several datasets and
# families, so that the gap between them can be seen to be a property of forests
# rather than of one dataset.
#
# Usage: Rscript _dev/where-mixing-lives.R [draws]
options(parallelly.availableCores.fallback = 4, parallelly.maxWorkers.localhost = Inf)
suppressMessages({
  library(bartisan); library(future)
})
plan(multisession, workers = 4)

ROOT   <- "/Users/NoahGreifer/Dropbox/Research/R/bartisan"
args   <- commandArgs(TRUE)
DRAWS  <- if (length(args) > 0) as.integer(args[[1]]) else 2000L
CHAINS <- 4L

data("lalonde", package = "cobalt")
data("rhc", package = "bartisan")

set.seed(11L)
n <- 800L
fried <- as.data.frame(matrix(stats::runif(n * 10L), n, 10L))
names(fried) <- paste0("x", 1:10)
fried[["y"]] <- 10 * sin(pi * fried$x1 * fried$x2) + 20 * (fried$x3 - 0.5)^2 +
  10 * fried$x4 + 5 * fried$x5 + stats::rnorm(n)

# `vignette("causal")` fits the lalonde model with `sparsity = FALSE`, which is
# the fit that will not settle, so it gets a cell of its own rather than being
# assumed to behave like the default.
cases <- list(
  list(name = "friedman, gaussian", f = y ~ ., d = fried, fam = gaussian()),
  list(name = "lalonde, gaussian",  f = re78 ~ ., d = lalonde, fam = gaussian()),
  list(name = "lalonde, dpm",       f = re78 ~ ., d = lalonde, fam = dpm()),
  list(name = "lalonde, dpm, no sparsity", f = re78 ~ ., d = lalonde,
       fam = dpm(), sparsity = FALSE),
  list(name = "rhc, binomial",      f = death ~ . - days, d = rhc, fam = binomial()))

`%||%` <- function(a, b) if (is.null(a)) b else a

stats <- function(x) bartisan:::diagnosis_stats(bartisan:::as_chains(x, CHAINS))

out <- list()
for (cs in cases) {
  set.seed(20260905)
  fit <- bartisan(cs[["f"]], data = cs[["d"]], family = cs[["fam"]],
                  chains = CHAINS, num_burn = DRAWS, num_draws = DRAWS,
                  sparsity = cs[["sparsity"]] %||% TRUE, verbose = FALSE)
  eta <- fit[["eta"]][[1L]]
  pc <- vapply(seq_len(ncol(eta)), function(j) stats(eta[, j]), numeric(4L))
  av <- stats(rowMeans(eta))
  out[[length(out) + 1L]] <- data.frame(
    case = cs[["name"]], n = ncol(eta), kept = nrow(eta),
    avg_rhat = av[["rhat"]], avg_ess = av[["ess_bulk"]],
    med_rhat = stats::median(pc[1L, ]), med_ess = stats::median(pc[3L, ]),
    w5_rhat = stats::quantile(pc[1L, ], 0.95, names = FALSE),
    w5_ess = stats::quantile(pc[3L, ], 0.05, names = FALSE),
    share_bad = mean(pc[1L, ] > 1.01, na.rm = TRUE))
  cat(sprintf("%-20s done\n", cs[["name"]])); utils::flush.console()
}
res <- do.call(rbind, out)
saveRDS(res, file.path(ROOT, "_dev/where-mixing-lives.rds"))

cat(sprintf("\n%d chains x %d draws, so %d kept.\n\n", CHAINS, DRAWS, CHAINS * DRAWS))
cat(sprintf("%-20s %6s | %8s %8s | %8s %8s | %8s %8s | %6s\n",
            "", "n", "avg", "avg", "median", "median", "worst5", "worst5", "% over"))
cat(sprintf("%-20s %6s | %8s %8s | %8s %8s | %8s %8s | %6s\n",
            "", "", "R-hat", "ESS", "R-hat", "ESS", "R-hat", "ESS", "1.01"))
for (i in seq_len(nrow(res)))
  cat(sprintf("%-20s %6d | %8.3f %8.0f | %8.3f %8.0f | %8.3f %8.0f | %5.0f%%\n",
              res[["case"]][i], res[["n"]][i], res[["avg_rhat"]][i], res[["avg_ess"]][i],
              res[["med_rhat"]][i], res[["med_ess"]][i], res[["w5_rhat"]][i],
              res[["w5_ess"]][i], 100 * res[["share_bad"]][i]))
cat("\nDONE\n")
