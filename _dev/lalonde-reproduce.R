# Does the lalonde fit that "mixes" mix again with a different seed?
#
# The two configurations are the ones actually tried: six chains at 20000 burn
# and 10000 draws, which failed, and four chains at 30000 and 30000, which
# passed. Each is run from several seeds, because the grid in
# `_dev/chains-vs-draws.R` found the per-observation statistics to be bimodal
# across seeds rather than smoothly improving with chain length: at four chains
# and 16000 draws, two seeds gave a worst-5% effective sample size of 8 and 6
# and the third gave 286. If that is what is happening here, a configuration
# that passes once is not a configuration that passes.
#
# The model is `vignette("causal")`'s `fit_earn`, including `sparsity = FALSE`.
#
# Usage: Rscript _dev/lalonde-reproduce.R [seeds]
options(parallelly.availableCores.fallback = 4, parallelly.maxWorkers.localhost = Inf)
suppressMessages({
  library(bartisan); library(future)
})
plan(multisession, workers = 4)

ROOT  <- "/Users/NoahGreifer/Dropbox/Research/R/bartisan"
args  <- commandArgs(TRUE)
SEEDS <- if (length(args) > 0) as.integer(args[[1]]) else 3L

data("lalonde", package = "cobalt")

# The third cell is the prediction the chains-versus-draws argument makes, and
# the reason it is here rather than in the grid: it keeps the total number of
# kept draws and the warmup identical to the configuration that passes and moves
# the draws into more chains. Effective sample size should not notice, and R-hat
# should rise, because twelve chains of 10000 have each seen a third of what
# four chains of 30000 have seen.
cells <- data.frame(chains = c(6L, 4L, 12L), burn = c(20000L, 30000L, 30000L),
                    draws = c(10000L, 30000L, 10000L))

args2 <- if (length(args) > 1) as.integer(strsplit(args[[2]], ",")[[1]]) else seq_len(nrow(cells))

out <- list()
for (i in args2) {
  for (s in seq_len(SEEDS)) {
    ch <- cells[["chains"]][i]
    set.seed(2000L + s)
    t0 <- Sys.time()
    fit <- bartisan(re78 ~ treat + age + educ + race + married + nodegree +
                      re74 + re75, data = lalonde, family = dpm(), chains = ch,
                    num_burn = cells[["burn"]][i], num_draws = cells[["draws"]][i],
                    sparsity = FALSE, verbose = FALSE)
    secs <- as.numeric(Sys.time() - t0, units = "secs")

    d <- diagnose(fit)
    tab <- d[["table"]]
    worst <- tab[grepl("worst 5% of", tab[["quantity"]], fixed = TRUE), ][1L, ]
    avg <- tab[grepl("average over", tab[["quantity"]], fixed = TRUE), ][1L, ]
    warned <- d[["checks"]][["check"]][d[["checks"]][["status"]] == "warn"]

    out[[length(out) + 1L]] <- data.frame(
      chains = ch, burn = cells[["burn"]][i], draws = cells[["draws"]][i],
      seed = 2000L + s, secs = secs,
      w5_rhat = worst[["rhat"]], w5_ess = worst[["ess_bulk"]],
      share_bad = worst[["rhat_bad"]],
      avg_rhat = avg[["rhat"]], avg_ess = avg[["ess_bulk"]],
      passed = length(warned) == 0L,
      warned = paste(warned, collapse = ", "))
    cat(sprintf("%d chains, %d burn, %d draws, seed %d: %.0fs  worst5 rhat %.3f ess %.0f | avg rhat %.3f ess %.0f | %s\n",
                ch, cells[["burn"]][i], cells[["draws"]][i], 2000L + s, secs,
                worst[["rhat"]], worst[["ess_bulk"]], avg[["rhat"]],
                avg[["ess_bulk"]],
                if (length(warned) == 0L) "PASSED" else paste("warned:", paste(warned, collapse = ", "))))
    utils::flush.console()
    rm(fit); invisible(gc())
  }
}
res <- do.call(rbind, out)
saveRDS(res, file.path(ROOT, sprintf("_dev/lalonde-reproduce-%s.rds",
                                     paste(args2, collapse = "-"))))
cat("\nDONE\n")
