# Reads the four cells written by _dev/did-bcf-size.R.
# Usage: Rscript _dev/did-bcf-size-report.R <dir>
dir <- commandArgs(TRUE)[1L]
f <- file.path(dir, paste0("size-", c("A0","B0","A1","B1"), ".rds"))
f <- f[file.exists(f)]
z <- lapply(f, readRDS)
cat("replicates:", paste(sprintf("%s %d/%d%s", sapply(z, `[[`, "cell"),
    sapply(z, `[[`, "done"), sapply(z, `[[`, "total"),
    ifelse(sapply(z, `[[`, "complete"), "", "*")), collapse = "  "), "\n\n")
d <- do.call(rbind, lapply(z, `[[`, "res"))

nullcell <- substr(d$cell, 2L, 2L) == "0"
cat("=== NULL CELLS: rejection rate of the ATT, nominal 0.05 ===\n")
for (cl in intersect(c("A0","B0"), unique(d$cell))) {
  x <- d[d$cell == cl, ]; n <- nrow(x)
  se <- function(p) sqrt(p * (1 - p) / n)
  cat(sprintf("\n%s  (n = %d)\n", cl, n))
  cat(sprintf("  bartisan, 95%% CI excludes 0      %.3f  (+-%.3f)\n",
      mean(x$bart_rej95), 1.96 * se(mean(x$bart_rej95))))
  cat(sprintf("  bartisan, paper's p_Bayes < .05   %.3f  (+-%.3f)   [a 90%% interval]\n",
      mean(x$bart_rej_pb), 1.96 * se(mean(x$bart_rej_pb))))
  cat(sprintf("  TWFE                              %.3f\n", mean(x$twfe_rej, na.rm = TRUE)))
  cat(sprintf("  Callaway-Sant'Anna (did)          %.3f\n", mean(x$cs_rej, na.rm = TRUE)))
  cat(sprintf("  Gardner two-stage (did2s)         %.3f\n", mean(x$g2_rej, na.rm = TRUE)))
  cat(sprintf("  bartisan ATT bias %+.4f   sd %.4f\n", mean(x$bart_att), sd(x$bart_att)))
  cat(sprintf("  share of per-observation CATT 95%% intervals excluding 0: %.4f\n",
      mean(x$catt_excl)))
}

cat("\n\n=== EFFECT CELLS: recovery and power ===\n")
for (cl in intersect(c("A1","B1"), unique(d$cell))) {
  x <- d[d$cell == cl, ]; n <- nrow(x)
  cat(sprintf("\n%s  (n = %d, true ATT %.3f)\n", cl, n, mean(x$att_true)))
  cat(sprintf("  bartisan  bias %+.3f  rmse %.3f  coverage %.2f  power(95%%) %.2f\n",
      mean(x$bart_att - x$att_true), sqrt(mean((x$bart_att - x$att_true)^2)),
      mean(x$bart_lo <= x$att_true & x$att_true <= x$bart_hi), mean(x$bart_rej95)))
  cat(sprintf("  TWFE      bias %+.3f  rmse %.3f                 power     %.2f\n",
      mean(x$twfe_att - x$att_true, na.rm = TRUE),
      sqrt(mean((x$twfe_att - x$att_true)^2, na.rm = TRUE)),
      mean(x$twfe_rej, na.rm = TRUE)))
  cat(sprintf("  CS (did)  bias %+.3f  rmse %.3f                 power     %.2f\n",
      mean(x$cs_att - x$att_true, na.rm = TRUE),
      sqrt(mean((x$cs_att - x$att_true)^2, na.rm = TRUE)),
      mean(x$cs_rej, na.rm = TRUE)))
  cat(sprintf("  did2s     bias %+.3f  rmse %.3f                 power     %.2f\n",
      mean(x$g2_att - x$att_true, na.rm = TRUE),
      sqrt(mean((x$g2_att - x$att_true)^2, na.rm = TRUE)),
      mean(x$g2_rej, na.rm = TRUE)))
  cat(sprintf("  CATT rmse %.3f\n", mean(x$catt_rmse)))
}
