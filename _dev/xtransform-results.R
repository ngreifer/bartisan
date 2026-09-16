SCRATCH <- "/private/tmp/claude-501/-Users-NoahGreifer-Dropbox-Research-R-bartisan/42a89c8f-bd1a-4a13-9717-48692000a78c/scratchpad/xt"
r <- readRDS(file.path(SCRATCH, "results.rds"))
res <- r$res
cat("complete:", r$complete, "  rows:", nrow(res),
    "  failed cells:", sum(is.na(res$transform)), "\n\n")
res <- res[!is.na(res$transform), ]

agg <- aggregate(cbind(rmse, cover, width) ~ scenario + n + transform, res, mean)
w <- reshape(agg, idvar = c("scenario", "n"), timevar = "transform",
             direction = "wide")

# Paired within replicate, which is how the two were run, so the comparison is
# not carrying the between-replicate variation.
pair <- reshape(res[c("scenario", "n", "rep", "transform", "rmse")],
                idvar = c("scenario", "n", "rep"), timevar = "transform",
                direction = "wide")
pair$ratio <- pair$rmse.range / pair$rmse.quantile

ord <- c("uniform", "lognormal", "outliers", "pareto", "outliers_fine",
         "sparse_tail", "bimodal_linear", "mixed")
w <- w[order(match(w$scenario, ord), w$n), ]

cat("RMSE of the posterior mean against the true regression function,\n")
cat("as a fraction of the signal's own SD. Lower is better.\n\n")
cat(sprintf("%-16s %5s | %8s %8s %7s | %s\n",
            "scenario", "n", "quantile", "range", "ratio", "paired: range better in"))
for (i in seq_len(nrow(w))) {
  s <- w$scenario[i]; nn <- w$n[i]
  p <- pair[pair$scenario == s & pair$n == nn, ]
  cat(sprintf("%-16s %5d | %8.4f %8.4f %7.2f | %2d of %2d reps\n",
              s, nn, w$rmse.quantile[i], w$rmse.range[i],
              w$rmse.range[i] / w$rmse.quantile[i],
              sum(p$ratio < 1, na.rm = TRUE), nrow(p)))
}

cat("\n\n95% credible interval coverage of the true function (nominal .95),\n")
cat("and mean interval width as a fraction of the signal SD.\n\n")
cat(sprintf("%-16s %5s | %-17s | %s\n", "scenario", "n",
            "coverage q / r", "width q / r"))
for (i in seq_len(nrow(w))) {
  cat(sprintf("%-16s %5d | %8.3f %8.3f | %7.3f %7.3f\n",
              w$scenario[i], w$n[i], w$cover.quantile[i], w$cover.range[i],
              w$width.quantile[i], w$width.range[i]))
}

cat("\n\nWhy each scenario is here:\n")
for (nm in ord) if (!is.null(r$scenarios[[nm]])) {
  cat(sprintf("  %-16s %s\n", nm, r$scenarios[[nm]]))
}
