SCRATCH <- "/private/tmp/claude-501/-Users-NoahGreifer-Dropbox-Research-R-bartisan/42a89c8f-bd1a-4a13-9717-48692000a78c/scratchpad/xt"
r <- readRDS(file.path(SCRATCH, "candidates.rds"))
res <- r$res[!is.na(r$res$arm), ]
cat("complete:", r$complete, "  rows:", nrow(res), "\n\n")

ARMS <- c("quantile", "range", "robust", "smoothcdf", "winsor")
ord <- c("uniform", "lognormal", "outliers", "pareto", "outliers_fine",
         "sparse_tail", "bimodal_linear", "mixed")

agg <- aggregate(cbind(rmse, cover) ~ scenario + n + arm, res, mean)

cat("RMSE against the true regression function, as a fraction of the signal SD.\n")
cat("Best in each row is starred.\n\n")
cat(sprintf("%-15s %5s", "scenario", "n"))
for (a in ARMS) cat(sprintf(" %10s", a)); cat("\n")
for (s in ord) for (nn in sort(unique(agg$n))) {
  v <- vapply(ARMS, function(a) {
    z <- agg$rmse[agg$scenario == s & agg$n == nn & agg$arm == a]
    if (length(z) == 1L) z else NA_real_
  }, 0)
  if (all(is.na(v))) next
  best <- which.min(v)
  cat(sprintf("%-15s %5d", s, nn))
  for (i in seq_along(ARMS)) {
    cat(sprintf(" %9.4f%s", v[i], if (i == best) "*" else " "))
  }
  cat("\n")
}

cat("\n\nThe decision-relevant summary: how bad each arm gets.\n")
cat("Ratio to the best arm in the same scenario and n, over all 16 cells.\n\n")
cat(sprintf("%-12s %8s %8s %8s   %s\n", "arm", "median", "worst", "mean cov", "worst case"))
wide <- reshape(agg[c("scenario", "n", "arm", "rmse")],
                idvar = c("scenario", "n"), timevar = "arm", direction = "wide")
best <- apply(wide[paste0("rmse.", ARMS)], 1L, min)
for (a in ARMS) {
  ratio <- wide[[paste0("rmse.", a)]] / best
  k <- which.max(ratio)
  cov <- mean(agg$cover[agg$arm == a])
  cat(sprintf("%-12s %8.2f %8.2f %8.3f   %s n=%d\n", a, median(ratio),
              max(ratio), cov, wide$scenario[k], wide$n[k]))
}

cat("\n\n95% interval coverage (nominal .95); below about .9 is a real miss.\n\n")
cat(sprintf("%-15s %5s", "scenario", "n"))
for (a in ARMS) cat(sprintf(" %10s", a)); cat("\n")
for (s in ord) for (nn in sort(unique(agg$n))) {
  v <- vapply(ARMS, function(a) {
    z <- agg$cover[agg$scenario == s & agg$n == nn & agg$arm == a]
    if (length(z) == 1L) z else NA_real_
  }, 0)
  if (all(is.na(v))) next
  cat(sprintf("%-15s %5d", s, nn))
  for (i in seq_along(ARMS)) cat(sprintf(" %10.3f", v[i]))
  cat("\n")
}
