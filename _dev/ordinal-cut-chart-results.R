# Reads `_dev/ordinal-cut-chart.rds` and prints the comparison the run was for:
# the effective sample size of the reported cutpoints against that of the
# sampled ones, on the same draws. See the header of `_dev/ordinal-cut-chart.R`
# for what is being tested and what each outcome would mean.

r <- readRDS("_dev/ordinal-cut-chart.rds")

cat(sprintf("%d of %d fits, complete = %s\n\n", r$done, r$total, r$complete))

rows <- r$rows

cat("== the reporting identity ==\n")
print(r$res[, c("rule", "n", "rep", "eta_mean_max", "eta_mean_sd", "cut1_sd")])

cat(sprintf("\n== effective sample size by cutpoint index, averaged over %d fits ==\n", r$done))

wide <- function(series) {
  s <- rows[rows$series == series & !is.na(rows$k), ]
  a <- aggregate(cbind(ess_bulk, rhat, sd, cor_cut1) ~ k, s, mean)
  a[order(a$k), ]
}

rep_a <- wide("reported")
sam_a <- wide("sampled")
gap_a <- wide("gap")

cat("\n k | reported: ess  rhat   sd   cor(.,cut1) | sampled: ess  rhat   sd | gap: ess\n")

for (k in rep_a$k) {
  s <- sam_a[sam_a$k == k, ]
  g <- gap_a[gap_a$k == k, ]

  cat(sprintf("%2d | %14.0f %6.3f %6.3f %8.3f | %11.0f %6.3f %5.3f | %8s\n",
              k, rep_a$ess_bulk[rep_a$k == k], rep_a$rhat[rep_a$k == k],
              rep_a$sd[rep_a$k == k], rep_a$cor_cut1[rep_a$k == k],
              if (nrow(s)) s$ess_bulk else NA,
              if (nrow(s)) s$rhat else NA,
              if (nrow(s)) s$sd else NA,
              if (nrow(g)) sprintf("%.0f", g$ess_bulk) else "-"))
}

sh <- rows[rows$series == "shift", ]
cat(sprintf("\nshift (-aux.cut1): ess %.0f, rhat %.3f, sd %.3f\n",
            mean(sh$ess_bulk), mean(sh$rhat), mean(sh$sd)))

cat("\n== the decisive contrast ==\n")
cat(sprintf("reported, worst index:  ess %.0f at k = %d\n",
            min(rep_a$ess_bulk), rep_a$k[which.min(rep_a$ess_bulk)]))
cat(sprintf("reported, best index:   ess %.0f at k = %d\n",
            max(rep_a$ess_bulk), rep_a$k[which.max(rep_a$ess_bulk)]))
cat(sprintf("sampled,  worst index:  ess %.0f at k = %d\n",
            min(sam_a$ess_bulk), sam_a$k[which.min(sam_a$ess_bulk)]))
cat(sprintf("sampled,  best index:   ess %.0f at k = %d\n",
            max(sam_a$ess_bulk), sam_a$k[which.max(sam_a$ess_bulk)]))
cat(sprintf("gap,      worst index:  ess %.0f at k = %d\n",
            min(gap_a$ess_bulk), gap_a$k[which.min(gap_a$ess_bulk)]))

cat("\n== by sample size and rule, at the two ends ==\n")

ends <- rows[rows$series %in% c("reported", "sampled") & rows$k %in% c(2L, 19L), ]
print(aggregate(ess_bulk ~ series + k + rule + n, ends, mean))

cat("\n== the dilution account ==\n")
cat("If the reported rows are one slow scalar plus a fast one, their effective\n")
cat("sample size should track the shift's share of their variance.\n\n")

v_shift <- mean(sh$sd)^2
rep_a$share <- v_shift / rep_a$sd^2
cat(sprintf("%2s %10s %10s %9s\n", "k", "var share", "cor(.,cut1)", "ess"))

for (k in rep_a$k) {
  i <- rep_a$k == k
  cat(sprintf("%2d %10.3f %10.3f %9.0f\n", k, rep_a$share[i], rep_a$cor_cut1[i],
              rep_a$ess_bulk[i]))
}

cat(sprintf("\nSpearman correlation of ess with the shift's variance share: %.3f\n",
            stats::cor(rep_a$share, rep_a$ess_bulk, method = "spearman")))
cat(sprintf("Spearman correlation of ess with cor(., cut1):               %.3f\n",
            stats::cor(rep_a$cor_cut1, rep_a$ess_bulk, method = "spearman")))

# The sharp version. Writing the reported cutpoint as `c_k - m`, the variance of
# its posterior mean is the sum of the two pieces' own, so
#
#   ess(c_k - m) = var(c_k - m) / [var(c_k)/ess(c_k) + var(m)/ess(m)]
#
# which predicts every reported row from two measured series and nothing else.
# The cross-covariance is dropped, so this is an account rather than an identity
# and the question is whether it lands on the right order.
cat("\n== predicting each reported row from the sampled one and the shift ==\n")
cat(" k   observed   predicted   ratio\n")

ess_m <- mean(sh$ess_bulk)
v_m <- mean(sh$sd)^2

for (k in 2:max(rep_a$k)) {
  i <- rep_a$k == k
  j <- sam_a$k == k
  v_rep <- rep_a$sd[i]^2
  v_c <- sam_a$sd[j]^2
  pred <- v_rep / (v_c / sam_a$ess_bulk[j] + v_m / ess_m)

  cat(sprintf("%2d %10.0f %11.0f %7.2f\n", k, rep_a$ess_bulk[i], pred,
              rep_a$ess_bulk[i] / pred))
}

cat("\n== R-hat against its own null, 1 + chains/ess ==\n")
cat("The package's own calibration. An excess at or below zero means the chains\n")
cat("agree as well as eight chains at that effective sample size ever do.\n\n")

for (k in rep_a$k) {
  i <- rep_a$k == k
  null <- 1 + 8 / rep_a$ess_bulk[i]
  cat(sprintf("k %2d  observed %.3f  null %.3f  excess %+.3f\n",
              k, rep_a$rhat[i], null, rep_a$rhat[i] - null))
}

cat("\n== the averaged-predictor row diagnose() prints ==\n")
avg <- r$tables[grepl("average over observations", r$tables$quantity), ]
cat(sprintf("ess_bulk: median %.0f, range %.0f to %.0f, on a series whose sd is %.1e\n",
            stats::median(avg$ess_bulk), min(avg$ess_bulk), max(avg$ess_bulk),
            mean(r$res$eta_mean_sd)))

cat("\n== harness against diagnose() on the reported cutpoints ==\n")
tb <- r$tables[grepl("^aux.cut", r$tables$quantity), ]
tb$k <- as.integer(sub("aux.cut", "", tb$quantity))
mine <- rows[rows$series == "reported", ]
m <- merge(tb, mine, by = c("rule", "n", "rep", "k"))
cat(sprintf("%d rows compared, max |difference| in ess %.3g, in rhat %.3g\n",
            nrow(m), max(abs(m$ess_bulk.x - m$ess_bulk.y)),
            max(abs(m$rhat.x - m$rhat.y))))
