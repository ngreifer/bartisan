# The one table `_dev/TASKS.md` carries from this run: the effective sample size
# of the reported cutpoints against the sampler's own and against the adjacent
# gaps, on the same draws. Printed as markdown so the entry quotes it rather
# than a transcription of it.

r <- readRDS("_dev/ordinal-cut-chart.rds")
rows <- r$rows

pick <- function(series) {
  s <- rows[rows$series == series & !is.na(rows$k), ]
  a <- aggregate(cbind(ess_bulk, cor_cut1) ~ k, s, mean)
  a[order(a$k), ]
}

rep_a <- pick("reported")
sam_a <- pick("sampled")
gap_a <- pick("gap")

at <- function(d, k, col = "ess_bulk") {
  v <- d[[col]][d$k == k]
  if (length(v) == 0L) NA_real_ else v
}

cat(sprintf("%d fits: 20 categories, hard and soft rules, n of 500, 2000 and 8000, eight chains of 1000 draws after 200 warmup.\n\n",
            r$done))

cat("| k | reported `aux.cutk` | sampled `c_k` | gap `c_k - c_{k-1}` | cor with `aux.cut1` |\n")
cat("|---|---|---|---|---|\n")

for (k in c(1L, 2L, 3L, 4L, 6L, 10L, 14L, 19L)) {
  g <- at(gap_a, k - 1L)

  cat(sprintf("| %d | %.0f | %s | %s | %.2f |\n", k, at(rep_a, k),
              if (is.na(at(sam_a, k))) "pinned at 0" else sprintf("%.0f", at(sam_a, k)),
              if (is.na(g)) "--" else sprintf("%.0f", g),
              at(rep_a, k, "cor_cut1")))
}

cat(sprintf("\nreported: worst %.0f at k = %d, best %.0f at k = %d\n",
            min(rep_a$ess_bulk), rep_a$k[which.min(rep_a$ess_bulk)],
            max(rep_a$ess_bulk), rep_a$k[which.max(rep_a$ess_bulk)]))
cat(sprintf("sampled:  worst %.0f at k = %d, best %.0f at k = %d\n",
            min(sam_a$ess_bulk), sam_a$k[which.min(sam_a$ess_bulk)],
            max(sam_a$ess_bulk), sam_a$k[which.max(sam_a$ess_bulk)]))
cat(sprintf("gaps:     worst %.0f at k = %d, best %.0f at k = %d\n",
            min(gap_a$ess_bulk), gap_a$k[which.min(gap_a$ess_bulk)],
            max(gap_a$ess_bulk), gap_a$k[which.max(gap_a$ess_bulk)]))

ends <- rows[rows$series %in% c("reported", "sampled") & rows$k == 19L, ]
e <- aggregate(ess_bulk ~ series, ends, function(x) c(min(x), max(x)))
cat(sprintf("\nat k = 19 across the six cells: reported %.0f to %.0f, sampled %.0f to %.0f\n",
            e$ess_bulk[e$series == "reported", 1], e$ess_bulk[e$series == "reported", 2],
            e$ess_bulk[e$series == "sampled", 1], e$ess_bulk[e$series == "sampled", 2]))
