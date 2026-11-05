# The report `_dev/sbc.R` prints at the end of a run, and `_dev/sbc-merge.R` at
# the end of a merge, so the two print the same thing for the same replicates.
#
# `out` has one row per replicate with its rank among `L` thinned draws, `L` + 1
# possible values in all. `level` says whether the level's ranks are a valid
# target, which they are only when the run cancelled the empirical anchor.

# The mean-rank and dispersion contrasts, as z statistics against the exact
# moments of a rank uniform on 0..L: the mean rank's departure from L / 2, and
# the mean squared distance from L / 2. A positive dispersion z is a U shape, a
# posterior too narrow; this is the statistic the 2026-10-02 table reports, to
# within 0.08.
rank_z <- function(r, L) {
  v <- 0:L
  dv <- (v - L / 2)^2
  c(mean = (mean(r) - L / 2) / sqrt(mean((v - L / 2)^2) / length(r)),
    dispersion = (mean((r - L / 2)^2) - mean(dv)) /
      sqrt(mean((dv - mean(dv))^2) / length(r)))
}

sbc_report <- function(out, L, label, level = FALSE) {
  # A rank takes the L + 1 values 0 to L, which ten bins split evenly only when
  # L + 1 is a multiple of ten. Otherwise some bins hold one more value than the
  # others, so the expected count is per value, not per bin. Assuming 60 a bin
  # made a uniform L = 20 run read as chi-square 29.7, when it was 9.6.
  bins <- 10L
  breaks <- seq(0, L + 1, length.out = bins + 1L)
  counts <- table(cut(out$rank, breaks = breaks, include.lowest = TRUE))
  values_per_bin <- tabulate(cut(0:L, breaks = breaks, include.lowest = TRUE,
                                 labels = FALSE), nbins = bins)
  expected <- nrow(out) * values_per_bin / (L + 1)
  chisq <- sum((as.numeric(counts) - expected)^2 / expected)

  cat(sprintf("\n%s: %d replicates, %d thinned draws each\n\n", label,
              nrow(out), L))
  cat("rank histogram, 10 bins (uniform is what a correct sampler gives):\n")
  cat(sprintf("  %s\n", paste(sprintf("%4d", as.numeric(counts)), collapse = "")))
  cat(sprintf("  %s expected\n", paste(sprintf("%4.0f", expected), collapse = "")))
  cat(sprintf("\nchi-square %.1f on %d df, p = %.3f\n", chisq, bins - 1L,
              stats::pchisq(chisq, bins - 1L, lower.tail = FALSE)))

  # A U shape is a posterior too narrow, a hump in the middle one too wide, and
  # a slope is a bias. Reporting the three separately says which, where a single
  # p-value does not.
  mid <- mean(out$rank > L * 0.25 & out$rank < L * 0.75)
  mid_uniform <- mean((0:L) > L * 0.25 & (0:L) < L * 0.75)
  cat(sprintf("\nshape: %.0f%% of ranks in the middle half (%.0f%% is uniform)\n",
              100 * mid, 100 * mid_uniform))
  cat(sprintf("       mean rank %.1f of %d (%.1f is uniform)\n",
              mean(out$rank), L, L / 2))
  cat(sprintf("\n95%% interval coverage of the contrast: %.3f (SE %.3f)\n",
              mean(out$covered),
              sqrt(mean(out$covered) * (1 - mean(out$covered)) / nrow(out))))

  z <- rank_z(out$rank, L)
  cat(sprintf("\ncontrast: mean-rank z %+.2f, dispersion z %+.2f\n",
              z[["mean"]], z[["dispersion"]]))

  # The level is a valid target only when the anchor is cancelled, since
  # otherwise the fit's prior on it is centered somewhere the generator's is not.
  if (level) {
    z <- rank_z(out$rank_level, L)
    cat(sprintf("level (predictor at A): mean-rank z %+.2f, dispersion z %+.2f\n",
                z[["mean"]], z[["dispersion"]]))
  }

  invisible(out)
}
