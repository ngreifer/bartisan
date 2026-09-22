# The stacked Callaway-Sant'Anna cell construction, shared by the scripts in
# this group and printed in vignette("causal").
#
# For each cohort g and period t, take the cohort-g units plus the never-treated
# (coded first = Inf), set the outcome to the long difference Y_it - Y_i,g-1,
# and take the covariates at g-1. Label the rows with cell = g:t, cohort, and
# since = t - g, and stack every cell into one data frame.
#
# `treat` means "is a cohort-g unit", NOT "is treated now". That is what makes
# t < g - 1 a placebo rather than a separate model, and it is the property the
# whole design rests on. `pre = TRUE` emits those pre-treatment cells.

stack_cells <- function(data, id, time, y, first, covs, pre = FALSE) {
  cells <- list()

  for (g in sort(unique(data[[first]][is.finite(data[[first]])]))) {
    base <- g - 1
    if (!base %in% data[[time]]) next

    keep <- data[[first]] == g | !is.finite(data[[first]])
    pre_rows <- data[keep & data[[time]] == base, ]

    periods <- sort(unique(data[[time]]))
    periods <- if (pre) periods[periods != base] else periods[periods >= g]

    for (t in periods) {
      post <- data[keep & data[[time]] == t, ]
      cell <- merge(pre_rows[c(id, first, covs)], post[c(id, y)], by = id)

      cell[["dY"]] <- cell[[y]] - pre_rows[[y]][match(cell[[id]], pre_rows[[id]])]
      cell[["treat"]] <- as.integer(is.finite(cell[[first]]))
      cell[["cell"]] <- paste(g, t, sep = ":")
      cell[["cohort"]] <- g
      cell[["since"]] <- t - g

      cells[[length(cells) + 1L]] <- cell
    }
  }

  out <- do.call(rbind, cells)
  out[["cell"]] <- factor(out[["cell"]])
  out[["cohort"]] <- factor(out[["cohort"]])
  out[["since"]] <- factor(out[["since"]], levels = sort(unique(out[["since"]])))
  out
}
