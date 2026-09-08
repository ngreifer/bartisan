# Reads `_dev/shared-topology-sim.R`'s raw output and reports it.
#
#   Rscript _dev/shared-topology-sim.R      # ~25 min on 8 workers
#   Rscript _dev/shared-topology-timing.R   # ~15 min, wants a quiet machine
#   Rscript _dev/shared-topology-results.R  # prints the tables
#
# Two things are reported per component, because the question is not only
# whether sharing costs accuracy but where any cost comes from. The error is
# decomposed at fixed test points:
#
#   bias^2    mean over points of (mean fitted value - truth)^2
#   variance  mean over points of the variance of the fitted value
#   mse       their sum, which is the mean squared error over replicates
#
# A restriction that is wrong shows up as bias; one that is merely tighter than
# it needs to be shows up as a fall in variance. Both together say what sharing
# is doing, which a ratio of root mean squared errors does not.

sim <- readRDS("_dev/shared-topology-sim.rds")
res <- sim$res
cells <- sim$cells

component_names <- list(
  gaussian_ls = c("mean", "log sd"),
  zi_poisson = c("log count", "logit zero"),
  vc = c("control fn", "coefficient")
)

# Everything for one cell and one setting of `shared`, as an array of test
# points by replicates, per component.
gather <- function(model, overlap, p, shared) {
  keep <- Filter(function(x) {
    x$model == model && x$overlap == overlap && x$p == p && x$shared == shared
  }, res)

  list(eta = lapply(seq_len(ncol(keep[[1]]$eta)), function(j) {
         vapply(keep, function(x) x$eta[, j], numeric(nrow(keep[[1]]$eta)))
       }),
       resp = vapply(keep, function(x) x$resp, numeric(length(keep[[1]]$resp))),
       splits = vapply(keep, function(x) x$splits, numeric(length(keep[[1]]$splits))),
       elapsed = vapply(keep, function(x) x$elapsed, numeric(1L)),
       reps = length(keep))
}

# The decomposition, plus the correlation between the averaged fit and the
# truth. The correlation is not interesting in itself and is reported anyway,
# because it is the number that catches a scoring bug: a fit compared against a
# different draw of the test design has the right marginal distribution and a
# correlation of zero, which no amount of staring at a mean squared error
# reveals. An earlier run of this script did exactly that.
decompose <- function(fitted, truth) {
  m <- rowMeans(fitted)
  v <- apply(fitted, 1L, stats::var)
  c(bias2 = mean((m - truth)^2), variance = mean(v),
    mse = mean((m - truth)^2) + mean(v),
    cor = stats::cor(m, truth))
}

fmt <- function(x, d = 4) formatC(x, format = "f", digits = d)

rows <- list()

for (i in seq_len(nrow(cells))) {
  model <- cells$model[i]
  overlap <- cells$overlap[i]
  p <- cells$p[i]
  key <- paste(model, overlap, p)
  tru <- sim$truth[[key]]

  apart <- gather(model, overlap, p, FALSE)
  together <- gather(model, overlap, p, TRUE)

  for (j in seq_along(component_names[[model]])) {
    a <- decompose(apart$eta[[j]], tru$eta[, j])
    b <- decompose(together$eta[[j]], tru$eta[, j])

    rows[[length(rows) + 1L]] <- data.frame(
      model = model, overlap = overlap, p = p,
      component = component_names[[model]][j],
      sep_bias2 = a[["bias2"]], sep_var = a[["variance"]], sep_mse = a[["mse"]],
      shr_bias2 = b[["bias2"]], shr_var = b[["variance"]], shr_mse = b[["mse"]],
      mse_ratio = a[["mse"]] / b[["mse"]],
      sep_cor = a[["cor"]], shr_cor = b[["cor"]],
      stringsAsFactors = FALSE
    )
  }

  a <- decompose(apart$resp, tru$response)
  b <- decompose(together$resp, tru$response)

  rows[[length(rows) + 1L]] <- data.frame(
    model = model, overlap = overlap, p = p, component = "response mean",
    sep_bias2 = a[["bias2"]], sep_var = a[["variance"]], sep_mse = a[["mse"]],
    shr_bias2 = b[["bias2"]], shr_var = b[["variance"]], shr_mse = b[["mse"]],
    mse_ratio = a[["mse"]] / b[["mse"]],
    sep_cor = a[["cor"]], shr_cor = b[["cor"]], stringsAsFactors = FALSE
  )
}

tab <- do.call(rbind, rows)

cat("Shared topology against separate forests.\n")
cat(sprintf("n = %d train, %d test, %d replicates, %d trees, %d burn, %d draws\n\n",
            sim$settings$n_train, sim$settings$n_test, sim$settings$reps,
            sim$settings$num_trees, sim$settings$num_burn,
            sim$settings$num_draws))
cat("mse_ratio > 1 favors sharing. bias2 and var are the decomposition of mse\n")
cat("at fixed test points over replicates.\n\n")

for (model in unique(tab$model)) {
  cat("==", model, "==\n")
  part <- tab[tab$model == model, ]
  out <- data.frame(
    predictors = part$p,
    truth = part$overlap,
    component = part$component,
    sep_bias2 = fmt(part$sep_bias2),
    sep_var = fmt(part$sep_var),
    shr_bias2 = fmt(part$shr_bias2),
    shr_var = fmt(part$shr_var),
    mse_ratio = fmt(part$mse_ratio, 3),
    cor = fmt(part$sep_cor, 2)
  )
  print(out, row.names = FALSE)
  cat("\n")
}

# Loud, because a low correlation here means the fits were scored against the
# wrong data and every number above is meaningless.
worst_cor <- min(c(tab$sep_cor, tab$shr_cor))

if (worst_cor < 0.3) {
  stop("the weakest fit correlates ", formatC(worst_cor, digits = 3),
       " with its truth, which is a scoring bug rather than a result")
}

cat(sprintf("Weakest fit-to-truth correlation across all cells: %.2f\n\n",
            worst_cor))

# Tree complexity, which says whether a shared topology ends up bigger or
# smaller than the separate ones it replaces.
cat("== average splits per draw, per forest ==\n")

splits <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  a <- gather(cells$model[i], cells$overlap[i], cells$p[i], FALSE)
  b <- gather(cells$model[i], cells$overlap[i], cells$p[i], TRUE)
  data.frame(model = cells$model[i], truth = cells$overlap[i], p = cells$p[i],
             separate = paste(fmt(rowMeans(a$splits), 1), collapse = " / "),
             shared = paste(fmt(rowMeans(b$splits), 1), collapse = " / "),
             stringsAsFactors = FALSE)
}))

print(splits, row.names = FALSE)

cat("\n== timings, from _dev/shared-topology-timing.R ==\n")

if (file.exists("_dev/shared-topology-timing.rds")) {
  tm <- readRDS("_dev/shared-topology-timing.rds")
  g <- tm$grid
  print(data.frame(model = g$model, n = g$n, p = g$p, trees = g$num_trees,
                   bw_every = g$bandwidth_every, gate = g$gate,
                   separate = fmt(g$separate, 2), shared = fmt(g$shared, 2),
                   speedup = fmt(g$speedup, 2)), row.names = FALSE)
  cat("\nvarying coefficients, by number of forests:\n")
  print(data.frame(forests = tm$vc_wide$forests,
                   separate = fmt(tm$vc_wide$separate, 2),
                   shared = fmt(tm$vc_wide$shared, 2),
                   speedup = fmt(tm$vc_wide$speedup, 2)), row.names = FALSE)
} else {
  cat("(not run yet)\n")
}

saveRDS(list(accuracy = tab, splits = splits),
        "_dev/shared-topology-results.rds")
