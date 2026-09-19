# Is the low-end degradation in `aux.cut1`..`aux.cut19` the sampler's, or the
# chart the cutpoints are reported in?
#
# WHAT IS BEING TESTED
#
# `_dev/TASKS.md` holds two statements that cannot both describe the same
# quantity. The diagnostics entry says an ordinal model's first cutpoint "is
# pinned at zero for identifiability, so there is nothing to diagnose", and that
# R-hat and the effective sample size "now return `NA`, silently, with a test".
# The 48 ordinal fits of `_dev/mixing-against-n.R` report `aux.cut1` with a
# finite R-hat of 1.31 on 23 effective draws, the worst row in every one of
# them.
#
# Reading the source settles the mechanism without fitting anything, and the
# reading is what this script is built on rather than what it tests:
#
#   `update_ordinal_cuts()` (src/family.cpp) loops `for (k = 1; k < num_cat - 1)`,
#   so `cuts(0)` is never written; `prepare_ordered()` (R/response.R) starts it
#   at `raw - raw[1L]`, so `cuts(0)` is zero for the whole run. The sampler's
#   first cutpoint really is pinned.
#
#   `Family::report_shift()` returns `mean(eta.row(0))` and
#   `aux_values_shifted()` returns `cuts - shift`, so the *reported* cutpoints
#   are the sampled ones minus the sample mean of the additive predictor.
#
# Those two together say `aux.cut1` is `-mean(eta)` exactly -- the negated level
# of the fitted function, carrying no threshold at all -- and that `aux.cutk` is
# the identified gap from the first threshold to the k-th, minus that same
# level. So the pinning note and the measurement are about two different charts,
# and the question left is what the gradient across the 19 rows is made of.
#
# The claim under test, which could be false: **the degradation toward the low
# end is one slow scalar, `-mean(eta)`, appearing in all 19 reported rows and
# diluted in each by that row's own variance -- not 19 thresholds mixing badly.**
#
# What makes it falsifiable cheaply is that the sampler's own cutpoints are
# recoverable from the reported draws exactly, with no refitting: with the first
# sampled cutpoint at zero, `cut_k(sampled) = aux.cutk - aux.cut1`. If the
# gradient is the chart's, those recovered series carry no low-end gradient. If
# it is the sampler's, they carry it too, and the TASKS.md entry on the
# tridiagonal joint update is the fix to revisit.
#
# WHAT THE OUTCOME MEASUREMENT IS
#
# The generating process, the seeds and the fit configuration are those of
# `_dev/mixing-against-n.R`, ordinal cells only: a 20-category response from a
# Friedman latent function on ten predictors, hard and soft rules, n of 500,
# 2000 and 8000, eight chains of 1000 draws after 200 warmup, `num_trees` at the
# default. Replicates 1 to 3, so 18 fits. The data are identical to that run's
# first three replicates; the chains are not, since it left the sampler
# unseeded, so the reported column should reproduce its profile in shape and
# level rather than cell for cell.
#
# Every statistic is `bartisan:::diagnosis_stats(as_chains(x, 8))`, the function
# `diagnose()` itself calls, applied to four series built from one `fit$aux`:
#
#   reported   aux[, k]                 k = 1..19   what `diagnose()` prints
#   sampled    aux[, k] - aux[, 1]      k = 2..19   the sampler's own cutpoints
#   gap        aux[, k+1] - aux[, k]    k = 1..18   chart-free and local
#   shift      -aux[, 1]                            the level of the predictor
#
# **The column to read is `ess_bulk` against the cutpoint index**, once for the
# reported series and once for the sampled ones. R-hat is recorded beside it and
# is the endpoint that will mislead here: over eight chains its null is about
# `1 + 8/ess`, so 23 effective draws put it at 1.35 and the observed 1.31 is not
# distinguishable from a fit with nothing wrong. The effective sample size is
# the decisive quantity; R-hat at these sizes is a restatement of it.
#
# Three things are recorded alongside, for interpretation rather than for the
# verdict. The posterior standard deviation of each series and
# `cor(aux[, k], aux[, 1])`, which test the dilution account quantitatively: if
# the reported rows are one slow scalar plus a fast one, their effective sample
# size should fall as the shift's share of their variance rises, and rise as the
# correlation with `aux.cut1` falls. `max(abs(rowMeans(fit$eta[[1]])))`, which
# tests the reading of `report_shift()` directly -- if the shift is the mean of
# the predictor, the recorded predictor has mean exactly zero in every draw.
# And `diagnose(fit)$table` whole, both to confirm this harness reproduces what
# `diagnose()` reports on the same draws and to read the `eta.eta (average over
# observations)` row, which that same identity predicts is rounding noise.
#
# WHAT EACH OUTCOME WOULD MEAN
#
#   The sampled cutpoints are flat in k and mix well while the reported ones
#   carry the gradient. The degradation is the chart. An ordinal fit has one
#   slow scalar in it -- the level of the fitted function against the first
#   threshold -- and the table shows 19 correlated copies of it while the row
#   that used to carry it, `eta.eta (average over observations)`, shows rounding
#   noise. That is a defect in what `diagnose()` reports rather than in the
#   sampler, the tridiagonal entry stays where it is, and TASKS.md and
#   `?diagnose` both have to name the chart they are speaking about.
#
#   The sampled cutpoints degrade toward low k as well. Then the low thresholds
#   genuinely mix worst and the chart only adds to it. That is a sampler
#   finding, the tridiagonal joint update is the entry to revisit, and it wants
#   a design of its own before any code changes.
#
#   Both are flat and the reported gradient does not reproduce. Then this
#   harness is reading something `diagnose()` does not, and the profile from the
#   48 fits has to be re-derived before anything is concluded from it.
#
#   The recorded predictor's mean is not zero to machine precision. Then the
#   reading of `report_shift()` above is wrong, every algebraic step resting on
#   it is unsafe, and the design goes back to the source before it is run again.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

SMOKE <- nzchar(Sys.getenv("CUTCHART_SMOKE"))

N_TEST <- 500L
P <- 10L
K_ORD <- 20L

CHAINS <- 8L
BURN <- 200L
DRAWS <- 1000L
NS <- c(500L, 2000L, 8000L)
REPS <- 3L
OUT <- "_dev/ordinal-cut-chart.rds"

if (SMOKE) {
  CHAINS <- 4L
  BURN <- 20L
  DRAWS <- 60L
  NS <- c(500L)
  REPS <- 1L
  OUT <- "_dev/ordinal-cut-chart-smoke.rds"
}

GATES <- c(hard = "hard", soft = "smoothstep")

# `_dev/mixing-against-n.R`'s generating process for the ordinal family, seeds
# included, so that replicates 1 to 3 are the same data that run used.
FAMILIES <- c("gaussian", "probit", "ordinal")

sim <- function(n, rep) {
  set.seed(52000 + 211 * rep + match("ordinal", FAMILIES))
  total <- n + N_TEST
  X <- matrix(stats::runif(total * P), total, P,
              dimnames = list(NULL, paste0("x", seq_len(P))))
  f <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
    10 * X[, 4] + 5 * X[, 5]
  eta <- (f - mean(f)) / stats::sd(f)

  latent <- eta + stats::rlogis(total)
  cuts <- stats::quantile(latent, seq_len(K_ORD - 1L) / K_ORD)
  y <- factor(findInterval(latent, cuts), levels = 0:(K_ORD - 1L),
              ordered = TRUE)

  dat <- cbind(data.frame(y = y), as.data.frame(X))

  list(train = dat[seq_len(n), ])
}

# One row per series, all four built from the same `fit$aux` and all four passed
# through the statistic `diagnose()` uses, so the reported rows here are the
# rows it prints and the others are strictly comparable to them.
series_rows <- function(aux, chains, rule, n, rep) {
  stat <- function(x) {
    s <- bartisan:::diagnosis_stats(bartisan:::as_chains(x, chains))

    data.frame(rhat = s[["rhat"]], ess_bulk = s[["ess_bulk"]],
               ess_tail = s[["ess_tail"]], sd = stats::sd(x),
               cor_cut1 = stats::cor(x, aux[, 1L]))
  }

  k_max <- ncol(aux)
  parts <- list()

  for (k in seq_len(k_max)) {
    parts[[length(parts) + 1L]] <-
      cbind(data.frame(series = "reported", k = k), stat(aux[, k]))
  }

  # The sampler's own cutpoints: the first is pinned at zero, so the rest come
  # back by subtracting the first reported one, which is the shift.
  for (k in 2:k_max) {
    parts[[length(parts) + 1L]] <-
      cbind(data.frame(series = "sampled", k = k),
            stat(aux[, k] - aux[, 1L]))
  }

  for (k in seq_len(k_max - 1L)) {
    parts[[length(parts) + 1L]] <-
      cbind(data.frame(series = "gap", k = k),
            stat(aux[, k + 1L] - aux[, k]))
  }

  parts[[length(parts) + 1L]] <-
    cbind(data.frame(series = "shift", k = NA_integer_), stat(-aux[, 1L]))

  out <- do.call(rbind, parts)

  cbind(data.frame(rule = rule, n = n, rep = rep), out)
}

run <- function(rule, n, rep) {
  d <- sim(n, rep)

  set.seed(9000 + 17 * rep + n)

  elapsed <- system.time(
    fit <- bartisan(y ~ ., data = d$train, family = ordinal(),
                    chains = CHAINS, num_burn = BURN, num_draws = DRAWS,
                    num_thin = 1L, gate = GATES[[rule]])
  )[["elapsed"]]

  aux <- fit[["aux"]]

  # The reporting identity, tested rather than assumed: if the shift really is
  # the mean of the predictor, the recorded one has mean zero in every draw.
  eta_mean <- rowMeans(fit[["eta"]][[1L]])

  tb <- diagnose(fit)[["table"]]

  rows <- series_rows(aux, CHAINS, rule, n, rep)

  table <- cbind(data.frame(rule = rule, n = n, rep = rep),
                 tb[, c("quantity", "rhat", "ess_bulk", "ess_tail")])

  summary <- data.frame(
    rule = rule, n = n, rep = rep,
    eta_mean_max = max(abs(eta_mean)),
    eta_mean_sd = stats::sd(eta_mean),
    cut1_sd = stats::sd(aux[, 1L]),
    seconds = elapsed)

  list(summary = summary, rows = rows, table = table)
}

cells <- expand.grid(rep = seq_len(REPS), n = NS, rule = names(GATES),
                     stringsAsFactors = FALSE)
cells <- cells[order(cells$rep, cells$rule, cells$n), ]

pr <- prog_init(total = nrow(cells), title = "Ordinal cutpoint chart",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

summaries <- list()
all_rows <- list()
all_tables <- list()

for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]
  label <- sprintf("rep %d: %s n=%d", cell$rep, cell$rule, cell$n)

  out <- prog_do(pr, i, function(.) {
    run(cell$rule, cell$n, cell$rep)
  }, label = label)

  if (!is.null(out)) {
    summaries[[length(summaries) + 1L]] <- out$summary
    all_rows[[length(all_rows) + 1L]] <- out$rows
    all_tables[[length(all_tables) + 1L]] <- out$table
  }

  # After every fit, so a killed run leaves every finished one behind.
  saveRDS(list(res = do.call(rbind, summaries),
               rows = do.call(rbind, all_rows),
               tables = do.call(rbind, all_tables), complete = FALSE,
               done = length(summaries), total = nrow(cells)), OUT)
}

saveRDS(list(res = do.call(rbind, summaries),
             rows = do.call(rbind, all_rows),
             tables = do.call(rbind, all_tables), complete = TRUE,
             done = length(summaries), total = nrow(cells)), OUT)

on.exit()
prog_end(pr, "done")

cat("done:", length(summaries), "of", nrow(cells), "\n")
