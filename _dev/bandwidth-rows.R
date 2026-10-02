# Should `diagnose()` report the gate bandwidth, and how should it be graded?
#
# What is tested. `_dev/bandwidth-mixing.R` found that under soft rules a
# tree's bandwidth carries a median bulk effective sample size of 68 over 1000
# draws with a Poisson likelihood and 190 with a logit, worst tree 9 against
# 114, and `_dev/sbc.R`'s soft-rule Poisson arm deviates from uniform ranks in
# a way that holding the bandwidth fixed removes. `diagnose()` reports every
# other drawn parameter and not this one. The claim: the bandwidth is a slow
# row in its own right across families and sample sizes, slow enough to belong
# in the table, and per tree rather than per observation, so that it wants the
# average-and-worst pair the predictor gets rather than a single number.
#
# What is measured. Soft rules throughout, since a hard rule has no bandwidth.
# `gaussian()`, `binomial()`, `poisson()` and `negbin()` on the Friedman
# function with the family's own response, n of 500, 2000 and 8000, two
# replicates, 50 trees, four chains of 200 warmup and 1000 kept draws. Per fit,
# split R-hat and bulk effective sample size, folded into chains the way
# `diagnose()` folds them, for:
#
#   bandwidth, average over trees      the mean bandwidth across the forest
#   bandwidth, worst 5% of trees       the 95th percentile of R-hat and the
#                                      5th of effective sample size over trees
#   eta, average and worst 5%          as `diagnose()` reports them
#   loglik, splits                     the other rows of the table
#
# and which of those rows carries the fit's worst R-hat.
#
# Read the bandwidth's two rows against the predictor's two, and the count of
# fits in which a bandwidth row is the worst.
#
# What each outcome would mean. A worst-5% bandwidth row that routinely leads
# the fit while its average row does not puts it in the table as a pair, graded
# apart from the reported quantities the way `splits.*` is, since a fit's
# reported quantities are integrals over the bandwidth rather than functions of
# it. A bandwidth that tracks the predictor's rows needs one row and no
# separate grading. A bandwidth that is slow only under the Poisson makes the
# row worth reporting and the advice worth naming the family in.
#
# Run with: Rscript _dev/bandwidth-rows.R
# Writes:   _dev/bandwidth-rows.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

CHAINS <- 4L
DRAWS <- 1000L

friedman <- function(X) {
  10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5]
}

simulate <- function(n, seed, family) {
  set.seed(seed)
  X <- matrix(runif(n * 5L), n, 5L, dimnames = list(NULL, paste0("x", 1:5)))
  f <- friedman(X)
  lp <- (f - mean(f)) / sd(f)
  d <- as.data.frame(X)
  d$y <- switch(family,
                gaussian = lp + rnorm(n, 0, 0.5),
                binomial = rbinom(n, 1L, plogis(1.5 * lp)),
                poisson = rpois(n, exp(lp)),
                negbin = rnbinom(n, mu = exp(lp), size = 3))
  d
}

families <- list(gaussian = gaussian(), binomial = binomial(),
                 poisson = poisson(), negbin = negbin())

grid <- expand.grid(family = names(families), n = c(500L, 2000L, 8000L),
                    rep = 1:2, stringsAsFactors = FALSE)

OUT <- "_dev/bandwidth-rows.rds"

# Resume. Every cell seeds itself before it draws anything, so picking up at
# the next index reproduces exactly what an uninterrupted run would have
# produced; a killed run loses at most the fit it was in. Delete the output
# file to start over.
rows <- list()
START_AT <- 1L

if (file.exists(OUT)) {
  prev <- readRDS(OUT)

  if (!isTRUE(prev$complete) && !is.null(prev$rows) && nrow(prev$rows) > 0L) {
    rows <- lapply(seq_len(nrow(prev$rows)),
                   function(i) prev$rows[i, , drop = FALSE])
    START_AT <- nrow(prev$rows) + 1L
    cat(sprintf("resuming at fit %d of %d\n", START_AT, nrow(grid)))
  }
}

todo <- if (START_AT > nrow(grid)) integer(0) else seq(START_AT, nrow(grid))

pr <- prog_init(total = max(length(todo), 1L),
                title = paste0("Bandwidth as a diagnose() row",
                               if (START_AT > 1L)
                                 sprintf(" (resumed at %d)", START_AT) else ""),
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

# `diagnose()` folds a vector of draws into chains and takes split R-hat and
# rank-normalized effective sample sizes; these reproduce that for one column.
fold <- function(x) matrix(x, ncol = CHAINS)
stat_rhat <- function(x) posterior::rhat(fold(x))
stat_ess <- function(x) posterior::ess_bulk(fold(x))

summarize_columns <- function(wide) {
  rh <- apply(wide, 2L, stat_rhat)
  es <- apply(wide, 2L, stat_ess)
  list(rhat_worst = quantile(rh[is.finite(rh)], 0.95, names = FALSE),
       ess_worst = quantile(es[is.finite(es)], 0.05, names = FALSE),
       rhat_mean = stat_rhat(rowMeans(wide)),
       ess_mean = stat_ess(rowMeans(wide)))
}

for (k in todo) {
  g <- grid[k, ]
  d <- simulate(g$n, 700L + g$rep, g$family)

  t0 <- Sys.time()
  set.seed(800L + g$rep)
  fit <- bartisan(y ~ ., data = d, family = families[[g$family]],
                  control = bartisan_control(num_trees = 50L, num_burn = 200L,
                                             num_draws = DRAWS, chains = CHAINS,
                                             gate = "smoothstep",
                                             verbose = FALSE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  bw <- summarize_columns(fit[["bandwidth"]])
  eta <- summarize_columns(fit[["eta"]][[1L]])
  splits <- rowSums(fit[["counts"]][[1L]])

  row <- data.frame(
    family = g$family, n = g$n, rep = g$rep,
    bw_rhat_mean = bw$rhat_mean, bw_ess_mean = bw$ess_mean,
    bw_rhat_worst = bw$rhat_worst, bw_ess_worst = bw$ess_worst,
    eta_rhat_mean = eta$rhat_mean, eta_ess_mean = eta$ess_mean,
    eta_rhat_worst = eta$rhat_worst, eta_ess_worst = eta$ess_worst,
    loglik_rhat = stat_rhat(fit[["loglik"]]),
    loglik_ess = stat_ess(fit[["loglik"]]),
    splits_rhat = stat_rhat(splits), splits_ess = stat_ess(splits),
    secs = secs)

  candidates <- c(`bandwidth average` = row$bw_rhat_mean,
                  `bandwidth worst 5%` = row$bw_rhat_worst,
                  `eta average` = row$eta_rhat_mean,
                  `eta worst 5%` = row$eta_rhat_worst,
                  loglik = row$loglik_rhat, splits = row$splits_rhat)
  row$worst_row <- names(candidates)[which.max(candidates)]

  rows[[k]] <- row
  prog_tick(pr, label = sprintf("%s / n %d / rep %d", g$family, g$n, g$rep),
            secs = secs)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid),
               chains = CHAINS, draws = DRAWS), OUT)
}

if (!length(todo)) {
  saveRDS(list(rows = do.call(rbind, rows), complete = TRUE), OUT)
}

on.exit()
prog_end(pr, "done")
cat("\nwritten to _dev/bandwidth-rows.rds; read with _dev/ not analyzed here\n")
