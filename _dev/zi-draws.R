# How many draws does a zero-inflated estimand need?
#
# What is tested. `_dev/zi-vc-shortfall.R` found that a `zi_negbin()` fit at
# n = 1000 with 50 trees and 300 warmup and 500 kept draws carries a bulk
# effective sample size of 6 to 14 for the mean over units of the count
# predictor's treatment contrast, and 2 to 66 for the zero predictor's, in
# every structure it was fit in. `vignette("families")` recommends the
# zero-inflated families on their accuracy and says nothing about how many
# draws their estimands need. The claim: the effective sample size of these
# contrasts rises with the number of draws rather than plateauing, so that a
# number can be given.
#
# What is measured. `zi_poisson()` and `zi_negbin()`, n = 1000, three uniform
# predictors, `lp = f(x) + 0.8 z` with f centered, a count mean of
# `exp(lp - 0.5)`, size 3 for the negative binomial, and a structural zero with
# probability 0.25 that does not depend on the treatment. The plain structure
# `y ~ z + x1 + x2 + x3`, 50 trees, one chain, warmup fixed at 300. Both gates,
# since the default is soft and the earlier measurement was hard. Kept draws of
# 500, 2000 and 8000, two data-and-fit seeds. Per fit:
#
#   ess_count, ess_zero    bulk effective sample size of the mean over units of
#                          each predictor's contrast between z set to 1 and 0
#   ess_eta, ess_loglik    the predictor averaged over observations, and the
#                          reported log likelihood
#   effect_count           the count contrast's posterior mean, truth 0.8
#   secs                   wall time, which is not comparable across rows here
#                          because these run beside other jobs
#
# What each outcome would mean. An effective sample size that rises roughly in
# proportion to the draws says the remedy is draws and `vignette("families")`
# can say how many for a stated target, say 400 effective draws. One that
# plateaus says something else binds, and the number of draws is not the advice
# to give; the next thing to read would then be the augmentation, since both
# families draw a latent indicator per observation and a chain that cannot move
# those cannot move the estimand either.
#
# Run with: Rscript _dev/zi-draws.R
# Writes:   _dev/zi-draws.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

simulate <- function(seed) {
  n <- 1000L
  set.seed(seed)
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z
  list(d = d, lp = lp)
}

families <- list(zi_poisson = zi_poisson(), zi_negbin = zi_negbin())

grid <- expand.grid(family = names(families), gate = c("hard", "smoothstep"),
                    draws = c(500L, 2000L, 8000L), seed = 1:2,
                    stringsAsFactors = FALSE)

OUT <- "_dev/zi-draws.rds"

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
                title = paste0("Zero-inflated draws against ESS",
                               if (START_AT > 1L)
                                 sprintf(" (resumed at %d)", START_AT) else ""),
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

ess1 <- function(x) posterior::ess_bulk(matrix(x, ncol = 1L))

for (k in todo) {
  g <- grid[k, ]
  sim <- simulate(900L + g$seed)
  d <- sim$d
  lp <- sim$lp
  d$y <- {
    counts <- if (identical(g$family, "zi_poisson")) rpois(nrow(d), exp(lp - 0.5))
              else rnbinom(nrow(d), mu = exp(lp - 0.5), size = 3)
    counts * rbinom(nrow(d), 1L, 0.75)
  }
  d1 <- d
  d1$z <- 1L
  d0 <- d
  d0$z <- 0L

  t0 <- Sys.time()
  set.seed(1000L + g$seed)
  fit <- bartisan(y ~ z + x1 + x2 + x3, data = d, family = families[[g$family]],
                  control = bartisan_control(num_trees = 50L, num_burn = 300L,
                                             num_draws = g$draws, gate = g$gate,
                                             verbose = FALSE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  e1 <- predict(fit, newdata = d1, type = "link", draws = TRUE)
  e0 <- predict(fit, newdata = d0, type = "link", draws = TRUE)
  count <- rowMeans(e1[[1L]] - e0[[1L]])
  zero <- rowMeans(e1[[2L]] - e0[[2L]])

  rows[[k]] <- data.frame(
    family = g$family, gate = g$gate, draws = g$draws, seed = g$seed,
    ess_count = ess1(count), ess_zero = ess1(zero),
    ess_eta = ess1(rowMeans(fit[["eta"]][[1L]])),
    ess_loglik = ess1(fit[["loglik"]]),
    effect_count = mean(count), effect_zero = mean(zero), secs = secs)
  prog_tick(pr, label = sprintf("%s / %s / %d draws / seed %d", g$family,
                                g$gate, g$draws, g$seed), secs = secs)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)), OUT)
}

if (!length(todo)) {
  saveRDS(list(rows = do.call(rbind, rows), complete = TRUE), OUT)
}

on.exit()
prog_end(pr, "done")
cat("\nwritten to _dev/zi-draws.rds; not analyzed here\n")
