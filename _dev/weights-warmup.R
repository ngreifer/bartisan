# Is the weighted fit's stalled log likelihood an unfinished warmup?
#
# What is tested. `_dev/contrast-mixing-cause.R` isolated frequency weights as
# what stops a binomial logit fit mixing: with weights of 1 to 3 on n = 1000,
# the bulk effective sample size of the reported log likelihood over 500
# retained draws is 1 or 2 under either gate, against 92 to 249 without them,
# and under hard rules the predictor and the treatment contrast fall with it,
# to 9 and 3. An effective sample size of 1 is the floor and usually means the
# chain has not moved in that coordinate, which for a log likelihood most often
# means it is still climbing. The claim: the default warmup of 200 sweeps, and
# the 300 used in those runs, is not enough once weights are present, because
# the weights multiply the information and the forest has further to grow.
#
# What is measured. Binomial logit, n = 1000, `y ~ z + x1 + x2 + x3`, weights
# of 1 to 3, 50 trees, one chain, 1000 kept draws, both gates, two
# data-and-fit seeds, with warmup at 300 and at 3000. Per fit:
#
#   ess_loglik    bulk effective sample size of the reported log likelihood
#   drift         its second-half mean minus its first-half mean, in units of
#                 its own standard deviation over the retained draws; a chain
#                 that is still climbing shows a large positive value
#   ess_contrast  bulk effective sample size of the treatment contrast
#   effect        the contrast's posterior mean
#
# What each outcome would mean. If the longer warmup removes the drift and
# restores the effective sample sizes, weights lengthen the burn-in transient,
# the fits above were read before they had converged, and the remedy is warmup,
# which belongs in the documentation of `weights` and in
# `vignette("diagnostics")`. If the drift persists at 3000 sweeps of warmup,
# the chain is not reaching stationarity on a practical budget with weights,
# which is a defect rather than a caution and needs chasing in the engine.
#
# Run with: Rscript _dev/weights-warmup.R
# Writes:   _dev/weights-warmup.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

simulate <- function(seed) {
  n <- 1000L
  set.seed(seed)
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  d$w <- sample(1:3, n, TRUE)
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z
  d$y <- rbinom(n, 1L, plogis(lp))
  d
}

grid <- expand.grid(burn = c(300L, 3000L), gate = c("hard", "smoothstep"),
                    seed = 1:2, stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "Weighted fits: warmup", unit = "fit",
                kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

ess1 <- function(x) posterior::ess_bulk(matrix(x, ncol = 1L))

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  d <- simulate(300L + g$seed)
  d1 <- d
  d1$z <- 1L
  d0 <- d
  d0$z <- 0L

  t0 <- Sys.time()
  set.seed(400L + g$seed)
  fit <- bartisan(y ~ z + x1 + x2 + x3, data = d, family = binomial(),
                  weights = d$w,
                  control = bartisan_control(num_trees = 50L, num_burn = g$burn,
                                             num_draws = 1000L, gate = g$gate,
                                             verbose = FALSE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ll <- fit[["loglik"]]
  half <- seq_len(length(ll) / 2L)
  drift <- (mean(ll[-half]) - mean(ll[half])) / sd(ll)
  contrast <- rowMeans(predict(fit, newdata = d1, type = "link", draws = TRUE) -
                         predict(fit, newdata = d0, type = "link", draws = TRUE))

  rows[[k]] <- data.frame(gate = g$gate, burn = g$burn, seed = g$seed,
                          ess_loglik = ess1(ll), drift = drift,
                          ess_contrast = ess1(contrast),
                          effect = mean(contrast), secs = secs)
  prog_tick(pr, label = sprintf("%s / burn %d / seed %d", g$gate, g$burn, g$seed),
            secs = secs)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/weights-warmup.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- do.call(rbind, lapply(split(res, list(res$gate, res$burn)), function(x) {
  if (!nrow(x)) return(NULL)
  data.frame(gate = x$gate[1L], burn = x$burn[1L],
             ess_loglik = round(median(x$ess_loglik), 1),
             drift_sd = round(median(x$drift), 2),
             ess_contrast = round(median(x$ess_contrast)),
             effect = round(mean(x$effect), 2), secs = round(mean(x$secs), 1))
}))
print(out[order(out$gate, out$burn), ], row.names = FALSE)
