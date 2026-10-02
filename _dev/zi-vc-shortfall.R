# Where does the zero-inflated count effect go under a varying coefficient?
#
# What is tested. In `_dev/recovery-matrix.R` the `zi_negbin()` row reads about
# 0.1 lower under `vc()` than in its other structures (0.54 to 0.56 against
# 0.64 to 0.73 for a count-predictor effect of 0.8), and `_dev/recovery-
# followup.R` found the gap falling from 0.11 at n = 1000 to 0.03 at n = 4000,
# which is what weak identification between the count and the zero process
# looks like. There is a second candidate the matrix cannot see. A `vc(z)` term
# on a two-predictor family builds four forests, `count`, `count:z`, `zero` and
# `zero:z`, and the matrix's generator makes a structural zero with a constant
# probability of 0.25, so `zero:z` has nothing real to fit. A coefficient
# forest with nothing to fit can still take rules, and what it takes comes out
# of the count side's share of the zeros. The claim: giving the zero process a
# real dependence on the treatment returns the count effect to its plain-
# structure value.
#
# What is measured. Zero-inflated negative binomial, three uniform predictors,
# `lp = f(x) + 0.8 z` with f centered, count mean `exp(lp - 0.5)` and size 3.
# Two generators: the matrix's, where a structural zero has probability 0.25
# whatever `z` is, and one where its log odds carry the same effect of 0.8, so
# that `zero:z` has signal. Two structures, `y ~ z + x1 + x2 + x3` and
# `y ~ x1 + x2 + x3 + vc(z)`. Sample sizes 1000 and 4000, three seeds, 50
# trees, 300 warmup and 500 kept draws, hard rules. Per fit, from
# `predict(type = "link", draws = TRUE)`:
#
#   effect_count, effect_zero   the mean over units of the contrast in each
#                               predictor between z set to 1 and to 0, read from
#                               the list of one draws-by-observations matrix per
#                               predictor that `predict()` returns here
#   ess_count, ess_zero         their bulk effective sample sizes
#
# The count effect is 0.8 under both generators; the zero effect is 0 under the
# first and 0.8 under the second.
#
# What each outcome would mean. A count effect that is short under `vc()` with
# the constant zero process and recovers with the dependent one puts the
# shortfall in the idle `zero:z` forest, and the matrix's dip is an artifact of
# a generator that leaves one forest nothing to do rather than a property of
# the family. A shortfall under both generators is the weak identification
# already measured, which n reduces and which nothing here would change.
#
# Run with: Rscript _dev/zi-vc-shortfall.R
# Writes:   _dev/zi-vc-shortfall.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

simulate <- function(n, seed, zero_on_z) {
  set.seed(seed)
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  f <- 1.2 * d$x1 - 1.2 * d$x2
  lp <- f - mean(f) + 0.8 * d$z
  # A structural zero with probability 0.25 at z = 0 either way; under the
  # second generator its log odds rise by 0.8 with the treatment.
  p_zero <- if (zero_on_z) plogis(qlogis(0.25) + 0.8 * d$z) else rep(0.25, n)
  d$y <- rnbinom(n, mu = exp(lp - 0.5), size = 3) * rbinom(n, 1L, 1 - p_zero)
  d
}

grid <- expand.grid(n = c(1000L, 4000L), zero_on_z = c(FALSE, TRUE),
                    structure = c("plain", "vc"), seed = 1:3,
                    stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "Zero-inflated effect under vc()",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

ess1 <- function(x) posterior::ess_bulk(matrix(x, ncol = 1L))

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  d <- simulate(g$n, 500L + g$seed, g$zero_on_z)
  d1 <- d
  d1$z <- 1L
  d0 <- d
  d0$z <- 0L

  form <- {
    if (identical(g$structure, "plain")) y ~ z + x1 + x2 + x3
    else y ~ x1 + x2 + x3 + vc(z)
  }
  environment(form) <- globalenv()

  t0 <- Sys.time()
  set.seed(600L + g$seed)
  fit <- bartisan(form, data = d, family = zi_negbin(),
                  control = bartisan_control(num_trees = 50L, num_burn = 300L,
                                             num_draws = 500L, gate = "hard",
                                             verbose = FALSE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  e1 <- predict(fit, newdata = d1, type = "link", draws = TRUE)
  e0 <- predict(fit, newdata = d0, type = "link", draws = TRUE)
  # A list with one draws-by-observations matrix per additive predictor; the
  # first is the count predictor and the second the zero one.
  count <- rowMeans(e1[[1L]] - e0[[1L]])
  zero <- rowMeans(e1[[2L]] - e0[[2L]])

  rows[[k]] <- data.frame(n = g$n, zero_on_z = g$zero_on_z,
                          structure = g$structure, seed = g$seed,
                          effect_count = mean(count), effect_zero = mean(zero),
                          ess_count = ess1(count), ess_zero = ess1(zero),
                          secs = secs)
  prog_tick(pr, label = sprintf("n %d / zero_on_z %s / %s / seed %d", g$n,
                                g$zero_on_z, g$structure, g$seed), secs = secs)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/zi-vc-shortfall.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- do.call(rbind, lapply(split(res, list(res$n, res$zero_on_z, res$structure)),
  function(x) {
    if (!nrow(x)) return(NULL)
    data.frame(n = x$n[1L], zero_on_z = x$zero_on_z[1L],
               structure = x$structure[1L],
               count = round(mean(x$effect_count), 3),
               zero = round(mean(x$effect_zero), 3),
               ess_count = round(median(x$ess_count)),
               ess_zero = round(median(x$ess_zero)))
  }))
out <- out[order(out$n, out$zero_on_z, out$structure), ]
print(out, row.names = FALSE)
