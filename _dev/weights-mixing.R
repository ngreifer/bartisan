# Is the weighted-against-replicated gap real, or chain noise?
#
# What is tested. `_dev/weights-check.R` found that for the binomial logit
# under hard rules the weighted fit and the fit on row-replicated data differ
# by +0.112, +0.184, +0.036 and -0.181 across four fit seeds, and concluded
# from the sign changes that the gap was chain-to-chain variation, so that
# frequency weights behave as documented. `_dev/augmented-mixing.R` then
# measured the same family on the same data without weights and found it mixes
# well: effective sample sizes of 332 to 557 out of 500 retained draws and a
# between-seed spread of the posterior mean of 0.02, which is 14% of a
# posterior standard deviation. If the weighted and replicated fits mix that
# well too, a gap of 0.112 is not noise and that conclusion was wrong. The
# claim being tested is the earlier conclusion itself.
#
# What is measured. The recovery matrix's data, `ybin ~ z + x1 + x2 + x3`,
# binomial logit with its default augmentation, 50 trees, 300 warmup and 500
# kept draws. Four cells: hard and smoothstep rules crossed with the weighted
# fit (n = 1000, weights 1 to 3) and the fit on the same rows replicated by
# those weights (n = 1996). Four fit seeds per cell. The effect is the mean
# over units of link(z = 1) - link(z = 0), averaged over the 1000 original
# rows in both cases so that the two routes report the same quantity. Per cell:
#
#   effect     its posterior mean, averaged over seeds
#   seed_sd    the spread of that mean across the four seeds
#   post_sd    the median posterior standard deviation
#   ess        the median bulk effective sample size of the effect
#
# and per gate the gap between the weighted and replicated cells, read against
# seed_sd.
#
# What each outcome would mean. A between-seed spread that is a small fraction
# of the posterior standard deviation, with effective sample sizes in the
# hundreds, makes the hard-rule gap a real difference between the two routes:
# frequency weights would then not be entering the hard-rule target the way row
# replication does, which is a defect, and `_dev/weights-check.R`'s conclusion
# would be withdrawn. A spread comparable to the gap leaves the conclusion
# standing and makes the mixing of these particular fits the finding.
#
# Run with: Rscript _dev/weights-mixing.R
# Writes:   _dev/weights-mixing.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

library(bartisan)

drawn <- function(k, expr) {
  set.seed(41L + k)
  expr
}
n <- 1000L
d <- drawn(0L, {
  d <- as.data.frame(matrix(runif(n * 3L), n, 3L))
  names(d) <- paste0("x", 1:3)
  d$z <- rbinom(n, 1L, 0.5)
  d$g <- factor(sample(letters[1:8], n, TRUE))
  d$w <- sample(1:3, n, TRUE)
  d
})
u <- drawn(1L, setNames(rnorm(8L, 0, 0.6), letters[1:8]))
f <- 1.2 * d$x1 - 1.2 * d$x2
lp <- f - mean(f) + 0.8 * d$z + u[as.character(d$g)]
d$ybin <- drawn(3L, rbinom(n, 1L, plogis(lp)))

expanded <- d[rep(seq_len(n), d$w), ]
d1 <- d
d1$z <- 1L
d0 <- d
d0$z <- 0L

grid <- expand.grid(gate = c("hard", "smoothstep"),
                    route = c("weighted", "replicated"), seed = 1:4,
                    stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "Weights against replication: mixing",
                unit = "fit", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  control <- bartisan_control(num_trees = 50L, num_burn = 300L, num_draws = 500L,
                              gate = g$gate, verbose = FALSE)
  t0 <- Sys.time()
  set.seed(200L + g$seed)
  fit <- {
    if (identical(g$route, "weighted"))
      bartisan(ybin ~ z + x1 + x2 + x3, data = d, family = binomial(),
               weights = d$w, control = control)
    else
      bartisan(ybin ~ z + x1 + x2 + x3, data = expanded, family = binomial(),
               control = control)
  }
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Averaged over the original rows either way, so the two routes report one
  # quantity rather than two.
  by_draw <- rowMeans(predict(fit, newdata = d1, type = "link", draws = TRUE) -
                        predict(fit, newdata = d0, type = "link", draws = TRUE))

  rows[[k]] <- data.frame(gate = g$gate, route = g$route, seed = g$seed,
                          effect = mean(by_draw),
                          ess = posterior::ess_bulk(matrix(by_draw, ncol = 1L)),
                          post_sd = sd(by_draw), secs = secs)
  prog_tick(pr, label = sprintf("%s / %s / seed %d", g$gate, g$route, g$seed),
            secs = secs)
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          "_dev/weights-mixing.rds")
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
cells <- do.call(rbind, lapply(split(res, list(res$gate, res$route)), function(x) {
  if (!nrow(x)) return(NULL)
  data.frame(gate = x$gate[1L], route = x$route[1L], effect = mean(x$effect),
             seed_sd = sd(x$effect), post_sd = median(x$post_sd),
             ess = median(x$ess))
}))
print(cells[order(cells$gate, cells$route), ], digits = 3, row.names = FALSE)

cat("\nper gate: weighted minus replicated, against the pooled seed spread\n")
for (gt in c("hard", "smoothstep")) {
  a <- cells[cells$gate == gt & cells$route == "weighted", ]
  b <- cells[cells$gate == gt & cells$route == "replicated", ]
  cat(sprintf("  %-10s gap %+.3f   seed_sd %.3f and %.3f   post_sd %.3f\n",
              gt, a$effect - b$effect, a$seed_sd, b$seed_sd, a$post_sd))
}
