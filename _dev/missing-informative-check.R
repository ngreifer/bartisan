# Which form of test-missing.R's "informative missingness is recovered" check
# passes on nearly every seed and fails when missingness carries nothing?
#
# What is tested. The test compared x1's mean split count with 1.25 times the
# busier noise variable's, from one chain, and failed at seed 60 after 963a3fe
# moved the default dpm() onto a different chain. The claim: some form of "more
# splits on x1 than on noise" passes on nearly all seeds with the signal present
# and on almost none without it.
#
# The measurement. The test's design (n = 400, x1 missing for 35%, y = 2 *
# missing + N(0, .3^2), 20 trees, 300 + 300) on seeds 1 to 200, under the
# default family, with one chain and with four, and an uninformative arm where y
# depends on an independent draw of the same rate. Per fit, x1's mean count and
# the noise variables' maximum and mean. READ: the pass shares of (a) x1 > 1.25
# max(noise) and (b) x1 > 1.25 mean(noise), by arm and chain count.
#
# What each outcome means. A form that passes on at least 99% of informative
# seeds and on few uninformative ones replaces the test's check. If none
# separates them that well, the check is too noisy at this size for a test.
#
# Result (2026-11-05, 963a3fe): (a) one chain .88 / .01, four chains .97 / .00;
# (b) one chain .97 / .03, four chains 1.00 / .00. The test takes (b), four
# chains. See _dev/TASKS.md.
#
# Run against an installed build: Rscript _dev/missing-informative-check.R [workers]
source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
library(bartisan)
args <- commandArgs(TRUE); out <- "_dev/missing-informative-check.rds"
grid <- expand.grid(seed = 1:200, informative = c(TRUE, FALSE), chains = c(1L, 4L))
pr <- prog_init(total = nrow(grid), title = "test-missing: candidate statistics",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)
one <- function(k) {
  g <- grid[k, ]
  set.seed(g$seed); n <- 400
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n), x3 = stats::runif(n))
  gone <- stats::runif(n) < 0.35
  signal <- if (g$informative) gone else stats::runif(n) < 0.35
  d$y <- 2 * signal + stats::rnorm(n, sd = 0.3)
  d$x1[gone] <- NA
  fit <- suppressMessages(bartisan(y ~ ., data = d, na.action = stats::na.pass,
    control = bartisan_control(chains = g$chains, num_trees = 20, num_burn = 300,
                               num_draws = 300, verbose = FALSE)))
  used <- colMeans(fit[["counts"]][["eta"]])
  data.frame(g, x1 = used[["x1"]], max_noise = max(used[c("x2", "x3")]),
             mean_noise = mean(used[c("x2", "x3")]))
}
res <- prog_parallel(pr, seq_len(nrow(grid)), one, workers = if (length(args)) as.integer(args[1]) else NULL,
                     label = function(k) paste(grid$seed[k], grid$informative[k], grid$chains[k]),
                     checkpoint = "_dev/missing-informative-check", packages = "bartisan")
res <- do.call(rbind, res)
saveRDS(list(res = res, complete = TRUE), out)
on.exit(); prog_end(pr, "done")
res$a <- res$x1 > 1.25 * res$max_noise
res$b <- res$x1 > 1.25 * res$mean_noise
res$r <- res$x1 / res$max_noise
print(aggregate(cbind(a, b) ~ informative + chains, res, mean))
print(aggregate(r ~ informative + chains, res, function(x) round(quantile(x, c(.01, .05, .5, .95, .99)), 2)))
