# Does Wooldridge's additive baseline fix the Category-B failure on dgp_A?
#
# What is being tested. Phase 14 (DID.md) found the Category-B model biased
# +0.51 with 0/3 coverage on dgp_A as published, because dgp_A carries a group
# level shift (0.75 * ever) that B's baseline was forbidden to see, so only a
# shrunk unit random intercept could carry it. Wooldridge's ETWFE puts the
# cohort in the baseline, additively with time: y ~ cohort + x + vc(one, ~ t +
# x), fit on the untreated rows and imputed onto the treated ones. The claim is
# that this is unbiased on the same draws. It could be false if the cohort
# forest shrinks the level shift too, leaving part of it for the treated rows.
#
# The measurement. Phase 14's seeds (dgp_A(seed = 4200 + r), r = 1..3), its two
# cells ((a) as published, (b) 0.75 * ever subtracted from the same noise
# draws) and its settings (hard gates, 2 chains, 250 + 500 draws). Overall ATT
# and 95% interval against tau = 3, read against phase 14's table: Category B
# +0.510 (0/3) and -0.011 (3/3), stacked -0.292. dgp_A redraws its covariates
# every period, so they enter both forests as time-varying covariates
# (Wooldridge 2025, sec. 10.1).
#
# What each outcome means. Unbiased in both cells: the additive baseline fixes
# the Category-B failure and imputation replaces B as the level-form
# recommendation. Biased in (a) only: forest shrinkage of a cohort level is
# itself the problem and the cohort should enter as free scalars
# (vc(cohort, ~ 1)) rather than through the forest.

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages(library(bartisan))
setwd("/Users/NoahGreifer/Dropbox/Research/R/bartisan")
source("_dev/did/did-bcf-dgp.R")

reps <- 3L
out_file <- Sys.getenv("OUT", "_dev/did/etwfe-dgpA.rds")
xs <- paste0("x", 1:7)
f <- stats::as.formula(paste(
  "y ~ cohort +", paste(xs, collapse = " + "),
  "+ vc(one, ~ t +", paste(xs, collapse = " + "), ") + (1 | id)"))

pr <- prog_init(total = 2L * reps, title = "ETWFE imputation on dgp_A",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)
rows <- list()

for (r in seq_len(reps)) {
  base <- dgp_A(seed = 4200L + r)
  for (cell in c("a", "b")) {
    d <- base
    if (cell == "b") {
      d$y <- d$y - 0.75 * d$ever
    }
    d$one <- 1
    ctl <- d[d$Dit == 0L, ]
    tr <- d[d$Dit == 1L, ]
    set.seed(1L)
    fit <- bartisan(f, data = ctl, family = gaussian(), chains = 2L,
                    gate = "hard", num_trees = c(100L, 50L), num_burn = 250L,
                    num_draws = 500L, verbose = FALSE)
    y0 <- predict(fit, newdata = tr, draws = TRUE)
    att <- rowMeans(sweep(-y0, 2L, tr$y, "+"))
    rows[[length(rows) + 1L]] <- data.frame(
      rep = r, cell = cell, est = mean(att),
      lo = unname(stats::quantile(att, .025)),
      hi = unname(stats::quantile(att, .975)))
    saveRDS(list(res = do.call(rbind, rows), complete = FALSE), out_file)
    prog_tick(pr, label = sprintf("rep %d cell %s", r, cell))
  }
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, complete = TRUE), out_file)
on.exit()
prog_end(pr, "done")

res$bias <- res$est - 3
res$covers <- res$lo <= 3 & 3 <= res$hi
res$width <- res$hi - res$lo
print(res, row.names = FALSE, digits = 3)
print(aggregate(cbind(bias, covers, width) ~ cell, res, mean), digits = 3)
