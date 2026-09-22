# Is the Category-B model of did_B_D_and_bart_design.md biased by a group level
# shift its baseline is forbidden to represent?
#
# What is tested: the plan's claim (§3.2-3.3) that a unit random intercept plus
# f_0(X, t) with no access to cohort/ever, and saturated free tau_gt, recovers
# ATT(g,t). A smoke fit returned +3.645 against a truth of 3. dgp_A contains
# 0.75 * ever_i, a pure level shift that f_0 cannot see by construction; the
# suspicion is that the shrunken alpha_i carries only part of it and tau takes
# the rest.
#
# Measured: overall ATT by the plan's own §3.4 formula (mean over treated rows
# of tau_{G,t} + f_tau) with its 95% interval, against tau = 3, over 3
# replicates. Decisive contrast is two cells identical but for the shift:
# (a) as published, (b) 0.75 * ever subtracted from y, same noise draws.
#
# If B is biased ~+0.75 in (a) and unbiased in (b) while the stacked
# long-difference model is unbiased in both, the prohibition on f_0 seeing
# cohort is what makes a group level shift unrepresentable and a random
# intercept too weak to carry it, which is disqualifying for ordinary DiD where
# treated and control differ in level. Biased in both means the diagnosis is
# wrong. Unbiased in both means the smoke fit was an artifact.
suppressPackageStartupMessages(library(bartisan))
SK <- "/Users/NoahGreifer/.claude/skills/live-progress/assets"
source(file.path(SK, "progress.R"))
source("_dev/did/catB-build.R")
setwd("/Users/NoahGreifer/Dropbox/Research/R/bartisan")
source("_dev/did/did-bcf-dgp.R"); source("_dev/did/did-stack-fun.R")

REPS <- 3L
out <- "_dev/did/catB-cmp.rds"
ctl <- function(...) bartisan_control(gate = "hard", num_burn = 250L,
                                      num_draws = 500L, chains = 2L,
                                      verbose = FALSE, ...)
jobs <- list(list(cell="a", model="B-plan"), list(cell="a", model="B-noOverlap"),
             list(cell="a", model="stacked"), list(cell="b", model="B-plan"),
             list(cell="b", model="stacked"))
pr <- prog_init(total = REPS * length(jobs), title = "Category B vs stacked",
                unit = "fit", workers = 1, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)
rows <- list()

for (r in seq_len(REPS)) {
  base <- dgp_A(seed = 4200L + r)
  for (j in jobs) {
    d <- base
    if (j$cell == "b") d$y <- d$y - 0.75 * d$ever      # same noise, no level shift
    res <- try({
      if (j$model == "stacked") {
        s <- stack_cells(d, "id", "t", "y", "G", paste0("x", 1:7))
        f <- dY ~ cell + x1+x2+x3+x4+x5+x6+x7 +
          vc(treat, ~ x1+x2+x3+x4+x5+x6+x7 + cohort + since)
        fit <- bartisan(f, data = s, family = gaussian(),
                        control = ctl(num_trees = c(200L, 10L)))
        a <- rowMeans(coef(fit, draws = TRUE)[["treat"]][, s$treat == 1, drop = FALSE])
      } else {
        dd <- prep_B(d)
        mods <- if (j$model == "B-plan") c("dot","cohort","since") else "dot"
        fit <- bartisan(form_B(dd, mods), data = dd, family = gaussian(),
                        control = ctl(num_trees = 50L))
        a <- as.numeric(att_B(fit, dd, "overall"))
      }
      c(est = mean(a), lo = unname(quantile(a, .025)), hi = unname(quantile(a, .975)))
    }, silent = TRUE)
    ok <- !inherits(res, "try-error")
    rows[[length(rows) + 1L]] <- data.frame(
      rep = r, cell = j$cell, model = j$model,
      est = if (ok) res[["est"]] else NA_real_,
      lo = if (ok) res[["lo"]] else NA_real_,
      hi = if (ok) res[["hi"]] else NA_real_,
      err = if (ok) NA_character_ else conditionMessage(attr(res, "condition")))
    saveRDS(list(res = do.call(rbind, rows), complete = FALSE,
                 done = length(rows), total = REPS * length(jobs)), out)
    prog_tick(pr, label = sprintf("rep %d %s/%s", r, j$cell, j$model))
  }
}
saveRDS(list(res = do.call(rbind, rows), complete = TRUE,
             done = length(rows), total = length(rows)), out)
on.exit(); prog_end(pr, "done")
