# Three routes to ATT(g,t) from the Category-B model, checked against each other.
#
# 1. coef(): sum the d_g_t coefficient (tau_gt) and the Dit coefficient (f_tau)
#    on each treated row. Exact, but hand-rolled.
# 2. predict(): g-computation. Zero every treatment-block column and difference
#    the posterior predictions on the response scale. Works for any link, which
#    is what §3.6 requires, and the random intercept cancels.
# 3. estimate_effect(): needs the whole block to hang off ONE variable, so it
#    works only on the "folded" model where tau_gt lives inside f_tau's forest
#    (cohort x since = cell, so the forest can represent it).
#
# If 1 and 2 agree on the plan's model, the manual route is validated. If 3
# agrees with 2 on the folded model, estimate_effect() is the supported path and
# folding costs only shrinkage across cells, not identification.
suppressPackageStartupMessages(library(bartisan))
source("_dev/did/catB-build.R")
setwd("/Users/NoahGreifer/Dropbox/Research/R/bartisan"); source("_dev/did/did-bcf-dgp.R")

gcomp <- function(fit, d, by = NULL) {
  tr <- d$Dit == 1
  d0 <- d; d0$Dit <- 0
  for (k in attr(d, "cell_cols")) d0[[k]] <- 0
  p1 <- predict(fit, newdata = d[tr, ], draws = TRUE)
  p0 <- predict(fit, newdata = d0[tr, ], draws = TRUE)
  g <- p1 - p0
  key <- if (is.null(by)) rep("all", sum(tr)) else as.character(d[[by]][tr])
  sapply(split(seq_len(ncol(g)), key), function(j) {
    z <- rowMeans(g[, j, drop = FALSE]); c(est = mean(z), lo = quantile(z,.025), hi = quantile(z,.975))
  })
}
say <- function(lab, m) cat(sprintf("%-34s %+.3f [%+.3f, %+.3f]\n", lab, m[1], m[2], m[3]))

d <- prep_B(dgp_A(seed = 11L))
ctl <- bartisan_control(gate = "hard", num_trees = 50L, num_burn = 250L,
                        num_draws = 500L, chains = 2L, verbose = FALSE)

cat("=== plan's model (separate tau_gt) ===\n")
set.seed(7); fitB <- bartisan(form_B(d), data = d, family = gaussian(), control = ctl)
a1 <- as.numeric(att_B(fitB, d, "overall"))
say("1. coef() sum", c(mean(a1), quantile(a1,.025), quantile(a1,.975)))
say("2. predict() g-computation", gcomp(fitB, d)[, 1])
cat("3. estimate_effect(): ")
e <- try(as.data.frame(estimate_effect(fitB, treat = "Dit", estimand = "ATT")), silent = TRUE)
if (inherits(e, "try-error")) cat("errors:", conditionMessage(attr(e,"condition")), "\n") else
  cat(sprintf("%+.3f  <- f_tau only, MISSES tau_gt\n", e$estimate[1]))

cat("\n=== folded model (tau_gt inside f_tau) ===\n")
ff <- y ~ x1+x2+x3+x4+x5+x6+x7 + t + vc(tf, ~ 1) +
  vc(Dit, ~ dot_x1+dot_x2+dot_x3+dot_x4+dot_x5+dot_x6 + cohort + sf) + (1 | id)
set.seed(7); fitF <- bartisan(ff, data = d, family = gaussian(),
                              control = bartisan_control(gate = "hard",
                                num_trees = c(50L, 25L), num_burn = 250L,
                                num_draws = 500L, chains = 2L, verbose = FALSE))
dF <- d; attr(dF, "cell_cols") <- character(0)
say("2. predict() g-computation", gcomp(fitF, dF)[, 1])
eF <- as.data.frame(estimate_effect(fitF, treat = "Dit", estimand = "ATT"))
say("3. estimate_effect(ATT)", c(eF$estimate[1], eF$lower[1], eF$upper[1]))
cat("\nby cohort, estimate_effect:\n")
print(as.data.frame(estimate_effect(fitF, treat="Dit", estimand="ATT", by=~cohort))[,1:5], row.names=FALSE)
cat("\nby cohort, g-computation:\n"); print(round(gcomp(fitF, dF, "cohort"), 3))
