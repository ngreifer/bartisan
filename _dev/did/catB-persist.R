# Is the Category-B model's bias specific to panels with no unit persistence?
#
# What is tested: B was biased +0.51 with 0/3 coverage on dgp_A, and the cause
# is that a shrunken alpha_i cannot carry the 0.75 * ever level shift. But
# dgp_A redraws its covariates every period, so it has almost no between-unit
# variance (0.105, implying sigma_alpha ~ 0.32) and the random intercept has
# nothing to anchor on. Real panels are the opposite: mpdta has level sd 1.51
# against innovation sd 0.16. If the bias is a shrinkage artifact, adding a
# unit effect makes sigma_alpha large, the RE stops shrinking, and B recovers.
#
# Measured: overall ATT by g-computation against tau = 3, 3 replicates, on
# dgp_A with a unit effect of sd 2 added (shrinkage factor rises from ~0.46 to
# ~0.97, so a bias of +0.51 should fall to roughly +0.02), beside the stacked
# long-difference model, which differences unit effects away and should be
# unmoved from its no-unit-effect numbers.
#
# If B is unbiased here, the verdict is conditional rather than fatal: B wins on
# panels with strong unit persistence and fails on panels without it, and the
# recommendation has to name which. If B stays biased at +0.5, the architecture
# cannot carry a group level shift under any persistence and differencing is
# the only route.
suppressPackageStartupMessages(library(bartisan))
source("_dev/did/catB-build.R")
setwd("/Users/NoahGreifer/Dropbox/Research/R/bartisan")
source("_dev/did/did-bcf-dgp.R"); source("_dev/did/did-stack-fun.R")

gc_att <- function(fit, d) {
  tr <- d$Dit == 1; d0 <- d; d0$Dit <- 0
  for (k in attr(d, "cell_cols")) d0[[k]] <- 0
  z <- rowMeans(predict(fit, newdata = d[tr, ], draws = TRUE) -
                predict(fit, newdata = d0[tr, ], draws = TRUE))
  c(mean(z), quantile(z, .025), quantile(z, .975))
}
ctl <- function(nt) bartisan_control(gate = "hard", num_trees = nt, num_burn = 250L,
                                     num_draws = 500L, chains = 2L, verbose = FALSE)
for (r in 1:3) {
  d <- dgp_A(seed = 4200L + r)
  set.seed(9000L + r)
  d$y <- d$y + rnorm(length(unique(d$id)), 0, 2)[d$id]     # strong unit persistence
  dd <- prep_B(d)
  set.seed(7)
  fb <- bartisan(form_B(dd), data = dd, family = gaussian(), control = ctl(50L))
  a <- gc_att(fb, dd)
  s <- stack_cells(d, "id", "t", "y", "G", paste0("x", 1:7))
  set.seed(7)
  fs <- bartisan(dY ~ cell + x1+x2+x3+x4+x5+x6+x7 +
                   vc(treat, ~ x1+x2+x3+x4+x5+x6+x7 + cohort + since),
                 data = s, family = gaussian(), control = ctl(c(200L, 10L)))
  z <- rowMeans(coef(fs, draws = TRUE)[["treat"]][, s$treat == 1, drop = FALSE])
  b <- c(mean(z), quantile(z, .025), quantile(z, .975))
  cat(sprintf("rep %d  B %+.3f [%+.3f,%+.3f] cov %-5s | stacked %+.3f [%+.3f,%+.3f] cov %s\n",
      r, a[1], a[2], a[3], a[2] <= 3 && 3 <= a[3],
      b[1], b[2], b[3], b[2] <= 3 && 3 <= b[3]))
  utils::flush.console()
}
