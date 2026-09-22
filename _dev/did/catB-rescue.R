# Can the Category-B model be rescued by giving f_0 a pre-treatment outcome?
#
# What is tested: the model is biased on dgp_A because 0.75 * ever_i is a unit
# level shift that f_0 is forbidden to represent and a shrunken alpha_i carries
# only part of. A pre-treatment mean outcome is a legitimate baseline covariate
# (every unit is untreated at t <= 3, since the earliest cohort is 4), it is not
# affected by treatment, and it proxies the unit level without naming cohort.
# Adding it to f_0's inputs is a one-line change that should let the baseline
# absorb the shift.
#
# Measured: overall ATT by g-computation against tau = 3 over 3 replicates, for
# the plan's model as written and the same model with ypre in the fixed part.
#
# If ypre removes the bias, the defect is the missing unit-level information
# rather than the architecture, and the fix is available -- but it must be
# reported with the caveat that conditioning on a pre-treatment outcome changes
# the identifying assumption away from parallel trends (regression to the mean,
# Ashenfelter dip), so it is not free. If the bias survives, the architecture
# cannot represent a group level shift at all and differencing is the only route.
suppressPackageStartupMessages(library(bartisan))
source("_dev/did/catB-build.R")
setwd("/Users/NoahGreifer/Dropbox/Research/R/bartisan"); source("_dev/did/did-bcf-dgp.R")

gc_att <- function(fit, d) {
  tr <- d$Dit == 1; d0 <- d; d0$Dit <- 0
  for (k in attr(d, "cell_cols")) d0[[k]] <- 0
  g <- predict(fit, newdata = d[tr, ], draws = TRUE) -
       predict(fit, newdata = d0[tr, ], draws = TRUE)
  z <- rowMeans(g); c(mean(z), quantile(z, .025), quantile(z, .975))
}
ctl <- function(nt) bartisan_control(gate = "hard", num_trees = nt, num_burn = 250L,
                                     num_draws = 500L, chains = 2L, verbose = FALSE)
for (r in 1:3) {
  d <- prep_B(dgp_A(seed = 4200L + r))
  yp <- tapply(d$y[d$t <= 3], d$id[d$t <= 3], mean)      # untreated for everyone
  d$ypre <- as.numeric(yp[as.character(d$id)])
  f0 <- form_B(d)
  f1 <- stats::update(f0, . ~ . + ypre)
  for (lab in c("plan", "plan + ypre")) {
    set.seed(7)
    fit <- bartisan(if (lab == "plan") f0 else f1, data = d,
                    family = gaussian(), control = ctl(50L))
    a <- gc_att(fit, d)
    cat(sprintf("rep %d  %-12s ATT %+.3f [%+.3f, %+.3f]  covers 3: %s\n",
                r, lab, a[1], a[2], a[3], a[2] <= 3 && 3 <= a[3]))
    utils::flush.console()
  }
}
