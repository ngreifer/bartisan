# Which regime is real data in? mpdta has level sd 1.51 against innovation sd
# 0.16, i.e. strong unit persistence, which the simulation says is where the
# Category-B model's random intercept anchors and its bias disappears. If so, B
# should agree with did on mpdta and be tighter than the vignette's stacked
# model. Measured: overall ATT and the group and dynamic aggregations from the
# FOLDED B model via estimate_effect(), against did::aggte(est_method = "reg")
# and against the vignette's stacked fit (-0.0442 [-0.0697, -0.0194]).
# Agreement with did plus a narrower interval makes B the better choice on real
# panels; disagreement means persistence was not the whole story.
suppressPackageStartupMessages({library(bartisan); library(did)})
data(mpdta, package = "did")
d <- mpdta
d$first <- ifelse(d$first.treat == 0, Inf, d$first.treat)
d$ever  <- as.integer(is.finite(d$first))
d$Dit   <- as.integer(is.finite(d$first) & d$year >= d$first)
d$cohort <- factor(ifelse(is.finite(d$first), as.character(d$first), "never"))
d$yf    <- factor(d$year)
d$sf    <- factor(ifelse(d$Dit == 1, d$year - d$first, -99))
m <- tapply(d$lpop, d$cohort, mean)
d$dot_lpop <- ifelse(d$ever == 1, d$lpop - m[as.character(d$cohort)], 0)

ff <- lemp ~ lpop + year + vc(yf, ~ 1) +
  vc(Dit, ~ dot_lpop + cohort + sf) + (1 | countyreal)
nt <- c(200L, rep(1L, nlevels(d$yf)), 10L)           # intercept, yf levels, Dit
set.seed(1234)
t0 <- Sys.time()
fit <- bartisan(ff, data = d, family = gaussian(),
                control = bartisan_control(num_trees = nt, num_burn = 500L,
                                           num_draws = 1000L, chains = 4L,
                                           verbose = FALSE))
cat(sprintf("fit: %.0f s\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
e <- as.data.frame(estimate_effect(fit, treat = "Dit", estimand = "ATT"))
cat(sprintf("B (folded) overall ATT  %+.4f [%+.4f, %+.4f]  width %.4f\n",
            e$estimate[1], e$lower[1], e$upper[1], e$upper[1] - e$lower[1]))
cat("vignette stacked        -0.0442 [-0.0697, -0.0194]  width 0.0503\n")
cs <- att_gt(yname="lemp", tname="year", idname="countyreal", gname="first.treat",
             xformla=~lpop, data=mpdta, control_group="nevertreated", est_method="reg")
s <- aggte(cs, type="simple", na.rm=TRUE)
cat(sprintf("did simple              %+.4f (se %.4f) -> [%+.4f, %+.4f]\n\n",
            s$overall.att, s$overall.se, s$overall.att-1.96*s$overall.se,
            s$overall.att+1.96*s$overall.se))
cat("by cohort (B):\n")
print(as.data.frame(estimate_effect(fit, treat="Dit", estimand="ATT", by=~cohort))[,c(1,3,4,5)], row.names=FALSE)
g <- aggte(cs, type="group", na.rm=TRUE)
cat(sprintf("did group g%s %+.4f (%.4f)\n", g$egt, g$att.egt, g$se.egt), sep="")
cat("\nby event time (B):\n")
print(as.data.frame(estimate_effect(fit, treat="Dit", estimand="ATT", by=~sf))[,c(1,3,4,5)], row.names=FALSE)
dy <- aggte(cs, type="dynamic", na.rm=TRUE)
cat(sprintf("did dyn e=%-2s %+.4f (%.4f)\n", dy$egt, dy$att.egt, dy$se.egt), sep="")
tb <- diagnose(estimate_effect(fit, treat="Dit", estimand="ATT"), verbose=FALSE)$table
cat("\ndiagnostics:\n"); print(tb[, c("quantity","rhat","ess_bulk")])
cat("\nranef sd (sigma_alpha) mean:", round(mean(unlist(fit$tau)), 3), "\n")
