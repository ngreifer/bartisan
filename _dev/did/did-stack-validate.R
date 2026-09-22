# The exact configuration vignettes/causal.Rmd will use, run end to end.
# Confirms (1) the contrast's ESS stays usable at a vignette-affordable number
# of draws now that the effect forest is 10 trees, (2) the three aggregations
# still agree with did, (3) the pre-treatment placebo still covers zero, and
# (4) coef(newdata=)[, "treat"] resolves for the heterogeneity plot. Also times
# both fits, since two BART fits are being added to a vignette CRAN has to build.
suppressPackageStartupMessages({library(bartisan); library(did)})
source("_dev/did/did-stack-fun.R")
data("mpdta", package = "did")
mpdta <- transform(mpdta, first = ifelse(first.treat == 0, Inf, first.treat))
cells <- stack_cells(mpdta, "countyreal", "year", "lemp", "first", "lpop")

ctl <- bartisan_control(num_trees = c(200L, 10L), num_burn = 500L,
                        num_draws = 1000L, chains = 4L, verbose = FALSE)
t0 <- Sys.time(); set.seed(1234)
fit_did <- suppressMessages(bartisan(dY ~ cell + lpop + vc(treat, ~ lpop + cohort + since),
                                     data = cells, family = gaussian(), control = ctl))
t1 <- Sys.time(); cat(sprintf("fit 1: %.0f s (%d rows)\n\n", as.numeric(difftime(t1,t0,units="secs")), nrow(cells)))

cat("== ATT ==\n"); print(estimate_effect(fit_did, treat = "treat", estimand = "ATT"))
cat("\n== by cohort ==\n"); print(as.data.frame(estimate_effect(fit_did, treat="treat", estimand="ATT", by=~cohort))[,1:5], row.names=FALSE)
cat("\n== by since ==\n"); print(as.data.frame(estimate_effect(fit_did, treat="treat", estimand="ATT", by=~since))[,1:5], row.names=FALSE)
cat("\n== diagnose ==\n"); print(diagnose(estimate_effect(fit_did, treat="treat", estimand="ATT"), verbose=FALSE)$table)

cs <- att_gt(yname="lemp", tname="year", idname="countyreal", gname="first.treat",
             xformla=~lpop, data=mpdta, control_group="nevertreated", est_method="reg")
s <- aggte(cs, type="simple", na.rm=TRUE); g <- aggte(cs, type="group", na.rm=TRUE); dy <- aggte(cs, type="dynamic", na.rm=TRUE)
cat(sprintf("\ndid simple %+.4f (se %.4f)\n", s$overall.att, s$overall.se))
cat(sprintf("did group  g%s %+.4f (%.4f)\n", g$egt, g$att.egt, g$se.egt), sep="")
cat(sprintf("did dyn    e=%-2s %+.4f (%.4f)\n", dy$egt, dy$att.egt, dy$se.egt), sep="")

tc <- subset(cells, treat == 1)
tc$tau <- coef(fit_did, newdata = tc)[, "treat"]
cat(sprintf("\ncoef()[, 'treat'] OK: n=%d range %+.3f to %+.3f\n", nrow(tc), min(tc$tau), max(tc$tau)))
ct <- estimate_effect(fit_did, treat="treat", newdata=tc, estimand="CATE")
cat("CATE rows:", nrow(as.data.frame(ct)), "\n")

cells_pre <- stack_cells(mpdta, "countyreal", "year", "lemp", "first", "lpop", pre = TRUE)
t2 <- Sys.time(); set.seed(1234)
fit_pre <- suppressMessages(bartisan(dY ~ cell + lpop + vc(treat, ~ lpop + cohort + since),
                                     data = cells_pre, family = gaussian(), control = ctl))
cat(sprintf("\nfit 2: %.0f s (%d rows)\n", as.numeric(difftime(Sys.time(),t2,units="secs")), nrow(cells_pre)))
cat("\n== placebo (pre) ==\n")
print(as.data.frame(estimate_effect(fit_pre, treat="treat", estimand="ATT", by=~since))[,1:5], row.names=FALSE)
csu <- att_gt(yname="lemp", tname="year", idname="countyreal", gname="first.treat",
              xformla=~lpop, data=mpdta, control_group="nevertreated", est_method="reg",
              base_period="universal")
du <- aggte(csu, type="dynamic", na.rm=TRUE)
cat(sprintf("did universal e=%-2s %+.4f (%.4f)\n", du$egt, du$att.egt, du$se.egt), sep="")
cat(sprintf("\nTOTAL %.0f s\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))
