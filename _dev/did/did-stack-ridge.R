# Can the cell/(cohort,since) ridge be parameterized away?
# On treated rows `cell` (fixed part) and `cohort`+`since` (modifiers) span the
# same space, so mu and tau trade off: Y[1] mixes well, Y[0] drifts, and the
# contrast inherits ESS ~75 at any chain length. Measured at 2500 draws, which
# is long enough for the drift to show: contrast rhat/ESS and the ATT, against
# did's -0.0420 (se 0.0118). Variant (c) cannot report group or dynamic effects
# and is a diagnostic reference, not a candidate. Hundreds of ESS in (b) or (d)
# means the modifier set is at fault and the vignette uses that one; all near
# 75 means the ridge is intrinsic and gets documented as a limitation.
suppressPackageStartupMessages(library(bartisan))
source("_dev/did/did-stack-fun.R")
data("mpdta", package = "did")
mpdta <- transform(mpdta, first = ifelse(first.treat == 0, Inf, first.treat))
cells <- stack_cells(mpdta, "countyreal", "year", "lemp", "first", "lpop")
tc <- subset(cells, treat == 1)

go <- function(lab, f, ...) {
  set.seed(1234)
  fit <- suppressMessages(bartisan(f, data = cells, family = gaussian(),
    control = bartisan_control(num_burn = 500L, num_draws = 2500L, chains = 4L,
                               verbose = FALSE, ...)))
  e <- estimate_effect(fit, treat = "treat", newdata = tc)
  tb <- diagnose(e, verbose = FALSE)$table; el <- as.data.frame(e)
  i <- grep("-", tb$quantity)[1]; j <- match("Y[0]", tb$quantity)
  cat(sprintf("%-34s ATT %+.4f [%+.4f,%+.4f] | contrast rhat %.3f ess %5.0f | Y0 ess %5.0f\n",
      lab, el$estimate[1], el$lower[1], el$upper[1], tb$rhat[i], tb$ess_bulk[i], tb$ess_bulk[j]))
  utils::flush.console()
}
go("(b) no lpop in tau", dY ~ cell + lpop + vc(treat, ~ cohort + since))
go("(c) no cohort/since in tau", dY ~ cell + lpop + vc(treat, ~ lpop))
go("(d) tau 10 trees", dY ~ cell + lpop + vc(treat, ~ lpop + cohort + since),
   num_trees = c(200L, 10L))
go("(e) cell only in tau", dY ~ cell + lpop + vc(treat, ~ lpop + cell))
