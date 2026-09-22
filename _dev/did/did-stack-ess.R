# Does ESS really fall as the chain gets longer, or does diagnose() mis-report it?
# ESS is non-decreasing in chain length for a stationary chain, so the observed
# 674 -> 75 across 800 -> 2500 draws is either genuine drift (the nuisance has
# not converged and only a long chain shows it) or a defect in our R-hat/ESS
# code. Measured: diagnose()'s rhat and ess_bulk for the ATT contrast and for
# Y[0], beside posterior::rhat()/ess_bulk() computed independently on the same
# draws reshaped to iterations x chains. Agreement plus a falling ESS means
# drift and the specification is not ready to document; disagreement means the
# bug is ours and gets fixed before the vignette claims anything.
suppressPackageStartupMessages({library(bartisan); library(posterior)})
source("_dev/did/did-stack-fun.R")
data("mpdta", package = "did")
mpdta <- transform(mpdta, first = ifelse(first.treat == 0, Inf, first.treat))
cells <- stack_cells(mpdta, "countyreal", "year", "lemp", "first", "lpop")
tc <- subset(cells, treat == 1)

for (nd in c(800L, 1600L, 2500L)) {
  set.seed(1234)
  fit <- suppressMessages(bartisan(
    dY ~ cell + lpop + vc(treat, ~ lpop + cohort + since), data = cells,
    family = gaussian(),
    control = bartisan_control(num_burn = 500L, num_draws = nd, chains = 4L,
                               verbose = FALSE)))
  e <- estimate_effect(fit, treat = "treat", newdata = tc)
  tab <- diagnose(e, verbose = FALSE)$table
  pd <- attr(e, "po_draws"); ch <- attr(e, "chains")
  cat(sprintf("\nnum_draws=%d  chains=%s  po_draws dims: %s\n", nd, ch,
              paste(sapply(pd, function(z) paste(dim(as.matrix(z)), collapse="x")), collapse=" ")))
  for (q in names(pd)) {
    v <- as.numeric(as.matrix(pd[[q]]))
    m <- matrix(v, ncol = ch)          # iterations x chains
    cat(sprintf("  %-12s diagnose rhat %.3f ess %6.0f | posterior rhat %.3f ess %6.0f  (n=%d)\n",
        q, tab$rhat[match(q, tab$quantity)], tab$ess_bulk[match(q, tab$quantity)],
        posterior::rhat(m), posterior::ess_bulk(m), length(v)))
  }
  cn <- as.numeric(as.matrix(pd[[2]])) - as.numeric(as.matrix(pd[[1]]))
  mc <- matrix(cn, ncol = ch)
  i <- grep("-", tab$quantity)[1]
  cat(sprintf("  %-12s diagnose rhat %.3f ess %6.0f | posterior rhat %.3f ess %6.0f\n",
      "contrast", tab$rhat[i], tab$ess_bulk[i], posterior::rhat(mc), posterior::ess_bulk(mc)))
  utils::flush.console()
}
