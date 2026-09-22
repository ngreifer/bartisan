# Where does the null-cell bias come from, and what removes it?
#
# Cell B0 (their DGP 5 / Setting 3 with the effect set to zero) returns an ATT
# of about +0.18 where the truth is exactly zero. An earlier probe ruled out the
# time trend's functional form -- period dummies did not help -- and found that
# tripling the control forest cut the bias to a third, which says the shortage
# is prognostic capacity rather than shape. Capacity can be bought two ways,
# and this asks which: more trees, or deeper ones. The branching probability is
# gamma * (1 + d)^(-beta), so a smaller `beta` on the control forest alone
# deepens its trees without adding any.
#
# Also asked: whether giving the period a random intercept, rather than a
# numeric column the forest must split on, does the same job more cheaply. That
# one leans on modifiers being allowed to name a variable the control function
# does not carry, since the period then reaches the effect forest and nothing
# else.
#
# Reading: whichever specification takes the bias toward zero without widening
# the interval is the vignette's advice. If depth matches trees, say so, since
# depth is the cheaper of the two. If the random intercept matches both, it is
# cheaper still and is the natural way to write a period effect on panel data.
#
# Usage: Rscript _dev/did-bcf-capacity.R <reps> <out.rds>
suppressPackageStartupMessages(library(bartisan))
setwd("/Users/NoahGreifer/Dropbox/Research/R/bartisan"); source("_dev/did/did-bcf-dgp.R")

args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args)) as.integer(args[1L]) else 12L
OUT <- if (length(args) > 1L) args[2L] else "_dev/did-bcf-capacity.rds"

mods <- paste(c(covs_B, "t", "cohort"), collapse = " + ")
fixed <- paste(c("cohort", "t", covs_B), collapse = " + ")
f_fixed <- stats::as.formula(sprintf("y ~ %s + vc(Dit, ~ %s)", fixed, mods))
# The period leaves the fixed part and becomes a group intercept; it still
# modifies the effect, which only works because a modifier may now name a
# variable the control function does not carry.
f_ranef <- stats::as.formula(sprintf("y ~ %s + vc(Dit, ~ %s) + (1 | t_f)",
                                     paste(c("cohort", covs_B), collapse = " + "), mods))

specs <- list(
  "50/25 trees (reference)"    = list(f = f_fixed, trees = c(50L, 25L)),
  "150/25 trees"               = list(f = f_fixed, trees = c(150L, 25L)),
  "50/25, control beta = 1"    = list(f = f_fixed, trees = c(50L, 25L), beta = c(1, 2)),
  "50/25, control beta = 0.5"  = list(f = f_fixed, trees = c(50L, 25L), beta = c(0.5, 2)),
  "50/25, (1 | period)"        = list(f = f_ranef, trees = c(50L, 25L), ranef = TRUE))

rows <- list()
for (nm in names(specs)) {
  sp <- specs[[nm]]
  for (r in seq_len(REPS)) {
    d <- dgp_B(seed = 5000L + r, tau_mult = 0)
    d$t_f <- factor(d$t)
    tr <- d$Dit == 1
    ctl <- do.call(bartisan_control,
      c(list(num_trees = sp$trees, num_burn = 400L, num_draws = 600L,
             chains = 2L, verbose = FALSE),
        if (!is.null(sp$beta)) list(beta = sp$beta)))
    fit <- suppressMessages(suppressWarnings(
      bartisan(sp$f, data = d, family = gaussian(), control = ctl)))
    att <- rowMeans(coef(fit, draws = TRUE)[[1L]][, tr, drop = FALSE])
    rows[[length(rows) + 1L]] <- data.frame(spec = nm, rep = r,
      att = mean(att), lo = unname(quantile(att, .025)),
      hi = unname(quantile(att, .975)), stringsAsFactors = FALSE)
    saveRDS(list(res = do.call(rbind, rows), complete = FALSE), OUT)
  }
  x <- do.call(rbind, rows); x <- x[x$spec == nm, ]
  cat(sprintf("%-28s bias %+.3f  sd %.3f  rej95 %.2f  width %.3f\n", nm,
      mean(x$att), sd(x$att), mean(x$lo > 0 | x$hi < 0), mean(x$hi - x$lo)))
  utils::flush.console()
}
saveRDS(list(res = do.call(rbind, rows), complete = TRUE), OUT)
