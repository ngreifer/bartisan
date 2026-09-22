# Build the Category-B model of did_B_D_and_bart_design.md §3.2 in bartisan:
#
#   Y_it = alpha_i + lambda_t + f_0(X_it, t) + D_it[tau_{G_i,t} + f_tau(Xdot, G, K)]
#
# Mapping to bartisan's formula language:
#   alpha_i        -> (1 | id)                     random intercept
#   lambda_t       -> vc(tf, ~ 1)                  factor -> one drawn scalar per period
#   f_0            -> the fixed part (x1..x7 + t), which IS the control forest
#   tau_{g,t}      -> vc(d_g_t, ~ 1) per treated cell; each d is 0 off its cell
#   f_tau          -> vc(Dit, ~ xdot... [+ cohort + since])
#
# f_0 never receives ever, cohort, or Dit, which is the plan's §3.3 restriction.
suppressPackageStartupMessages(library(bartisan))

prep_B <- function(d) {
  d$tf <- factor(d$t)
  d$since <- ifelse(is.finite(d$G), d$t - d$G, -99)
  d$sf <- factor(d$since)
  xn <- paste0("x", 1:6)                      # x7 is a factor; center the numerics
  for (j in xn) {
    m <- tapply(d[[j]], d$cohort, mean)
    d[[paste0("dot_", j)]] <- ifelse(d$ever == 1, d[[j]] - m[as.character(d$cohort)], 0)
  }
  # one 0/1 column per treated (g, t) cell
  cells <- unique(d[d$Dit == 1, c("G", "t")])
  cells <- cells[order(cells$G, cells$t), ]
  nm <- character(0)
  for (r in seq_len(nrow(cells))) {
    k <- sprintf("d_%s_%s", cells$G[r], cells$t[r])
    d[[k]] <- as.integer(d$G == cells$G[r] & d$t == cells$t[r])
    nm <- c(nm, k)
  }
  attr(d, "cell_cols") <- nm
  attr(d, "cells") <- cells
  d
}

form_B <- function(d, tau_mods = c("dot", "cohort", "since")) {
  nm <- attr(d, "cell_cols")
  taus <- paste(sprintf("vc(%s, ~ 1)", nm), collapse = " + ")
  mods <- character(0)
  if ("dot" %in% tau_mods) mods <- c(mods, paste0("dot_x", 1:6))
  if ("cohort" %in% tau_mods) mods <- c(mods, "cohort")
  if ("since" %in% tau_mods) mods <- c(mods, "sf")
  stats::as.formula(paste0(
    "y ~ ", paste(c(paste0("x", 1:7), "t"), collapse = " + "),
    " + vc(tf, ~ 1)",                       # lambda_t, free per period
    " + ", taus,
    " + vc(Dit, ~ ", paste(mods, collapse = " + "), ")",
    " + (1 | id)"))
}

# ATT(g,t) = mean over cohort-g units of [tau_gt + f_tau(row)].
# coef() gives every coefficient per row per draw, so tau_gt is the d_g_t column
# (constant) and f_tau is the Dit column (per row). No estimate_effect() call can
# reach this, because the block hangs off more than one variable.
att_B <- function(fit, d, by = c("overall", "cohort", "since")) {
  by <- match.arg(by)
  cl <- coef(fit, newdata = d, draws = TRUE)
  tr <- which(d$Dit == 1)
  ftau <- cl[["Dit"]][, tr, drop = FALSE]
  tot <- ftau
  for (k in attr(d, "cell_cols")) {
    on <- d[[k]][tr] == 1
    if (any(on)) tot[, on] <- tot[, on] + cl[[k]][, tr[on], drop = FALSE]
  }
  key <- switch(by, overall = rep("all", length(tr)),
                cohort = as.character(d$cohort[tr]),
                since  = as.character(d$since[tr]))
  out <- lapply(split(seq_along(tr), key), function(j) rowMeans(tot[, j, drop = FALSE]))
  do.call(cbind, out)
}
