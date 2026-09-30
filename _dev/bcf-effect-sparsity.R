# Should bcf() turn the sparsity prior on for the effect forest by default?
#
# What is being tested. The package default became `sparsity = FALSE` on
# 2026-09-29 because a variable-selection prior on a predictor whose contrast is
# the estimand puts an atom at zero in the posterior of the effect. In bcf() the
# treatment is the coefficient rather than a predictor the forest splits on, so
# the prior on the effect forest selects among the moderators and cannot drop the
# treatment. The claim at stake: with the control function kept at
# `sparsity = FALSE`, turning the prior on for the effect forest alone
# (`sparsity = c(FALSE, TRUE)`) recovers the conditional effect at least as well
# as the default when there are nuisance moderators, and costs nothing when there
# are none, so it should be bcf()'s default. It could be false if the prior
# over-shrinks real heterogeneity toward a constant effect, or if it costs
# average-effect accuracy or coverage, or mixing. An earlier measurement
# (_dev/TASKS.md, "sparsity in bcf(), now that it can be asked for") found it
# helped on one design; this one varies the presence of nuisance moderators.
#
# The measurement. Gaussian outcome, binary treatment confounded by prognostic
# covariates, n = 500, p in {5, 20} uniform covariates, a five-term prognostic
# function, and three effect surfaces: constant (every moderator is nuisance),
# one linear moderator (x2), and two moderators (x2 and x5). The moderator set
# is every covariate (nuisance present unless p is small and every covariate
# moderates, which happens in no cell here) or, for the two-moderator surface,
# the true two alone (no nuisance). Hard gates, 200 warmup and 500 kept draws,
# the propensity score fit once per dataset and supplied to both settings so
# they differ only in the effect forest's prior. Ten replicates per cell. Read
# `cate_rmse` first (root mean squared error of the conditional effect against
# the true tau(x_i) over the sample), then `ate_bias` and `ate_cover` (the
# mixed average effect, 95% interval), `cate_cover`, `share_true` (the share of
# the effect forest's splitting rules spent on the true moderators), `ess_ate`
# and `seconds`. `atom` is the share of draws with the contrast exactly zero,
# which should be zero everywhere.
#
# What each outcome means. Sparsity on the effect forest lower or equal in
# `cate_rmse` in every cell, with `ate_bias` and both coverages no worse and no
# atom: it becomes bcf()'s default. Worse in the no-nuisance cells (the
# two-moderator surface with the true moderators alone) or in the constant
# cell: the prior over-shrinks real heterogeneity or costs something where it
# has nothing to select, and the default stays off with the option documented.
# Better with nuisance and worse without: the choice depends on the analyst's
# moderator set, which is an argument for leaving it to them and saying so.

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages(library(bartisan))

reps <- as.integer(Sys.getenv("REPS", "10"))
n <- 500L
out_file <- Sys.getenv("OUT", "_dev/bcf-effect-sparsity.rds")

tau_fun <- list(
  none = function(x) rep(1, nrow(x)),
  one  = function(x) 1 + 1.5 * x$x2,
  few  = function(x) 0.5 + 1.5 * x$x2 + 1 * x$x5)
true_mods <- list(none = character(), one = "x2", few = c("x2", "x5"))

cells <- rbind(
  expand.grid(p = c(5L, 20L), tau = c("none", "one", "few"), modset = "all",
              stringsAsFactors = FALSE),
  expand.grid(p = c(5L, 20L), tau = "few", modset = "true",
              stringsAsFactors = FALSE))
settings <- c(off = FALSE, on = TRUE)
total <- reps * nrow(cells) * length(settings)

simulate <- function(p, tau, seed) {
  set.seed(seed)
  x <- as.data.frame(matrix(stats::runif(n * p), n, p,
                            dimnames = list(NULL, paste0("x", seq_len(p)))))
  mu <- 2 * sin(pi * x$x1 * x$x2) + 4 * (x$x3 - 0.5)^2 + 2 * x$x4 + x$x5
  ps <- stats::plogis(-0.5 + 1.5 * x$x1 + 1 * x$x3 - 1 * x$x4)
  x$z <- stats::rbinom(n, 1L, ps)
  tau <- tau_fun[[tau]](x)
  x$y <- mu + tau * x$z + stats::rnorm(n)
  list(data = x, tau = tau)
}

pr <- prog_init(total = total, title = "Sparsity on the bcf() effect forest",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)
rows <- list()

for (r in seq_len(reps)) {
  for (k in seq_len(nrow(cells))) {
    cell <- cells[k, ]
    s <- simulate(cell$p, cell$tau, seed = 1000L * r + k)
    d <- s$data
    covs <- paste0("x", seq_len(cell$p))

    # One propensity score per dataset, shared by both settings.
    ps_fit <- bartisan(stats::reformulate(covs, "z"), data = d,
                       family = stats::binomial(), gate = "hard",
                       num_burn = 200L, num_draws = 300L, verbose = FALSE)
    ps <- as.numeric(stats::fitted(ps_fit))

    mods <- {
      if (cell$modset == "all") NULL
      else stats::reformulate(true_mods[[cell$tau]])
    }

    for (nm in names(settings)) {
      t0 <- Sys.time()
      fit <- bcf(stats::reformulate(covs, "y"), treat = ~ z, data = d,
                 family = stats::gaussian(), gate = "hard",
                 sparsity = c(FALSE, settings[[nm]]), moderators = mods,
                 propensity = ps, num_burn = 200L, num_draws = 500L,
                 verbose = FALSE)
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      ate <- estimate_effect(fit)
      cate <- estimate_effect(fit, estimand = "CATE")
      truth_ate <- mean(s$tau)

      # One vector of draws per contrast; a binary treatment has one contrast.
      ate_draws <- attr(ate, "draws")
      atom <- if (is.null(ate_draws)) NA_real_ else mean(ate_draws[[1L]] == 0)

      ess <- tryCatch({
        tb <- diagnose(ate)$table
        tb$ess_bulk[grepl("-", tb$quantity, fixed = TRUE)][1L]
      }, error = function(e) NA_real_)

      imp <- variable_importance(fit)
      eff <- imp[imp$predictor == "z", ]
      share_true <- {
        if (length(true_mods[[cell$tau]]) == 0L) NA_real_
        else sum(eff$prop_splits[eff$variable %in% true_mods[[cell$tau]]])
      }

      rows[[length(rows) + 1L]] <- data.frame(
        rep = r, p = cell$p, tau = cell$tau, modset = cell$modset,
        sparsity = nm, seconds = round(secs, 2),
        ate_bias = ate$estimate - truth_ate,
        ate_cover = ate$lower <= truth_ate && truth_ate <= ate$upper,
        ate_width = ate$upper - ate$lower,
        cate_rmse = sqrt(mean((cate$estimate - s$tau)^2)),
        cate_cover = mean(cate$lower <= s$tau & s$tau <= cate$upper),
        cate_width = mean(cate$upper - cate$lower),
        cate_sd = stats::sd(cate$estimate),
        share_true = share_true, ess_ate = ess, atom = atom)
      prog_tick(pr, label = sprintf("rep %d p=%d %s/%s %s", r, cell$p,
                                    cell$tau, cell$modset, nm))
    }
  }
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE,
               done = length(rows), total = total), out_file)
}

res <- do.call(rbind, rows)
saveRDS(list(res = res, complete = TRUE, done = nrow(res), total = total),
        out_file)
on.exit()
prog_end(pr, "done")

a <- aggregate(cbind(cate_rmse, cate_cover, cate_width, ate_bias, ate_cover,
                     ate_width, share_true, ess_ate, atom, seconds) ~
                 p + tau + modset + sparsity, res, mean, na.action = na.pass)
a <- a[order(a$p, a$tau, a$modset, a$sparsity), ]
print(a, row.names = FALSE, digits = 3)
