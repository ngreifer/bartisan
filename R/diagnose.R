# Convergence and mixing, in one call.
#
# Everything here is computed from the stored draws, so the report needs no
# package beyond this one. Three statistics do the work and they answer different
# questions, which is what lets the advice at the end be specific rather than a
# list of everything a reader might try:
#
#   * R-hat compares chains against each other. Elevated, it says the chains
#     disagree, and says nothing about why.
#   * The same R-hat on the second half of the draws alone says why. Discarding
#     the early draws is what a longer warmup would have done, so if that fixes
#     it warmup was the problem, and if it does not the chains have settled in
#     different places and more warmup will not help.
#   * Effective sample size says how much information the draws carry. Low with
#     R-hat fine is the benign case and wants only more draws.

#' Check whether a fit converged and mixed
#'
#' Reports the convergence and mixing diagnostics worth looking at before
#' anything is read off a fit, and says what to do about whichever of them fall
#' short. Everything is computed from the stored draws, so no other package is
#' needed.
#'
#' This is where the diagnostics are computed, and the only place: [bartisan()]
#' does not carry a table of its own, because the per-observation statistics cost
#' more than the sampling and would then be recomputed here. With a \CRANpkg{future}
#' plan in place the pass is spread over the workers, and it reports progress
#' through \CRANpkg{progressr} the way the sampler does.
#'
#' @param object a `<bartisan_fit>` object; the output of a call to [bartisan()].
#' @param rhat_max `numeric`; the largest R-hat treated as acceptable. Default
#'   is 1.01, the threshold of Vehtari et al. (2021); 1.1 was the older
#'   convention and is now considered too permissive.
#' @param ess_min `numeric`; the smallest effective sample size treated as
#'   acceptable, for the bulk and the tail alike. Default is 400, Vehtari et
#'   al.'s recommendation of 100 per chain at four chains, which is about what it
#'   takes for the Monte Carlo error of an interval endpoint to be small next to
#'   the posterior's own width.
#'
#' @returns
#' A `<bartisan_diagnosis>` object, a list with a `print()` method that shows the
#' table, the checks and the advice. Its components are
#' \describe{
#'   \item{`table`}{a data frame with one row per quantity: `rhat`, `rhat_late`
#'     (the same statistic on the second half of the draws alone), `ess_bulk`,
#'     `ess_tail`, and `ess_frac`, the bulk effective sample size as a fraction
#'     of the draws kept.}
#'   \item{`checks`}{a data frame of `check`, `status` (`"ok"`, `"warn"` or
#'     `"note"`) and `detail`.}
#'   \item{`advice`}{a character vector, most important first, empty when
#'     everything passed.}
#'   \item{`chains`,`draws`}{how many chains, and how many draws were kept in
#'     total.}
#' }
#'
#' @details
#' ## What Is Reported
#'
#' One row per scalar the sampler draws (the log likelihood, the nuisance
#' parameters of the family, and the scale of each random-effect term), plus one
#' row for the additive predictor and one for each set of group intercepts,
#' summarized over their worst 5% of observations or levels rather than averaged,
#' since an average over a thousand observations hides the ones that have not
#' converged.
#'
#' `rhat` is split-R-hat (Gelman and Rubin, as revised in Gelman et al. 2013):
#' every chain is halved and the halves are compared, so drift inside a chain
#' counts as disagreement rather than hiding inside a chain mean. With one chain
#' it is computed by splitting that chain into segments, which detects drift but
#' cannot detect two chains settling in different places; that is why one chain
#' draws a warning of its own.
#'
#' `rhat_late` is that same statistic computed on the second half of the retained
#' draws alone, and it separates the two reasons chains disagree by running the
#' experiment rather than by testing for it. Discarding the early retained draws
#' is exactly what a longer warmup would have done, so if R-hat is high overall
#' and acceptable late, warmup ended too early. If it stays high late, the chains
#' have each settled somewhere different and a longer warmup will not help.
#'
#' A within-chain drift statistic would answer that question more directly and
#' cannot be made to work at BART's autocorrelation. Three versions were
#' calibrated against stationary autoregressive series, where by construction
#' there is nothing to find: taking each half's Monte Carlo error from that half
#' alone fires 31% of the time at an autocorrelation of 0.995 against a nominal
#' 5%; taking it from the whole chain holds specificity under 3% but then misses
#' a linear trend of six standard deviations three times in four; batch means
#' catch everything and fire 92% of the time on a chain that has converged. A
#' forest is sticky enough to sit where all three fail, so there is no threshold
#' to pick and the statistic is not offered.
#'
#' The rows summarized over observations or levels report the **worst 5%**
#' boundary rather than the single worst column, and the checks are keyed to the
#' *share* of columns that failed rather than to that boundary: the worst of a
#' thousand values is extreme even when every chain has converged, so a threshold
#' applied to a maximum would condemn every fit.
#'
#' `ess_bulk` and `ess_tail` are the rank-normalized effective sample sizes of
#' Vehtari et al. (2021). The tail one is reported separately because a chain can
#' be ample for a posterior mean and nowhere near enough for an interval
#' endpoint.
#'
#' The forest itself is checked too, through the total number of splitting rules
#' in it at each draw. Chains that disagree about how large the forest is are
#' exploring different tree structures, and no generic MCMC diagnostic can see
#' that, because none of them looks at the forest.
#'
#' ## The Leaf Scale Is Left Out
#'
#' `sigma_mu` is deliberately absent from the table. It mixes badly
#' and not for a reason this package can fix. It is still in `fit$sigma_mu` and still reaches
#' [as_draws()][bartisan-interop] for anyone who wants to look.
#'
#' ## What to Do About Poor Mixing
#'
#' The advice the print method gives follows from which statistic failed, and the
#' order matters because the fixes are not interchangeable.
#'
#' **One chain** comes first, because nothing else can be diagnosed properly
#' until there are several. `chains = 4` is the setting to reach for, and with
#' \CRANpkg{future} installed and a parallel backend used, the chains run in parallel, so it usually costs
#' little wall clock.
#'
#' **R-hat elevated but acceptable on the late draws** says warmup ended too
#' early, so `num_burn` is the one to raise; raising `num_draws` instead adds
#' draws from a distribution the sampler has not reached yet. **R-hat elevated on
#' the late draws too** says the chains have each settled somewhere different.
#' Raise `num_burn` and `num_draws` together, and if that does not settle it,
#' reduce `num_trees` (a smaller forest has fewer ways to represent the same fit,
#' so the sampler has less room to wander between them) and check the family,
#' because a likelihood that fits badly can produce a posterior with no single
#' place to be.
#'
#' **An effective sample size that is low while R-hat is fine** is the benign
#' case, and it wants only more draws. Note that thinning does not help:
#' `num_thin` discards draws that were already paid for, so it lowers the
#' effective sample size per unit of time and is worth it only when storing the
#' draws is the binding constraint. **A low tail effective sample size with the
#' bulk fine** says the posterior mean is fine and the interval endpoints are
#' not, so raising `num_draws` is warranted when intervals are what gets
#' reported.
#'
#' @references
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin,
#' D. B. (2013). *Bayesian Data Analysis* (3rd ed.). Chapman and Hall/CRC.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Buerkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: an improved \eqn{\hat{R}} for
#' assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667--718.
#' \doi{10.1214/20-BA1221}
#'
#' @seealso
#' [bartisan_control()] for the settings the advice names;
#' [as_draws()][bartisan-interop] for handing the draws to \CRANpkg{bayesplot}
#' or \CRANpkg{posterior}; `vignette("diagnostics")` for the fuller treatment,
#' including posterior predictive checks
#'
#' @examples
#' data("rhc")
#' set.seed(123)
#'
#' # Two chains, both deliberately short, so that there is something to report
#' fit <- bartisan(death ~ . - days, data = rhc, num_trees = 10, chains = 2,
#'                 num_burn = 50, num_draws = 50, verbose = FALSE)
#'
#' # The table, the checks, and what to do about whichever of them failed
#' diagnose(fit)
#'
#' # A stricter effective sample size, which is what an interval endpoint needs
#' # and a posterior mean does not
#' diagnose(fit, ess_min = 1000)
#'
#' @export
diagnose <- function(object, rhat_max = 1.01, ess_min = 400) {

  if (!inherits(object, "bartisan_fit")) {
    arg::err("{.arg object} must be a fit from {.fn bartisan}")
  }

  arg::arg_number(rhat_max)
  arg::arg_gte(rhat_max, 1)
  arg::arg_number(ess_min)
  arg::arg_gte(ess_min, 1)

  chains <- object[["chains"]] %or% 1L
  draws <- nrow(object[["sigma_mu"]])

  # The pass is the slow part of a multi-chain fit now that `bartisan()` no
  # longer runs one of its own, so it reports progress the way the sampler does:
  # nothing is shown unless a handler is active.
  columns <- diagnosis_column_count(object)
  ticks <- min(PROGRESS_DIAG_TICKS, max(columns, 0))

  # A budget rather than one stepper, because the pass is split by forest and
  # then by chunk, and each piece needs its own share of the bar; see
  # `progress_budget()`.
  budget <- progress_budget(diagnosis_reporter(ticks), columns, ticks)

  table <- diagnosis_table(object, chains, rhat_max, budget)
  checks <- diagnosis_checks(table, chains, draws, rhat_max, ess_min)

  out <- list(table = table,
              checks = checks,
              advice = diagnosis_advice(checks, object[["control"]]),
              chains = chains,
              draws = draws,
              rhat_max = rhat_max,
              ess_min = ess_min)

  class(out) <- "bartisan_diagnosis"
  out
}

# The draws of one quantity as draws by chains, which is the shape every
# statistic below wants. The engine stacks the chains, so the reshape is a fold
# rather than a computation.
as_chains <- function(x, chains) {
  per <- length(x) %/% chains
  matrix(x[seq_len(per * chains)], nrow = per, ncol = chains)
}

# One row per quantity. The reductions over observations and over levels take the
# worst rather than the average, because an average over a thousand observations
# hides the one that has not converged.
# How many columns the pass will walk, which is what its share of a progress bar
# is spread over.
diagnosis_column_count <- function(object) {
  sum(vapply(object[["eta"]], ncol, numeric(1L))) +
    sum(vapply(object[["ranef"]] %or% list(), ncol, numeric(1L)))
}

# Every scalar the sampler draws, flattened to one named vector per quantity.
# Lives here rather than beside the fit because the convergence pass is the only
# thing that reads it: `bartisan()` stopped computing a diagnostics table of its
# own, so this is the one caller.
scalar_draws <- function(object) {
  out <- list(loglik = object[["loglik"]])

  for (h in seq_len(ncol(object[["sigma_mu"]]))) {
    nm <- sprintf("sigma_mu.%s", colnames(object[["sigma_mu"]])[h])
    out[[nm]] <- object[["sigma_mu"]][, h]
  }

  if (!is_null(object[["aux"]])) {
    for (nm in colnames(object[["aux"]])) {
      out[[sprintf("aux.%s", nm)]] <- object[["aux"]][, nm]
    }
  }

  # The scale of each random-effect term is a scalar worth diagnosing; the
  # intercepts themselves are summarized like the predictor, over the worst
  # level, since there is one per level.
  for (h in seq_along(object[["tau"]])) {
    for (r in seq_len(ncol(object[["tau"]][[h]]))) {
      nm <- sprintf("tau.%s.%s", names(object[["tau"]])[h],
                    colnames(object[["tau"]][[h]])[r])
      out[[nm]] <- object[["tau"]][[h]][, r]
    }
  }

  out
}

diagnosis_table <- function(object, chains, rhat_max, budget = NULL) {
  scalars <- scalar_draws(object)

  # Left out for the reason in the documentation: it mixes badly in every
  # implementation and nothing downstream depends on it.
  scalars <- scalars[!startsWith(names(scalars), "sigma_mu.")]

  # The size of the forest at each draw, which is the check a generic diagnostic
  # cannot do: a forest still growing through the retained draws means warmup
  # ended too early, whatever the log likelihood looks like.
  for (h in seq_along(object[["counts"]])) {
    scalars[[sprintf("splits.%s", names(object[["counts"]])[h])]] <-
      rowSums(object[["counts"]][[h]])
  }

  rows <- lapply(names(scalars), function(nm) {
    diagnosis_row(nm, as_chains(scalars[[nm]], chains), rhat_max)
  })

  rows <- c(rows, diagnosis_worst_rows(object, chains, rhat_max, budget))

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# The four numbers, from one quantity's draws by chains.
#
# R-hat and the effective sample sizes need at least two columns to compare, so
# a single chain is folded into its own two halves first -- which is how R-hat is
# defined for one chain, and detects drift even though it cannot detect two
# chains settling in different places. Drift is computed on the unfolded draws,
# because the gap it is looking for is exactly the one folding would hide.
diagnosis_stats <- function(x) {
  wide <- if (ncol(x) < 2L) fold_halves(x) else x

  # R-hat and the bulk effective sample size are both computed from the
  # rank-normalized draws, and ranking them is the single most expensive step in
  # the pass, so it happens once here and both are handed the result. Worth
  # about a twelfth of the pass.
  normalized <- rank_normalize(wide)

  c(rhat = rhat_rank(wide, normalized),
    rhat_late = rhat_late(wide),
    ess_bulk = ess_from(normalized),
    ess_tail = ess_tail(wide))
}

fold_halves <- function(x) {
  half <- nrow(x) %/% 2L

  if (half < 2L) {
    return(x)
  }

  cbind(x[seq_len(half), 1L], x[nrow(x) - half + seq_len(half), 1L])
}

diagnosis_row <- function(quantity, x, rhat_max) {
  stats <- diagnosis_stats(x)

  data.frame(quantity = quantity,
             rhat = stats[["rhat"]],
             rhat_late = stats[["rhat_late"]],
             ess_bulk = stats[["ess_bulk"]],
             ess_tail = stats[["ess_tail"]],
             ess_frac = stats[["ess_bulk"]] / length(x),
             # One component, so the fraction failing is 0 or 1 and the rule
             # below reduces to the plain threshold.
             rhat_bad = as.numeric(isTRUE(stats[["rhat"]] > rhat_max)),
             late_bad = as.numeric(isTRUE(stats[["rhat_late"]] > rhat_max)))
}

# The additive predictor and the group intercepts have one column per observation
# or per level, so each contributes one row summarized over its worst column.
# The statistics for one column of draws each, which is where this whole pass
# spends its time: one rank-normalization and one autocovariance per column, and
# there is one column per observation.
#
# Split over a `future` plan when there is one and there are enough columns to
# pay for the hand-off. The columns are independent and no random numbers are
# drawn, so the chunks are pure functions of their inputs and the result does not
# depend on how many workers ran them; `cut()` gives contiguous chunks in
# ascending order, so `cbind()` puts the columns back where they were. Without
# \CRANpkg{future.apply} it is the same `vapply()` it always was, which is the
# same choice `run_chains()` makes for the chains themselves.
diagnosis_columns <- function(wide, chains, budget) {
  columns <- ncol(wide)

  block <- function(part, step) {
    vapply(seq_len(ncol(part)), function(j) {
      out <- diagnosis_stats(as_chains(part[, j], chains))
      step()
      out
    }, numeric(4L))
  }

  # The hand-off is cheap enough that splitting pays well below the sizes that
  # motivated it: measured on four workers, 200 columns went 3.0x and even 50
  # went 2.1x. The floor is not there because the split stops paying, then, but
  # because below it the whole pass is a few hundredths of a second and the
  # first parallel call in a session still has to start the workers.
  workers <- if (rlang::is_installed("future.apply")) future::nbrOfWorkers() else 1L

  if (columns < 100L || !isTRUE(workers > 1L)) {
    return(block(wide, budget(columns)))
  }

  chunks <- split(seq_len(columns),
                  cut(seq_len(columns), workers, labels = FALSE))

  # One block per worker, each worker slicing its own columns out of `wide`.
  # That makes `wide` a global, which `future` exports to every worker, and by
  # size that looks alarming: 102 MB at 8000 observations, so 819 MB across
  # eight workers. Measured, it costs nothing. Cutting the blocks in this
  # session first, so that only a worker's own columns cross, was within noise
  # of this at every worker count tried (+0.9% at eight, +5.7% at four) and
  # holds a second copy of the draws for the length of the pass. Both forms are
  # bit-identical to the sequential result.
  #
  # What limits the speedup is how many fast cores there are, and nothing in
  # here. Given the same 1000 columns, a worker takes 1.17s alone and 1.97s when
  # eight run at once; per-worker time is nearly flat to four workers and climbs
  # after. That knee is the machine: the M4 this was measured on has four
  # performance cores and six efficiency ones, and a compute-bound loop over
  # 8 KB of data, where memory bandwidth cannot come into it, gives the same
  # curve (1.16x at four workers, 1.43x at six, 1.65x at eight). So expect about
  # 3.5x here and better where there are more than four equal cores.
  # One stepper per chunk, drawn here before anything is dispatched, so that the
  # shares add up to the whole bar however many chunks there are. A single
  # stepper is copied to each worker, and the copies each count from zero, which
  # left the bar short of full by two ticks at four workers and more above that.
  steppers <- lapply(chunks, function(js) budget(length(js)))

  if (rlang::is_installed("future.apply")) {
    parts <- future.apply::future_lapply(
      seq_along(chunks),
      function(k) block(wide[, chunks[[k]], drop = FALSE], steppers[[k]]),
      future.seed = FALSE,
      future.packages = "bartisan")
  }
  else {
    parts <- lapply(seq_along(chunks),
                    function(k) block(wide[, chunks[[k]], drop = FALSE],
                                      steppers[[k]]))
  }

  do.call(cbind, parts)
}

diagnosis_worst_rows <- function(object, chains, rhat_max, budget = NULL) {
  out <- list()

  parts <- list(list(draws = object[["eta"]], stem = "eta", over = "observations"),
                list(draws = object[["ranef"]], stem = "ranef", over = "levels"))

  budget <- budget %or% function(n) function() invisible(NULL)

  for (part in parts) {
    for (h in seq_along(part[["draws"]])) {
      wide <- part[["draws"]][[h]]

      per_column <- diagnosis_columns(wide, chains, budget)

      # The worst 5% boundary rather than the single worst column, because the
      # worst of a thousand values is extreme even when every chain has
      # converged. `high()` and `low()` keep the direction straight: a large
      # R-hat is bad and a small effective sample size is.
      out[[length(out) + 1L]] <- data.frame(
        quantity = sprintf("%s.%s (worst 5%% of %s)", part[["stem"]],
                           names(part[["draws"]])[h], part[["over"]]),
        rhat = high(per_column[1L, ]),
        rhat_late = high(per_column[2L, ]),
        ess_bulk = low(per_column[3L, ]),
        ess_tail = low(per_column[4L, ]),
        ess_frac = low(per_column[3L, ]) / (nrow(wide)),
        # What the checks are keyed to. The percentiles above are for reading;
        # the fraction is what can be thresholded, since the worst of a thousand
        # values is extreme even when every chain has converged.
        rhat_bad = mean(per_column[1L, ] > rhat_max, na.rm = TRUE),
        late_bad = mean(per_column[2L, ] > rhat_max, na.rm = TRUE))
    }
  }

  out
}

# The 95th and 5th percentiles, ignoring the columns where a statistic could not
# be computed. `worst()` handles the all-missing case for the scalar rows; these
# two do the same for the rows summarized over many columns.
high <- function(x) {
  if (!any(is.finite(x))) {
    return(NA_real_)
  }

  stats::quantile(x[is.finite(x)], 0.95, names = FALSE)
}

low <- function(x) {
  if (!any(is.finite(x))) {
    return(NA_real_)
  }

  stats::quantile(x[is.finite(x)], 0.05, names = FALSE)
}

# The same R-hat, computed on the second half of the retained draws.
#
# This is what separates the two reasons chains disagree, and it does so by
# running the experiment rather than by testing for it: if warmup ended too
# early, the early retained draws are the contaminated ones and throwing them
# away is what more burn-in would have done, so R-hat falls. If the chains have
# each settled somewhere different, discarding the early draws changes nothing.
#
# A within-chain drift statistic would answer the same question more directly and
# cannot be made to work here. Three ways of writing one were calibrated against
# stationary autoregressive series: taking each half's Monte Carlo error from
# that half rejects 31% of the time at an autocorrelation of 0.995 against a
# nominal 5%; taking it from the whole chain holds specificity under 3% but then
# misses a linear trend of six standard deviations three quarters of the time;
# batch means catch everything and reject 92% of the time on a stationary chain.
# A BART forest is sticky enough to sit in the region where all three fail, so
# there is no threshold to pick. Re-running a statistic that is already
# calibrated avoids the problem entirely.
rhat_late <- function(x) {
  half <- nrow(x) %/% 2L

  if (half < 4L) {
    return(NA_real_)
  }

  rhat_rank(x[nrow(x) - half + seq_len(half), , drop = FALSE])
}

# The checks, in the order they are worth reading. Each is a statement about the
# fit rather than about a number, because it is the statement the advice below
# is keyed to.
# More than this share of a row's components failing is taken as real. For a
# scalar row the share is 0 or 1, so the rule is the plain threshold; for a row
# over a thousand observations it is the noise floor -- about 5% of them exceed a
# 95% critical value even when every chain is stationary -- with room to spare.
FAIL_SHARE <- 0.2

diagnosis_checks <- function(table, chains, draws, rhat_max, ess_min) {
  rows <- list()

  add <- function(rows, check, status, detail) {
    rows[[length(rows) + 1L]] <- data.frame(check = check, status = status,
                                            detail = detail)

    rows
  }

  if (chains < 2L) {
    rows <- add(rows, "chains", "warn",
                sprintf("one chain, so R-hat can only compare it with itself; %s",
                        "set `chains = 4`"))
  }
  else {
    rows <- add(rows, "chains", "ok",
                sprintf("%d chains, %d draws kept in total", chains,
                        draws))
  }

  worst_at <- function(column, f) {
    v <- table[[column]]

    if (!any(is.finite(v))) {
      return(NULL)
    }

    i <- which(v == f(v, na.rm = TRUE) & is.finite(v))[1L]
    list(value = v[i], quantity = table[["quantity"]][i])
  }

  # Keyed to the share of each row's components that failed, not to the
  # percentile shown in the table; see `FAIL_SHARE`.
  worst_share <- function(column) {
    v <- table[[column]]

    if (!any(is.finite(v))) {
      return(NULL)
    }

    i <- which(v == max(v, na.rm = TRUE) & is.finite(v))[1L]
    list(share = v[i], quantity = table[["quantity"]][i])
  }

  bad_rhat <- worst_share("rhat_bad")

  if (is_null(bad_rhat)) {
    rows <- add(rows, "rhat", "note", "not available")
  }
  else if (bad_rhat[["share"]] > FAIL_SHARE) {
    rows <- add(rows, "rhat", "warn",
                sprintf("above %.2f for %s%s", rhat_max, bad_rhat[["quantity"]],
                        share_suffix(bad_rhat[["share"]])))
  }
  else {
    rows <- add(rows, "rhat", "ok", sprintf("below %.2f throughout", rhat_max))
  }

  bad_late <- worst_share("late_bad")

  if (is_null(bad_late) || is_null(bad_rhat)) {
    rows <- add(rows, "warmup", "note", "not available")
  }
  else if (bad_rhat[["share"]] <= FAIL_SHARE) {
    rows <- add(rows, "warmup", "ok", "long enough, since R-hat is already fine")
  }
  else if (bad_late[["share"]] <= FAIL_SHARE) {
    rows <- add(rows, "warmup", "warn",
                sprintf("too short: R-hat is fine on the second half of the draws alone, which is what more `num_burn` would have given"))
  }
  else {
    rows <- add(rows, "warmup", "note",
                "not the whole story: R-hat stays high on the second half of the draws alone, so the chains disagree rather than merely start badly")
  }

  # The forest's own size gets its own line, because it is the one signal that
  # points at warmup rather than at the number of draws and a reader will not
  # think to look for it.
  forest <- table[startsWith(table[["quantity"]], "splits."), , drop = FALSE]

  if (nrow(forest) > 0L && any(is.finite(forest[["rhat"]]))) {
    at <- which.max(forest[["rhat"]])

    if (isTRUE(forest[["rhat_bad"]][at] > 0)) {
      rows <- add(rows, "forest size", "warn",
                  sprintf("the chains disagree about how many splitting rules the forest has (R-hat %.2f)",
                          forest[["rhat"]][at]))
    }
    else {
      rows <- add(rows, "forest size", "ok",
                  "the chains agree about the size of the forest")
    }
  }

  for (which in c("ess_bulk", "ess_tail")) {
    lo <- worst_at(which, min)
    label <- if (identical(which, "ess_bulk")) "bulk ESS" else "tail ESS"

    if (is_null(lo)) {
      rows <- add(rows, label, "note", "not available")
    }
    else if (lo[["value"]] < ess_min) {
      rows <- add(rows, label, "warn",
                  sprintf("%.0f for %s, below %.0f", lo[["value"]], lo[["quantity"]],
                          ess_min))
    }
    else {
      rows <- add(rows, label, "ok", sprintf("at least %.0f, above %.0f", lo[["value"]],
                                             ess_min))
    }
  }

  lo_frac <- worst_at("ess_frac", min)

  if (!is_null(lo_frac) && lo_frac[["value"]] < 0.05) {
    rows <- add(rows, "autocorrelation", "note",
                sprintf("%s carries %.1f effective draws per hundred kept",
                        lo_frac[["quantity"]], 100 * lo_frac[["value"]]))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ", for 34% of them" when a row has many components and only some failed;
# nothing when the row is one number, where saying "for 100% of it" would be odd.
share_suffix <- function(share) {
  if (isTRUE(all.equal(share, 1))) {
    return("")
  }

  sprintf(", for %.0f%% of them", 100 * share)
}

# What to do, in the order to try it. Keyed to which check failed rather than to
# the numbers, so that the two reasons chains disagree get the two different
# fixes instead of one list of everything.
diagnosis_advice <- function(checks, control = NULL) {
  failed <- function(name) {
    any(checks[["check"]] == name & checks[["status"]] == "warn")
  }

  # What the fit used, so that "raise `num_draws`" names a number the reader can
  # act on rather than sending them back to the call to find out what it was.
  # A clause rather than a sentence, so it reads as an aside where it lands. It
  # comes out empty when the setting is not a single number, which `num_trees`
  # is not when it differs by forest; naming one of several would be worse than
  # naming none.
  had <- function(...) {
    nms <- c(...)
    got <- control[nms]

    if (!all(vapply(got, function(v) length(v) == 1L, logical(1L)))) {
      return("")
    }

    values <- paste(sprintf("`%s`", vapply(got, format, character(1L))),
                    collapse = " and ")

    sprintf(", %s %s", if (length(nms) > 1L) "which were" else "which was",
            values)
  }

  out <- character()

  if (failed("chains")) {
    out <- c(out, paste(
      "Refit with `chains = 4`. R-hat compares chains against each other, and",
      "one chain can only be compared with itself, so nothing below is",
      "reliable until there are several. With *future* installed the chains run",
      "in parallel."))
  }

  warmup <- failed("warmup")

  if (warmup) {
    out <- c(out, paste(
      paste0("Raise `num_burn`", had("num_burn"), "."),
      "R-hat is already acceptable on the second half of the",
      "retained draws on their own, which is what a longer warmup would have",
      "given, so it is the early draws the chains disagree about."))
  }

  if (failed("rhat") && !warmup) {
    out <- c(out, paste(
      paste0("Raise `num_burn` and `num_draws` together",
             had("num_burn", "num_draws"), "."),
      "R-hat stays high even on the",
      "second half of the draws alone, so the chains have each settled",
      "somewhere different rather than merely started badly."))
    out <- c(out, paste(
      paste0("If that does not settle it, reduce `num_trees`",
             had("num_trees"), "."),
      "A smaller forest has",
      "fewer ways to represent the same fit, so the sampler has less room to",
      "move between them."))
    out <- c(out, paste(
      "Then check the family. A likelihood that fits the data badly can give a",
      "posterior with no single place to be; `pp_check()` is the diagnostic."))
  }

  if ((failed("bulk ESS") || failed("tail ESS")) && !failed("rhat")) {
    out <- c(out, paste(
      paste0("Raise `num_draws`", had("num_draws"), "."),
      "The chains agree and are stationary, so they simply",
      "have not run long enough. Do not reach for `num_thin`: thinning",
      "discards draws already paid for and lowers the effective sample size",
      "per unit of time."))
  }

  if (failed("tail ESS") && !failed("bulk ESS")) {
    out <- c(out, paste(
      "The tail is the binding constraint, so a posterior mean is already fine",
      paste0("and an interval endpoint is not. Raise `num_draws`",
             had("num_draws"), " if intervals are"),
      "what gets reported."))
  }

  out
}

#' @export
print.bartisan_diagnosis <- function(x, digits = 3L, ...) {
  cli_cat("{.strong Convergence and mixing}")
  cli::cat_line()

  show <- x[["table"]]
  show[c("ess_frac", "rhat_bad", "late_bad")] <- NULL

  for (nm in c("rhat", "rhat_late")) {
    show[[nm]] <- round(show[[nm]], digits)
  }

  for (nm in c("ess_bulk", "ess_tail")) {
    show[[nm]] <- round(show[[nm]])
  }

  print(show, row.names = FALSE)
  cli::cat_line()

  mark <- c(ok = "v", warn = "x", note = "i")

  for (i in seq_len(nrow(x[["checks"]]))) {
    row <- x[["checks"]][i, ]

    list(row[["detail"]]) |>
      setNames(mark[[row[["status"]]]]) |>
      cli::cli_bullets()
  }

  if (is_null(x[["advice"]])) {
    cli::cat_line()
    cli::cli_alert_success("Nothing to change.")
    return(invisible(x))
  }

  cli::cat_line()
  cli_cat("{.strong What to do}")
  cli::cat_line()

  for (i in seq_along(x[["advice"]])) {
    list(x[["advice"]][i]) |>
      setNames(as.character(i)) |>
      cli::cli_bullets()
  }

  invisible(x)
}
