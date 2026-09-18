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
#' @param object a `<bartisan_fit>` object, the output of a call to [bartisan()],
#'   or a `<bartisan_effect>` object, the output of a call to
#'   [estimate_effect()].
#' @param ... ignored.
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
#' ## Diagnosing an Estimand
#'
#' Called on the output of [estimate_effect()], the table has one row per
#' reported quantity rather than one per sampled parameter, and everything else
#' reads the same way.
#'
#' It is worth doing rather than inferred from the fit's own table, because the
#' two can disagree. An estimand is a contrast, and a contrast can mix badly
#' where the function it is a contrast of mixes well: under the default splitting
#' prior a draw that gives the treatment no rule puts the contrast at exactly
#' zero, and the sampler can stay there for a long run while the other predictors
#' keep the fitted function moving. The check reports the share of draws sitting
#' at the atom when there is one. Chains that have settled on different fitted
#' functions can still agree about an average over them, so this complements the
#' fit's own diagnosis rather than replacing it.
#'
#' @details
#' ## The Rows of the Table
#'
#' One row per scalar the sampler draws (the log likelihood, the nuisance
#' parameters of the family, and the scale of each random-effect term), plus two
#' rows for the additive predictor and two for each set of group intercepts: one
#' summarizing the worst 5% of observations or levels, and one for their average.
#'
#' Both are reported because they routinely disagree, and which one binds depends
#' on what is being reported. A forest settles the level of the fitted function
#' within a sweep or two and takes much longer to settle which observation gets
#' which share of it, so an average or a contrast of averages is governed by the
#' average row and a prediction for one observation by the worst 5%.
#'
#' `rhat` is split-R-hat, so drift inside a chain counts as disagreement rather
#' than hiding inside a chain mean. `rhat_late` is that same statistic on the
#' second half of the retained draws alone, which separates the two reasons
#' chains disagree: if R-hat is high overall and acceptable late, warmup ended
#' too early, and if it stays high late, the chains have each settled somewhere
#' different. `ess_bulk` and `ess_tail` are rank-normalized effective sample
#' sizes, reported separately because a chain can be ample for a posterior mean
#' and nowhere near enough for an interval endpoint. The forest itself is checked
#' through the total number of splitting rules at each draw, since chains that
#' disagree about how large the forest is are exploring different tree
#' structures.
#'
#' The leaf scale `sigma_mu` is left out of the table, and is in `fit$sigma_mu`
#' and [`as_draws()`][bartisan-interop] for anyone who wants to look.
#'
#' ## Setting `rhat_max` and `ess_min`
#'
#' R-hat is a ratio of two variance estimates taken from the same draws, so with
#' few effective draws it sits above 1 whether or not anything is wrong, and how
#' far above depends on how many chains are being compared. `rhat_max` is
#' therefore not a threshold a quantity can be held to at any effective sample
#' size: four chains need about 400 effective draws before 1.01 is even the
#' average of R-hat's null, which is where the pairing of the two defaults comes
#' from, and sixteen chains need about 1600 for the same 1.01. Where a quantity
#' fails R-hat while carrying fewer than that, the checks report that the number
#' cannot be read yet and send the reader to the effective sample size instead.
#'
#' ## Remedies for Poor Mixing
#'
#' The advice the print method gives follows from which statistic failed, and the
#' order matters because the fixes are not interchangeable.
#'
#' One chain comes first, since nothing else can be diagnosed properly until
#' there are several; `chains = 4` is the setting to reach for, and with
#' \CRANpkg{future} installed and a parallel backend in use it usually costs
#' little wall clock. R-hat elevated but acceptable on the late draws says warmup
#' ended too early, so `num_burn` is the one to raise. R-hat elevated on the late
#' draws too says the chains have each settled somewhere different: raise
#' `num_burn` and `num_draws` together, and failing that reduce `num_trees` and
#' check the family, since a likelihood that fits badly can produce a posterior
#' with no single place to be.
#'
#' Effective sample size depends on the total number of draws rather than on how
#' they are divided between chains, where R-hat compares chains against each
#' other, so a fit failing on R-hat wants longer chains rather than more of them
#' and a fit failing only on effective sample size can have either. A low
#' effective sample size with R-hat fine is the benign case and wants only more
#' draws; a low tail effective sample size with the bulk fine says the posterior
#' mean is sound and the interval endpoints are not.
#'
#' `vignette("diagnostics")` works all of this through on a fit.
#'
#' @seealso
#' [bartisan_control()] for the settings the advice names;
#' [`as_draws()`][bartisan-interop] for handing the draws to \CRANpkg{bayesplot}
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
diagnose <- function(object, rhat_max = 1.01, ess_min = 400, ...) {
  UseMethod("diagnose")
}

#' @rdname diagnose
#' @export
diagnose.default <- function(object, rhat_max = 1.01, ess_min = 400, ...) {
  arg::err("{.arg object} must be a fit from {.fn bartisan} or the output of
            {.fn estimate_effect}")
}

#' @rdname diagnose
#' @export
diagnose.bartisan_fit <- function(object, rhat_max = 1.01, ess_min = 400, ...) {

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

#' @rdname diagnose
#' @export
diagnose.bartisan_effect <- function(object, rhat_max = 1.01, ess_min = 400,
                                     ...) {

  arg::arg_number(rhat_max)
  arg::arg_gte(rhat_max, 1)
  arg::arg_number(ess_min)
  arg::arg_gte(ess_min, 1)

  draws <- attr(object, "draws")

  if (is_null(draws)) {
    arg::err(c("This effect carries no posterior draws to diagnose.",
               i = "It was not produced by {.fn estimate_effect}, or the draws
                    were dropped."))
  }

  chains <- attr(object, "chains") %or% 1L

  # A vector per reported quantity, except under `estimand = "CATE"`, where it
  # is one column per unit and the reduction below is the one the fit's own
  # table uses over observations.
  width <- function(d) if (is.matrix(d)) ncol(d) else 1L
  depth <- function(d) if (is.matrix(d)) nrow(d) else length(d)

  n_draws <- depth(draws[[1L]])

  columns <- sum(vapply(draws, width, integer(1L)))
  ticks <- min(PROGRESS_DIAG_TICKS, max(columns, 0))
  budget <- progress_budget(diagnosis_reporter(ticks), columns, ticks)

  rows <- lapply(names(draws), function(nm) {
    d <- draws[[nm]]

    if (!is.matrix(d)) {
      return(diagnosis_row(nm, as_chains(d, chains), rhat_max))
    }

    per_column <- diagnosis_columns(d, chains, budget)

    data.frame(quantity = sprintf("%s (worst 5%% of %d units)", nm, ncol(d)),
               rhat = high(per_column[1L, ]),
               rhat_late = high(per_column[2L, ]),
               ess_bulk = low(per_column[3L, ]),
               ess_tail = low(per_column[4L, ]),
               ess_frac = low(per_column[3L, ]) / nrow(d),
               rhat_bad = mean(per_column[1L, ] > rhat_max, na.rm = TRUE),
               late_bad = mean(per_column[2L, ] > rhat_max, na.rm = TRUE))
  })

  # The two averages the contrast is built from, reported under it. A contrast
  # that mixes badly can have one of them to blame rather than both, and the
  # contrast alone does not say which.
  po <- attr(object, "po_draws")

  if (!is_null(po)) {
    rows <- c(rows, lapply(names(po), function(nm) {
      diagnosis_row(nm, as_chains(po[[nm]], chains), rhat_max)
    }))
  }

  table <- unrowname(do_rbind(rows))
  checks <- diagnosis_checks(table, chains, n_draws, rhat_max, ess_min)
  advice <- diagnosis_advice(checks, attr(object, "control"))

  # A contrast under the sparsity prior has an atom at zero: in a draw where the
  # prior gives the treatment no splitting rule, both potential outcomes are the
  # same number and the contrast is exactly zero. The sampler can stay there for
  # long runs, which costs effective draws on the estimand while leaving the
  # fitted function's own mixing untouched, since the other predictors carry it.
  # So this is the one failure the fit's table cannot show.
  stuck <- vapply(draws, function(d) mean(d == 0, na.rm = TRUE), numeric(1L))

  if (any(stuck > 0.01, na.rm = TRUE)) {
    share <- round(100 * max(stuck, na.rm = TRUE))

    checks <- rbind(
      checks,
      data.frame(check = "atom", status = "note",
                 detail = sprintf(paste("%d%% of draws put the contrast at",
                                        "exactly zero, which is the splitting",
                                        "prior dropping the treatment"),
                                  share)))

    advice <- c(advice, i = paste(
      "Note the atom at zero. The splitting prior drops the treatment in some",
      "draws, and the sampler can stay there for a long run, which costs",
      "effective draws here without costing them in the fit. If the effect is",
      "the quantity being reported, `sparsity = FALSE` removes the atom, and",
      "`bcf()` gives the treatment a forest the prior cannot take it out of;",
      "`vignette(\"causal\")` covers both."))
  }

  out <- list(table = table,
              checks = checks,
              advice = advice,
              chains = chains,
              draws = n_draws,
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

  rows |>
    do_rbind() |>
    unrowname()
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
# The statistics for a block of columns, which is where this whole pass spends
# its time: one rank-normalization and one autocovariance per column, and there
# is one column per observation.
#
# Top level rather than a closure written inside `diagnosis_columns()`, for a
# reason that is invisible until it is measured. A closure carries the frame it
# was written in, and that frame holds `wide`, so handing the closure to a worker
# hands over the whole draw matrix with it however few columns that worker was
# asked for. Written here it carries the namespace and nothing else, and a worker
# receives its own columns and no more.
diagnosis_block <- function(part, chains, step) {
  vapply(seq_len(ncol(part)), function(j) {
    out <- diagnosis_stats(as_chains(part[, j], chains))
    step()
    out
  }, numeric(4L))
}

# Split over a *future* plan when there is one and there are enough columns to
# pay for the hand-off. The columns are independent and no random numbers are
# drawn, so the chunks are pure functions of their inputs and the result does not
# depend on how many workers ran them; `cut()` gives contiguous chunks in
# ascending order, so `cbind()` puts the columns back where they were. Without
# *future* it is the same `vapply()` it always was, which is the same choice
# `run_chains()` makes for the chains themselves.
#
# What limits the speedup is how many fast cores there are, and nothing in here.
# Given the same 1000 columns, a worker takes 1.17s alone and 1.97s when eight
# run at once; per-worker time is nearly flat to four workers and climbs after.
# That knee is the machine: the M4 this was measured on has four performance
# cores and six efficiency ones, and a compute-bound loop over 8 KB of data,
# where memory bandwidth cannot come into it, gives the same curve. So expect
# about 3.5x here and better where there are more than four equal cores.
diagnosis_columns <- function(wide, chains, budget) {
  columns <- ncol(wide)

  # The hand-off is cheap enough that splitting pays well below the sizes that
  # motivated it: measured on four workers, 200 columns went 3.0x and even 50
  # went 2.1x. The floor is not there because the split stops paying, then, but
  # because below it the whole pass is a few hundredths of a second and the
  # first parallel call in a session still has to start the workers.
  workers <- if (rlang::is_installed("future")) future::nbrOfWorkers() else 1L

  # A plan that does not fix a worker count cannot say how many pieces to cut
  # into, and infinitely many is not a question `cut()` can answer.
  if (!isTRUE(is.finite(workers))) {
    workers <- 1L
  }

  if (columns < 100L || !use_future()) {
    return(diagnosis_block(wide, chains, budget(columns)))
  }

  # Cut the columns here rather than on the worker. A worker that slices `wide`
  # itself has to be sent `wide` to slice, and that send is what the bar waits
  # on. One piece per worker is what `future_lapply()` schedules by default
  # anyway, so this decides what travels rather than how the work is spread.
  #
  # Each piece carries its own stepper, drawn before anything is dispatched so
  # that the shares add up to the whole bar however many pieces there are. One
  # stepper copied to every worker would have each copy counting from zero,
  # which left the bar short of full by two ticks at four workers and more above
  # that. The reports themselves come back from the workers as they are made, as
  # they do for the chains of a fit.
  chunks <- split(seq_len(columns),
                  cut(seq_len(columns), workers, labels = FALSE))

  parts <- lapply(chunks, function(js) {
    list(draws = wide[, js, drop = FALSE], step = budget(length(js)))
  })

  # `future_lapply()` rather than a loop over `future()`: it owns the
  # scheduling, and it ships the globals the piece needs without being told,
  # which a bare `future()` call did not. Written as one, \pkg{future} read
  # `diagnosis_block` as belonging to this package, dropped it from what it
  # sent, and left the worker to find an unexported function on a search path
  # that carries only exports.
  blocks <- future.apply::future_lapply(
    parts,
    function(part) diagnosis_block(part[["draws"]], chains, part[["step"]]),
    future.packages = "bartisan", future.seed = FALSE)

  do_cbind(blocks)
}

diagnosis_worst_rows <- function(object, chains, rhat_max, budget = NULL) {
  parts <- list(list(draws = object[["eta"]], stem = "eta", over = "observations"),
                list(draws = object[["ranef"]], stem = "ranef", over = "levels"))

  # Two rows per forest, an average and a worst, over however many forests each
  # part has. A part with nothing in it contributes none.
  out <- vector("list", 2L * sum(lengths(pluck(parts, "draws"))))
  at <- 0L

  budget <- budget %or% function(n) function() invisible(NULL)

  for (part in parts) {
    for (h in seq_along(part[["draws"]])) {
      wide <- part[["draws"]][[h]]

      per_column <- diagnosis_columns(wide, chains, budget)

      # The average first, because the two rows are meant to be read against
      # each other. A forest can disagree with itself about every observation
      # and still agree about their average, and it usually does: the level of
      # the fitted function is settled in a sweep or two, while which
      # observation gets which share of it is the slow part. Which of those is
      # happening decides both what to do about it and whether it matters for
      # what is being reported.
      at <- at + 1L
      out[[at]] <- diagnosis_row(
        sprintf("%s.%s (average over %s)", part[["stem"]],
                names(part[["draws"]])[h], part[["over"]]),
        as_chains(rowMeans(wide), chains), rhat_max)

      # The worst 5% boundary rather than the single worst column, because the
      # worst of a thousand values is extreme even when every chain has
      # converged. `high()` and `low()` keep the direction straight: a large
      # R-hat is bad and a small effective sample size is.
      at <- at + 1L
      out[[at]] <- data.frame(
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

# "bulk ESS" reads as a subject at the head of a sentence only with its first
# letter raised, and the labels are written lowercase because they are also keys.
upper_first <- function(x) {
  sub("^(.)", "\\U\\1", x, perl = TRUE)
}

diagnosis_checks <- function(table, chains, draws, rhat_max, ess_min) {
  rows <- list()

  add <- function(rows, check, status, detail) {
    rows[[length(rows) + 1L]] <- data.frame(check = check, status = status,
                                            detail = detail)

    rows
  }

  if (chains < 2L) {
    rows <- add(rows, "chains", "warn",
                sprintf("Only one chain, so R-hat can only compare it with itself; %s",
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
    list(share = v[i], quantity = table[["quantity"]][i],
         ess = table[["ess_bulk"]][i])
  }

  bad_rhat <- worst_share("rhat_bad")

  if (is_null(bad_rhat)) {
    rows <- add(rows, "rhat", "note", "R-hat is not available")
  }
  else if (bad_rhat[["share"]] > FAIL_SHARE) {
    rows <- add(rows, "rhat", "warn",
                sprintf("R-hat is above %.2f for %s%s", rhat_max,
                        bad_rhat[["quantity"]],
                        share_suffix(bad_rhat[["share"]])))
  }
  else {
    rows <- add(rows, "rhat", "ok",
                sprintf("R-hat is below %.2f throughout", rhat_max))
  }

  # R-hat is a ratio of two variance estimates taken from the same draws, so with
  # few effective draws it sits above 1 whether or not anything is wrong, and how
  # far above depends on how many chains there are. Calibrated against a
  # stationary autoregressive series, where every chain has the same distribution
  # by construction and there is nothing whatever to find, R-hat averages
  # `1 + chains / ess`: at 80 effective draws it is 1.026 over two chains, 1.050
  # over four, 1.103 over eight and 1.205 over sixteen, against 1.025, 1.050,
  # 1.100 and 1.200 from the formula, and it tracks the same way along the other
  # axis, giving 1.162 at 25 effective draws and 1.010 at 400 over four chains.
  #
  # So the threshold a fit has to clear before `rhat_max` means anything is the
  # effective sample size at which R-hat's null reaches it, and that is what is
  # compared here rather than `ess_min`. At the two defaults the two agree, since
  # four chains need 400; at sixteen chains the same 1.01 needs 1600, which is
  # why moving draws into more chains makes R-hat look worse while leaving the
  # information in the draws alone.
  null_rhat <- function(ess) 1 + chains / ess

  unreadable <- !is_null(bad_rhat) && bad_rhat[["share"]] > FAIL_SHARE &&
    isTRUE(null_rhat(bad_rhat[["ess"]]) > rhat_max)

  if (unreadable) {
    rows <- add(rows, "rhat readable", "warn",
                sprintf("That R-hat rests on only %.0f effective draws, where %d chains average %.3f even when they agree",
                        bad_rhat[["ess"]], chains, null_rhat(bad_rhat[["ess"]])))
  }

  bad_late <- worst_share("late_bad")

  if (is_null(bad_late) || is_null(bad_rhat)) {
    rows <- add(rows, "warmup", "note", "Warmup cannot be judged here")
  }
  else if (bad_rhat[["share"]] <= FAIL_SHARE) {
    rows <- add(rows, "warmup", "ok",
                "Warmup was long enough, since R-hat is already fine")
  }
  else if (bad_late[["share"]] <= FAIL_SHARE) {
    rows <- add(rows, "warmup", "warn",
                sprintf("Warmup was too short: R-hat is fine on the second half of the draws alone, which is what more `num_burn` would have given"))
  }
  else if (unreadable) {
    # The stronger reading is withheld here for the same reason the check above
    # exists: R-hat staying high is not evidence of disagreement at this many
    # effective draws. What it does rule out is warmup, which is all this line
    # claims.
    rows <- add(rows, "warmup", "note",
                "A longer warmup is not the fix: R-hat stays high on the second half of the draws alone as well")
  }
  else {
    rows <- add(rows, "warmup", "note",
                "A longer warmup is not the whole story: R-hat stays high on the second half of the draws alone, so the chains disagree rather than merely start badly")
  }

  # The forest's own size gets its own line, because it is the one signal that
  # points at warmup rather than at the number of draws and a reader will not
  # think to look for it.
  forest <- table[startsWith(table[["quantity"]], "splits."), , drop = FALSE]

  if (nrow(forest) > 0L && any(is.finite(forest[["rhat"]]))) {
    at <- which.max(forest[["rhat"]])

    if (isTRUE(forest[["rhat_bad"]][at] > 0)) {
      rows <- add(rows, "forest size", "warn",
                  sprintf("The chains disagree about how many splitting rules the forest has (R-hat %.2f)",
                          forest[["rhat"]][at]))
    }
    else {
      rows <- add(rows, "forest size", "ok",
                  "The chains agree about the size of the forest")
    }
  }

  for (which in c("ess_bulk", "ess_tail")) {
    lo <- worst_at(which, min)
    label <- if (identical(which, "ess_bulk")) "bulk ESS" else "tail ESS"

    if (is_null(lo)) {
      rows <- add(rows, label, "note",
                  sprintf("%s is not available", upper_first(label)))
    }
    else if (lo[["value"]] < ess_min) {
      rows <- add(rows, label, "warn",
                  sprintf("%s is %.0f for %s, below %.0f", upper_first(label),
                          lo[["value"]], lo[["quantity"]], ess_min))
    }
    else {
      rows <- add(rows, label, "ok",
                  sprintf("%s is at least %.0f everywhere, above %.0f",
                          upper_first(label), lo[["value"]], ess_min))
    }
  }

  # Not a fault, so a note: it says which coordinate the trouble is in. The
  # rows summarized over observations come in pairs, and the pair disagreeing is
  # the whole point of reporting both.
  averaged <- table[grepl("(average over", table[["quantity"]], fixed = TRUE), ,
                    drop = FALSE]
  worst <- table[grepl("(worst 5% of", table[["quantity"]], fixed = TRUE), ,
                 drop = FALSE]

  if (nrow(averaged) > 0L && nrow(worst) > 0L &&
        any(worst[["rhat_bad"]] > FAIL_SHARE, na.rm = TRUE) &&
        all(averaged[["rhat_bad"]] == 0, na.rm = TRUE)) {
    rows <- add(rows, "where it is", "note",
                sprintf("The chains disagree about individual observations and agree about their average (R-hat %.2f, %.0f effective draws)",
                        max(averaged[["rhat"]], na.rm = TRUE),
                        min(averaged[["ess_bulk"]], na.rm = TRUE)))
  }

  lo_frac <- worst_at("ess_frac", min)

  if (!is_null(lo_frac) && lo_frac[["value"]] < 0.05) {
    rows <- add(rows, "autocorrelation", "note",
                sprintf("Per-draw efficiency is lowest for %s, which carries %.1f effective draws per hundred kept",
                        lo_frac[["quantity"]], 100 * lo_frac[["value"]]))
  }

  rows |>
    do_rbind() |>
    unrowname()
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

  noted <- function(name) {
    any(checks[["check"]] == name & checks[["status"]] == "note")
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

    # The length check covers a missing `control` as well as a setting that is
    # not one number: subsetting `NULL` gives back nothing rather than a list of
    # nothings, and the clause would otherwise come out as "which was ."
    if (length(got) != length(nms) || !all(lengths(got) == 1L)) {
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
    if (failed("rhat readable")) {
      out <- c(out, paste(
        paste0("Raise `num_draws`", had("num_draws"), "."),
        "R-hat is above the threshold for a quantity that carries too few",
        "effective draws for the threshold to mean anything: with this many",
        "chains it would sit about where it does even if the chains agreed",
        "exactly, as the check above reports. Effective sample size is what",
        "makes it readable, and that grows with the total number of draws;",
        "using fewer chains lowers the bar as well, since R-hat's null rises",
        "with the number of chains being compared."))
    }
    else {
      out <- c(out, paste(
        paste0("Raise `num_burn` and `num_draws` together",
               had("num_burn", "num_draws"), "."),
        "R-hat stays high even on the",
        "second half of the draws alone, so the chains have each settled",
        "somewhere different rather than merely started badly."))
    }

    out <- c(out, paste(
      paste0("If that does not settle it, reduce `num_trees`",
             had("num_trees"), "."),
      "A smaller forest has",
      "fewer ways to represent the same fit, so the sampler has less room to",
      "move between them."))
    out <- c(out, paste(
      "Then check the family. A likelihood that fits the data badly can give a",
      "posterior with no single place to be; `bayesplot::pp_check()` is the",
      "diagnostic."))
  }

  if (noted("where it is")) {
    out <- c(out, paste(
      "Note that the chains disagree about the fitted values of individual",
      "observations and not about their average, which is the usual shape of",
      "this in a forest. What that means for an estimand cannot be read off",
      "this table either way, since an estimand is a contrast and a contrast",
      "can mix badly where the function it contrasts mixes well. Compute it:",
      "`diagnose()` takes the output of `estimate_effect()`, and",
      "`posterior::as_draws()` hands the draws to `posterior::summarise_draws()`",
      "for anything else."))
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

  out |>
    setNames(rep.int("i", length(out)))
}

#' @export
print.bartisan_diagnosis <- function(x, digits = 3L, ...) {
  cli_cat("{.underline Convergence and mixing}")
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
      cli_bullets_cat()
  }

  if (is_null(x[["advice"]])) {
    cli::cat_line()
    cli_bullets_cat(c(v = "Nothing to change."))
    return(invisible(x))
  }

  cli::cat_line()
  cli_cat("{.underline What to do}")
  cli::cat_line()

  # A bullet per step rather than a run-on block. The name has to be one cli
  # recognizes: a number is not, and a numbered name printed as unmarked text,
  # which ran the steps together.
  for (i in seq_along(x[["advice"]])) {
    list(x[["advice"]][i]) |>
      setNames("*") |>
      cli_bullets_cat()
  }

  invisible(x)
}
