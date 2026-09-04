# Progress reporting, delegated to whatever the caller has set up.
#
# The sampler is one long C++ call, so the only way it can report progress is to
# call back into R. What it calls is built here: a function of no arguments that
# signals one step's worth of progress and nothing else. Rendering is *progressr*'s
# business, and the caller's -- this package signals unconditionally and shows
# nothing unless a handler is active, which is what makes the reporting opt-in
# without an argument for it.

# How many times the sampler should report, per chain. Bounded so that reporting
# stays cheap: the callback itself is nothing next to a sweep, but a handler that
# redraws a bar is not, and a tick per sweep would have the reporting cost more
# than the sampling.
PROGRESS_TICKS <- 50L

progress_ticks <- function(control) {
  total <- (control[["num_burn"]] %or% 0L) +
    (control[["num_draws"]] %or% 0L) * (control[["num_thin"]] %or% 1L)

  if (!isTRUE(total > 0)) {
    return(0L)
  }

  as.integer(min(PROGRESS_TICKS, total))
}

# The reporter the engine calls, and how many times it will call it.
#
# `envir` is the frame the progressor belongs to, and it has to be the frame of
# the fit rather than of this function: a progressor finalizes when its frame
# exits, so one owned by this function would be finished before the sampler
# started.
#
# How many times the convergence pass reports, for the same reason the sampler is
# capped: it walks one column per observation, and a handler redrawing a bar
# thousands of times would cost more than the statistics.
PROGRESS_DIAG_TICKS <- 50L

# The reporter a wrapper has claimed for the whole of what a caller asked for.
#
# `bcf()` fits a propensity model and then an outcome model, so two `bartisan()`
# calls happen inside one thing a caller asked for, and each would otherwise
# build a progressor of its own and show a bar of its own. A wrapper claims one
# instead, in `the$claimed_progress`, and both inner fits report into it. It also
# covers the retry in `bcf()`, which calls the outcome fit a second time when a
# drawn coding turns out not to apply.

# One reporter for a call that fits several models in sequence, sized for all of
# them together, so the bar fills once over the whole call rather than once per
# fit. `specs` gives the chain count and the tick count of each fit that will
# actually run, as `progress_reporter()` would have sized them alone.
#
# The tick count is the smallest across the fits, because the same reporter goes
# to all of them and a fit reports at most once per sweep: asking a 20-sweep fit
# for 50 reports would leave the bar 30 short. Taking the minimum instead means
# every chain of every fit spends exactly `ticks` of the total.
#
# The claim is the caller's to release. It takes the frame that owns the bar to
# know when the bar is over, so `bcf()` sets `the$claimed_progress` and drops
# it again in its own `on.exit()`.
shared_reporter <- function(specs, envir = parent.frame()) {
  chains <- sum(vapply(specs, `[[`, numeric(1L), "chains"))
  ticks <- min(vapply(specs, `[[`, numeric(1L), "ticks"))

  if (!isTRUE(chains > 0) || !isTRUE(ticks > 0) ||
        !rlang::is_installed("progressr")) {
    return(list(report = NULL, ticks = 0L))
  }

  p <- progressr::progressor(steps = chains * ticks, envir = envir)

  list(report = function() p(), ticks = as.integer(ticks))
}

# What `progress_reporter()` would size a fit at, without building it: how many
# chains it will run and how many times each reports. `args` is the arguments a
# `bartisan()` call will be given, whose `control` and loose settings merge the
# way that call would merge them.
progress_spec <- function(args) {
  none <- list(chains = 0, ticks = 0)
  control <- args[["control"]] %or% bartisan_control()

  if (!inherits(control, "bartisan_control")) {
    return(none)
  }

  # An unusable setting is the fit's to complain about, with the fit's wording
  # and at the point the caller expects it, so sizing a bar does not raise it
  # first. Nothing is lost by going without a bar for a call that is about to
  # stop anyway.
  control <- tryCatch(merge_control(control, args[names(args) != "control"]),
                      error = function(e) NULL)

  if (is_null(control)) {
    return(none)
  }

  list(chains = control[["chains"]] %or% 1L, ticks = progress_ticks(control))
}

# Sized for every chain at once, so a multi-chain fit fills one bar once instead
# of restarting the count per chain. That the reporter is built here, in the
# calling session, is also what makes it work under *future.apply*: the closure is
# captured by the engine, sent to each worker, and the conditions it signals are
# relayed back as they arrive.
#
# The convergence pass is not part of this, because `bartisan()` no longer runs
# one: it is `diagnose()`'s work now and carries its own reporter.
progress_reporter <- function(chains, control, envir = parent.frame()) {
  # A wrapper may already own the bar for a call that fits more than one model,
  # in which case this fit reports into that one rather than starting a second.
  if (!is_null(the$claimed_progress)) {
    return(the$claimed_progress)
  }

  ticks <- progress_ticks(control)

  if (ticks == 0L || !rlang::is_installed("progressr")) {
    return(list(report = NULL, ticks = 0L))
  }

  p <- progressr::progressor(steps = chains * ticks, envir = envir)

  # A function of no arguments, because one more slice of the run being done is
  # all the sampler has to say.
  report <- function() {
    p()
  }

  list(report = report, ticks = ticks)
}

# The convergence pass's reporter. Separate from the one above because the pass
# is measured in columns rather than in sweeps, and because it is now reached
# through `diagnose()` rather than as part of a fit.
#
# `envir` has to be the caller's frame for the same reason it does above: a
# progressor finishes when its frame exits, so one owned by this function would
# be over before the pass started.
diagnosis_reporter <- function(steps, envir = parent.frame()) {
  if (!isTRUE(steps > 0) || !rlang::is_installed("progressr")) {
    return(NULL)
  }

  p <- progressr::progressor(steps = steps, envir = envir)

  function() {
    p()
  }
}

# Turn `total` units of work into at most `steps` reports, for a phase that runs
# in R rather than in the sampler. Returns a function called once per unit; the
# shortfall from a total that does not divide evenly is left to *progressr*,
# which finishes the bar when the progressor's frame exits.
progress_stepper <- function(report, total, steps) {
  if (is_null(report) || !isTRUE(total > 0) || !isTRUE(steps > 0)) {
    return(function() invisible(NULL))
  }

  every <- max(1L, as.integer(ceiling(total / steps)))
  seen <- 0L
  fired <- 0L

  function() {
    seen <<- seen + 1L

    if (seen %% every == 0L && fired < steps) {
      fired <<- fired + 1L
      report()
    }

    invisible(NULL)
  }
}
