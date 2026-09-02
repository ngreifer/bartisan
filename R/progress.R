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
# Sized for every chain at once, so a multi-chain fit fills one bar once instead
# of restarting the count per chain. That the reporter is built here, in the
# calling session, is also what makes it work under *future.apply*: the closure is
# captured by the engine, sent to each worker, and the conditions it signals are
# relayed back as they arrive.
progress_reporter <- function(chains, control, envir = parent.frame()) {
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
