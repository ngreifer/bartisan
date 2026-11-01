#' Add chains to a fit
#'
#' Runs further chains of the model a fit was made from and adds their draws to
#' the fit, leaving the chains it already has as they are.
#'
#' @param object a `<bartisan_fit>` object, the output of a call to [bartisan()]
#'   or [bcf()].
#' @param chains `numeric`; the number of chains to add. Default is 1. With
#'   \CRANpkg{future.apply} installed and a \CRANpkg{future} plan in place, the
#'   new chains run in parallel, as the chains of the original fit do.
#'
#' @returns
#' An object of the same class as `object`, with the draws of the new chains
#' after those of the existing ones and its call changed to request the total
#' number of chains.
#'
#' @details
#' The new chains are run on the random number streams that a fit with more
#' chains would have used next, so the result is the fit that setting `chains`
#' to the total would have produced from the same seed. For example, a fit made
#' with one chain after `set.seed(123)` and then given three more with
#' `add_chains()` has the same draws as a fit made with `chains = 4` after
#' `set.seed(123)`, and calling `add_chains()` twice with `chains = 1` is the
#' same as calling it once with `chains = 2`. The seed in effect when
#' `add_chains()` is called plays no part. Because the call stored in the result
#' requests the total, calling `update()` on it after the original seed
#' reproduces it as well.
#'
#' The model is rebuilt from the original call, so every setting the fit was
#' made with (e.g., `num_burn`, `num_draws`, and the prior) applies to the new
#' chains, and each new chain runs its own warmup. The data is found the way
#' [`kfold()`][bartisan-interop] finds it: the `data` argument of the original
#' call is evaluated in the environment the formula was created in, so it has to
#' be reachable from there and unchanged since the fit was made, and an error is
#' raised when it is not. Prior weights and an offset are taken from the fit
#' itself. For a [bcf()] fit, the propensity score is not estimated again; the
#' new chains use the score the original chains were given.
#'
#' Adding chains raises the total number of draws, on which the effective sample
#' size depends, so it is a substitute for longer chains when
#' [diagnose()] reports an effective sample size that is too small and R-hat
#' is fine. It costs little wall clock time when there are idle workers to run
#' the new chains on. It does not help when R-hat is the problem: R-hat compares
#' chains against each other, and its value when the chains agree rises with the
#' number of chains being compared, so a fit failing on R-hat needs longer
#' chains, which `update()` with a larger `num_draws` provides. The other common
#' use is a fit made with the default single chain, to which `add_chains()` can
#' add the chains that R-hat needs to compare without running the first chain
#' again.
#'
#' @seealso
#' [diagnose()] for deciding whether more draws are needed;
#' [bartisan_control()] for the `chains` setting itself;
#' `vignette("diagnostics")` for an example and for the other remedies for poor
#' mixing
#'
#' @examples
#' data("rhc")
#' set.seed(123)
#'
#' model <- death ~ rhc + age + sex + aps + meanbp + surv2m
#'
#' # One short chain, as the defaults would give
#' fit <- bartisan(model, data = rhc, family = binomial(),
#'                 num_trees = 10, num_burn = 50,
#'                 num_draws = 50, verbose = FALSE)
#'
#' # Three more, so that R-hat has chains to compare
#' fit4 <- add_chains(fit, chains = 3)
#' fit4
#'
#' # The same draws as four chains from the start
#' set.seed(123)
#' fit4_direct <- bartisan(model, data = rhc, family = binomial(),
#'                         num_trees = 10, num_burn = 50,
#'                         num_draws = 50, chains = 4,
#'                         verbose = FALSE)
#'
#' stopifnot(all.equal(fit4$eta, fit4_direct$eta))
#'
#' diagnose(fit4)
#'
#' @export
add_chains <- function(object, chains = 1L) {
  arg::arg_is(object, "bartisan_fit")
  arg::arg_whole_number(chains)
  arg::arg_gte(chains, 1)

  chains <- as.integer(chains)
  start <- object[["next_stream"]]

  if (is_null(start)) {
    arg::err(c("This fit does not record where its chains' random number
                streams end, so chains cannot be added to it.",
               i = "Fit the model again with this version of {.pkg bartisan}."))
  }

  data <- kfold_data(object, "add_chains")
  rows <- rownames(object[["model"]])

  call <- refit_call(object)
  call[["data"]] <- data[rows, , drop = FALSE]
  call[["chains"]] <- chains

  # The family the fit resolved to, so that a default family is not announced a
  # second time.
  if (!"family" %in% names(call)) {
    call["family"] <- list(object[["family"]])
  }

  # The weights and the offset are taken from the model frame, since either may
  # have lived in the caller's workspace rather than in `data`; see
  # `kfold.bartisan_fit()`, which does the same.
  weights <- stats::model.weights(object[["model"]])

  if (!is_null(weights)) {
    call[["weights"]] <- weights
  }

  supplied <- object[["model"]][["(offset)"]]

  if (!is_null(supplied)) {
    call[["offset"]] <- supplied
  }

  # A fitted propensity score would be fitted again from a stream of its own,
  # and the new chains would then condition on a different covariate from the
  # old ones. The score the fit kept is handed over instead, by row.
  score <- object[["bcf"]][["propensity"]]

  if (!is_null(score)) {
    call[["propensity"]] <- ss(score, match(rows, rownames(data)))
  }

  old <- the$streams
  on.exit(the$streams <- old, add = TRUE)
  the$streams <- stream_sequence(start, chains + 1L)

  more <- eval(call)

  if (!isTRUE(all.equal(as.list(more[["model"]]), as.list(object[["model"]]),
                        check.attributes = FALSE))) {
    arg::err(c("The data this fit was made from is not the data that name now
                reaches, so the new chains would be fitted to different data.",
               i = "{.fun add_chains} refits from the original call, so the data
                    has to be what it was."))
  }

  combine_fits(object, more)
}

# A fit with the chains of `more` after those of `object`. `combine_chains()`
# does the stacking, as it does for the chains of one call; what it leaves to be
# done here is what `bartisan()` does to the engine's output afterwards.
combine_fits <- function(object, more) {
  out <- combine_chains(list(object, more))

  total <- object[["chains"]] + more[["chains"]]

  out[["chains"]] <- total
  out[["control"]][["chains"]] <- total

  # As `merge_control()` would have recorded it had the call asked for the total.
  supplied <- attr(out[["control"]], "supplied") %or% list()
  supplied[["chains"]] <- total
  attr(out[["control"]], "supplied") <- supplied

  out[["next_stream"]] <- more[["next_stream"]]
  out[["loglik"]] <- c(object[["loglik"]], more[["loglik"]])

  for (nm in c("eta", "counts", "ranef", "tau")) {
    if (!is_null(out[[nm]])) {
      names(out[[nm]]) <- names(object[[nm]])
    }
  }

  # The call asks for the total, so that `update()` and `kfold()` refit the model
  # this fit now is.
  out[["call"]][["chains"]] <- total
  out[["call_values"]][["chains"]] <- total

  out[["fitted"]] <- fitted_from_eta(out, vc_combine(out, out[["eta"]], NULL),
                                     average = TRUE)

  out
}
