# Saved results for the vignettes whose fits take too long to run on CRAN.
#
# A vignette sources this from a hidden setup chunk and sets the chunk option
# `saved = TRUE` on the chunks that follow. Each such chunk is then rendered in
# one of two modes:
#
#   run     The code runs as usual and nothing is saved.
#   replay  The code is not run. The markdown saved in
#           `results/<vignette>/<label>.md` -- the echoed code, the printed
#           output and the links to its figures -- takes the chunk's place, so
#           the page is the one the run that saved it produced, figures
#           included.
#
# The mode comes from the environment variable `BARTISAN_VIGNETTES`, "run" or
# "replay". Unset, it is "run" when `NOT_CRAN` is "true", as `devtools::check()`
# sets it, and "replay" otherwise, which is what CRAN does. Nothing here writes
# a file in either mode.
#
# The saved results are made by code that is not part of the package. The
# script that refreshes them, `_dev/refresh-vignettes.R` in the package's
# repository, supplies it as a function in the option `bartisan.save_results`,
# which is called here with this file's environment in place of either mode.
#
# Replaying a chunk whose header or code has changed since its results were
# saved is an error when `BARTISAN_VIGNETTES` is set to "replay", as
# `_dev/check.sh` sets it, so that a local check catches results that were not
# refreshed. Where the mode is only implied, as on CRAN, it is a warning and the
# saved output is still used: a false alarm there would fail the build, which
# costs more than a stale chunk would.
#
# Inline code that needs a fitted object is wrapped in `saved_value("key",
# expr)`: a run evaluates `expr`, and a replay returns the value saved under
# `key` without evaluating anything.

local({
  name <- sub("[.][Rr]md$", "", basename(knitr::current_input()))
  dir <- file.path("results", name)
  index_file <- file.path(dir, "saved.rds")

  mode <- Sys.getenv("BARTISAN_VIGNETTES")
  strict <- identical(mode, "replay")
  if (!nzchar(mode)) {
    mode <- {
      if (identical(tolower(Sys.getenv("NOT_CRAN")), "true")) "run"
      else "replay"
    }
  }

  if (!mode %in% c("run", "replay")) {
    stop(sprintf("BARTISAN_VIGNETTES must be \"run\" or \"replay\", not \"%s\".",
                 mode), call. = FALSE)
  }

  # The chunk header as written and the code, so that a change to either one
  # invalidates the saved output. The header is known only once the chunk's own
  # hooks run, not yet when the option hooks do, so it is read there.
  #
  # Compared as text, in UTF-8 and without line endings or trailing space, and
  # never through a file: an MD5 of the lines as `writeLines()` wrote them was
  # different on Windows, where a file written in text mode ends its lines in
  # "\r\n", so every chunk saved here read as changed there.
  chunk_source <- function(code, options) {
    x <- enc2utf8(c(options$params.src, code))
    Encoding(x) <- "UTF-8"
    x <- sub("[[:space:]]+$", "", x)
    x[nzchar(x)]
  }

  save_results <- getOption("bartisan.save_results")

  if (is.function(save_results)) {
    save_results(environment())
  }
  else if (identical(mode, "run")) {
    assign("saved_value", function(key, expr) expr, envir = knitr::knit_global())
  }
  else {
    if (!file.exists(index_file)) {
      stop(sprintf("No saved results for vignette '%s' in '%s'. Refresh the saved results first.",
                   name, dir), call. = FALSE)
    }

    saved <- readRDS(index_file)

    first_difference <- function(saved, now) {
      n <- max(length(saved), length(now))
      i <- which(vapply(seq_len(n), function(k) {
        !identical(saved[k], now[k])
      }, logical(1L)))[1L]
      sprintf("line %d, saved \"%s\", now \"%s\"", i, saved[i], now[i])
    }

    read_snippet <- function(file) {
      out <- readChar(file, file.size(file), useBytes = TRUE)
      Encoding(out) <- "UTF-8"
      out
    }

    # The chunk is not run: its code is replaced by code that prints the saved
    # markdown as is, and the chunk hook below hands that markdown on untouched.
    # A chunk that had no output gets none.
    knitr::opts_hooks$set(saved = function(options) {
      if (!isTRUE(options$saved)) {
        return(options)
      }

      label <- options$label

      if (is.null(saved$chunks[[label]])) {
        stop(sprintf("No saved output for chunk '%s' of vignette '%s'. Refresh the saved results.",
                     label, name), call. = FALSE)
      }

      file <- file.path(dir, paste0(label, ".md"))
      options$saved_file <- file
      options$saved_code <- options$code
      options$code <- sprintf("cat(readChar(%s, file.size(%s), useBytes = TRUE))",
                              deparse(file), deparse(file))
      options$eval <- TRUE
      options$echo <- FALSE
      options$results <- "asis"
      options$include <- isTRUE(options$include) && file.size(file) > 0
      options$fig.show <- "hide"
      options$message <- FALSE
      options$warning <- FALSE
      options$error <- FALSE
      options
    })

    # A chunk that is no longer the one that made its saved output is reported
    # before that output is printed.
    knitr::knit_hooks$set(saved = function(before, options, envir) {
      if (!isTRUE(options$saved) || !before) {
        return(NULL)
      }

      label <- options$label
      now <- chunk_source(options$saved_code, options)

      if (!identical(saved$chunks[[label]], now)) {
        msg <- sprintf("The saved output for chunk '%s' of vignette '%s' was made from different code (%s). Refresh the saved results.",
                       label, name, first_difference(saved$chunks[[label]], now))
        if (strict) {
          stop(msg, call. = FALSE)
        }
        warning(msg, call. = FALSE)
      }

      NULL
    })

    chunk_hook <- knitr::knit_hooks$get("chunk")

    knitr::knit_hooks$set(chunk = function(x, options) {
      if (!isTRUE(options$saved)) {
        return(chunk_hook(x, options))
      }

      read_snippet(options$saved_file)
    })

    saved_value <- function(key, expr) {
      if (is.null(saved$values[[key]])) {
        stop(sprintf("No saved value '%s' for vignette '%s'. Refresh the saved results.",
                     key, name), call. = FALSE)
      }

      saved$values[[key]]
    }

    assign("saved_value", saved_value, envir = knitr::knit_global())
  }
})
