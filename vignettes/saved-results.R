# Saved results for the vignettes whose fits take too long to run on CRAN.
#
# A vignette sources this from a hidden setup chunk and sets the chunk option
# `saved = TRUE` on the chunks that follow. Each such chunk is then rendered in
# one of two modes:
#
#   run     The code runs as usual. The chunk's finished markdown -- the echoed
#           code, the printed output and the links to its figures -- is written
#           to `results/<vignette>/<label>.md`, its figures are copied beside it,
#           and the chunk's header and code are recorded.
#   replay  The code is not run. The saved markdown takes the chunk's place, so
#           the page is the one the run produced, figures included.
#
# The mode comes from the environment variable `BARTISAN_VIGNETTES`, "run" or
# "replay". Unset, it is "run" when `NOT_CRAN` is "true", as `devtools::check()`
# sets it, and "replay" otherwise, which is what CRAN does.
#
# `R CMD build` and `R CMD check` knit a temporary copy of the package, so a run
# there checks that the code works but its results are thrown away with the
# copy. To refresh the saved results, render in place: `_dev/refresh-vignettes.R`
# does, and checks that a replay reproduces the run's page byte for byte.
#
# Replaying a chunk whose header or code has changed since its results were
# saved is an error when `BARTISAN_VIGNETTES` is set to "replay", as
# `_dev/check.sh` sets it, so that a local check catches results that were not
# refreshed. Where the mode is only implied, as on CRAN, it is a warning and the
# saved output is still used: a false alarm there would fail the build, which
# costs more than a stale chunk would.
#
# Inline code that needs a fitted object is wrapped in `saved_value("key",
# expr)`: a run evaluates `expr` and stores the value under `key`, and a replay
# returns the stored value without evaluating anything.

local({
  name <- sub("[.][Rr]md$", "", basename(knitr::current_input()))
  dir <- file.path("results", name)
  index_file <- file.path(dir, "saved.rds")

  mode <- Sys.getenv("BARTISAN_VIGNETTES")
  run <- {
    if (nzchar(mode)) identical(mode, "run")
    else identical(tolower(Sys.getenv("NOT_CRAN")), "true")
  }
  strict <- identical(mode, "replay")

  # Each chunk's header and code, and the inline values. A run starts from
  # nothing, so that whatever a chunk no longer produces is not carried over.
  state <- new.env()
  state$saved <- {
    if (!run && file.exists(index_file)) readRDS(index_file)
    else list(chunks = list(), values = list())
  }
  state$written <- character()

  if (run) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  else if (!file.exists(index_file)) {
    stop(sprintf("No saved results for vignette '%s' in '%s'. Render it with BARTISAN_VIGNETTES=run first.",
                 name, dir), call. = FALSE)
  }

  # The chunk header as written and the code, so that a change to either one
  # invalidates the saved output. The header is known only once the chunk's own
  # hooks run, not yet when the option hooks do, so both modes read it there.
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

  first_difference <- function(saved, now) {
    n <- max(length(saved), length(now))
    i <- which(vapply(seq_len(n), function(k) {
      !identical(saved[k], now[k])
    }, logical(1L)))[1L]
    sprintf("line %d, saved \"%s\", now \"%s\"", i, saved[i], now[i])
  }

  save_index <- function() {
    saveRDS(state$saved, index_file)
  }

  read_snippet <- function(file) {
    out <- readChar(file, file.size(file), useBytes = TRUE)
    Encoding(out) <- "UTF-8"
    out
  }

  # A chunk's figures are written under `fig.path`, which rmarkdown deletes
  # after rendering. They are copied into the results directory and the
  # markdown is pointed at the copies, in a run and so in its replay alike.
  keep_figures <- function(out, options) {
    prefix <- options$fig.path
    if (!nzchar(prefix) || !grepl(prefix, out, fixed = TRUE)) {
      return(out)
    }

    # Each mention of the prefix is followed by a file name, which runs to the
    # bracket, quote or space that closes the link.
    after <- strsplit(out, prefix, fixed = TRUE)[[1L]][-1L]
    files <- unique(paste0(prefix, sub("[)\"' \n].*", "", after)))

    for (f in files) {
      to <- file.path(dir, basename(f))
      file.copy(f, to, overwrite = TRUE)
      state$written <- c(state$written, basename(to))
    }

    gsub(prefix, paste0(dir, "/"), out, fixed = TRUE)
  }

  # In a replay the chunk is not run: its code is replaced by code that prints
  # the saved markdown as is, and the chunk hook below hands that markdown on
  # untouched. A chunk that had no output gets none.
  knitr::opts_hooks$set(saved = function(options) {
    if (!isTRUE(options$saved) || run) {
      return(options)
    }

    label <- options$label

    if (is.null(state$saved$chunks[[label]])) {
      stop(sprintf("No saved output for chunk '%s' of vignette '%s'. Render it with BARTISAN_VIGNETTES=run.",
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

  # In a run, every chunk with the option set is recorded with its source and
  # an empty snippet, which the chunk hook overwrites when there is output; it
  # is called after this. In a replay, a chunk that is no longer the one that
  # made its saved output is reported before that output is printed.
  knitr::knit_hooks$set(saved = function(before, options, envir) {
    if (!isTRUE(options$saved)) {
      return(NULL)
    }

    label <- options$label

    if (run && !before) {
      file <- file.path(dir, paste0(label, ".md"))
      file.create(file)
      state$written <- c(state$written, basename(file))
      state$saved$chunks[[label]] <- chunk_source(options$code, options)
      save_index()
    }
    else if (!run && before) {
      saved <- state$saved$chunks[[label]]
      now <- chunk_source(options$saved_code, options)

      if (!identical(saved, now)) {
        msg <- sprintf("The saved output for chunk '%s' of vignette '%s' was made from different code (%s). Render it with BARTISAN_VIGNETTES=run.",
                       label, name, first_difference(saved, now))
        if (strict) {
          stop(msg, call. = FALSE)
        }
        warning(msg, call. = FALSE)
      }
    }

    NULL
  })

  chunk_hook <- knitr::knit_hooks$get("chunk")

  knitr::knit_hooks$set(chunk = function(x, options) {
    if (!isTRUE(options$saved)) {
      return(chunk_hook(x, options))
    }

    if (!run) {
      return(read_snippet(options$saved_file))
    }

    # knitr passes an excluded chunk through here too and drops it from the
    # page afterwards, so its snippet stays empty.
    out <- keep_figures(chunk_hook(x, options), options)
    if (isTRUE(options$include)) {
      writeChar(out, file.path(dir, paste0(options$label, ".md")), eos = NULL,
                useBytes = TRUE)
    }
    out
  })

  # At the end of a run, whatever this run did not write is left over from an
  # earlier one, a chunk since removed or a figure no longer drawn.
  document_hook <- knitr::knit_hooks$get("document")

  knitr::knit_hooks$set(document = function(x) {
    if (run) {
      save_index()
      kept <- c(basename(index_file), unique(state$written))
      stale <- setdiff(list.files(dir), kept)
      unlink(file.path(dir, stale))
    }
    document_hook(x)
  })

  saved_value <- function(key, expr) {
    if (run) {
      state$saved$values[[key]] <- expr
      save_index()
      return(expr)
    }

    if (is.null(state$saved$values[[key]])) {
      stop(sprintf("No saved value '%s' for vignette '%s'. Render it with BARTISAN_VIGNETTES=run.",
                   key, name), call. = FALSE)
    }

    state$saved$values[[key]]
  }

  assign("saved_value", saved_value, envir = knitr::knit_global())
})
