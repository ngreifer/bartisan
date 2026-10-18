# The saving half of `vignettes/saved-results.R`, which runs a vignette's chunks
# and writes their output to `vignettes/results/<vignette>/` for the replay that
# CRAN sees.
#
# This is kept out of the package because CRAN does not allow a vignette to
# write into the package directory, and a vignette rendered in place with this
# file does exactly that. It is not meant to be used directly:
# `_dev/refresh-vignettes.R` renders each vignette in `vignettes/` with the
# option `bartisan.save_results` set to a function that sources this file into
# the environment it is given, and `saved-results.R` calls that function with
# its own environment, where `name`, `dir`, `index_file` and `chunk_source()`
# are already defined.
#
# For each chunk with `saved = TRUE`, the finished markdown -- the echoed code,
# the printed output and the links to its figures -- is written to
# `<label>.md`, its figures are copied beside it, and the chunk's header and
# code are recorded in `saved.rds` along with the inline values from
# `saved_value()`.

# A run starts from nothing, so that whatever a chunk no longer produces is not
# carried over.
state <- new.env()
state$saved <- list(chunks = list(), values = list())
state$written <- character()

dir.create(dir, recursive = TRUE, showWarnings = FALSE)

save_index <- function() {
  saveRDS(state$saved, index_file)
}

# A chunk's figures are written under `fig.path`, which rmarkdown deletes
# after rendering. They are copied into the results directory and the markdown
# is pointed at the copies, in a run and so in its replay alike.
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

# Every chunk with the option set is recorded with its source and an empty
# snippet, which the chunk hook overwrites when there is output; it is called
# after this.
knitr::knit_hooks$set(saved = function(before, options, envir) {
  if (!isTRUE(options$saved) || before) {
    return(NULL)
  }

  label <- options$label
  file <- file.path(dir, paste0(label, ".md"))
  file.create(file)
  state$written <- c(state$written, basename(file))
  state$saved$chunks[[label]] <- chunk_source(options$code, options)
  save_index()

  NULL
})

chunk_hook <- knitr::knit_hooks$get("chunk")

knitr::knit_hooks$set(chunk = function(x, options) {
  if (!isTRUE(options$saved)) {
    return(chunk_hook(x, options))
  }

  # knitr passes an excluded chunk through here too and drops it from the page
  # afterwards, so its snippet stays empty.
  out <- keep_figures(chunk_hook(x, options), options)
  if (isTRUE(options$include)) {
    writeChar(out, file.path(dir, paste0(options$label, ".md")), eos = NULL,
              useBytes = TRUE)
  }
  out
})

# At the end, whatever this run did not write is left over from an earlier
# one, a chunk since removed or a figure no longer drawn.
document_hook <- knitr::knit_hooks$get("document")

knitr::knit_hooks$set(document = function(x) {
  save_index()
  kept <- c(basename(index_file), unique(state$written))
  stale <- setdiff(list.files(dir), kept)
  unlink(file.path(dir, stale))
  document_hook(x)
})

assign("saved_value", function(key, expr) {
  state$saved$values[[key]] <- expr
  save_index()
  expr
}, envir = knitr::knit_global())
