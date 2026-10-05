# Rerun the vignettes that keep saved results, save the results again, and check
# that replaying them reproduces each page.
#
# `vignettes/saved-results.R` holds the mechanism. A vignette that sources it
# runs its chunks and saves their output in "run" mode, and prints the saved
# output without running anything in "replay" mode, which is what CRAN sees.
# `R CMD build` and `R CMD check` knit a temporary copy of the package, so the
# saved results can only be refreshed by rendering in place, which is what this
# does: each vignette is rendered in run mode in `vignettes/`, which rewrites
# `vignettes/results/<name>/`, and then again in replay mode, and the two pages
# have to be the same bytes. That comparison is what makes "a reader cannot
# tell the difference" a property of the shipped vignettes rather than a hope.
#
# Each render is its own R process with the mode in its environment, and writes
# its page to a temporary directory, so that nothing but the results is left in
# `vignettes/`. Every package in Suggests has to be installed, because a chunk
# that skips itself for a missing package would save its code with no output,
# and that is what would ship.
#
# Run with: Rscript _dev/refresh-vignettes.R [name ...]
#   with no arguments, every vignette that sources saved-results.R.

A <- path.expand("~/.claude/skills/live-progress/assets")
if (file.exists(file.path(A, "progress.R"))) source(file.path(A, "progress.R"))

vdir <- normalizePath("vignettes")
stopifnot(dir.exists(vdir))

uses_saved <- function(f) {
  any(grepl("saved-results.R", readLines(f, warn = FALSE), fixed = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
sources <- {
  if (length(args)) file.path(vdir, paste0(args, ".Rmd"))
  else Filter(uses_saved, sort(Sys.glob(file.path(vdir, "*.Rmd"))))
}

missing <- sources[!file.exists(sources)]
if (length(missing)) {
  stop("no such vignette: ", paste(basename(missing), collapse = ", "))
}

plain <- sources[!vapply(sources, uses_saved, logical(1L))]
if (length(plain)) {
  stop("these vignettes do not keep saved results: ",
       paste(basename(plain), collapse = ", "))
}

suggests <- read.dcf("DESCRIPTION", fields = "Suggests")[1L, 1L] |>
  strsplit(",") |>
  unlist() |>
  trimws()
suggests <- sub("[ (].*", "", suggests[nzchar(suggests)])
absent <- suggests[!vapply(suggests, requireNamespace, logical(1L),
                           quietly = TRUE)]
if (length(absent)) {
  stop("install these Suggests first, or their chunks would be saved without ",
       "output: ", paste(absent, collapse = ", "))
}

render <- function(src, mode) {
  out <- tempfile(fileext = ".html")
  script <- sprintf("rmarkdown::render(%s, output_file = %s, quiet = TRUE)",
                    deparse(src), deparse(out))
  t0 <- Sys.time()
  log <- suppressWarnings(
    system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(script)),
            env = paste0("BARTISAN_VIGNETTES=", mode), stdout = TRUE,
            stderr = TRUE))

  if (!is.null(attr(log, "status")) || !file.exists(out)) {
    cat(utils::tail(log, 20L), sep = "\n")
    stop(sprintf("%s failed in %s mode", basename(src), mode), call. = FALSE)
  }

  list(page = out,
       secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

have_progress <- exists("prog_init")
pr <- if (have_progress) {
  prog_init(total = 2L * length(sources), title = "Refreshing saved vignette results",
            unit = "render", kind = "build")
}

rows <- list()
for (src in sources) {
  name <- sub("[.]Rmd$", "", basename(src))

  ran <- render(src, "run")
  if (have_progress) prog_tick(pr, label = sprintf("%s, run", name))

  replayed <- render(src, "replay")
  if (have_progress) prog_tick(pr, label = sprintf("%s, replay", name))

  same <- identical(unname(tools::md5sum(ran$page)),
                    unname(tools::md5sum(replayed$page)))
  files <- list.files(file.path(vdir, "results", name), full.names = TRUE)

  rows[[name]] <- data.frame(vignette = name, run = round(ran$secs),
                             replay = round(replayed$secs, 1),
                             files = length(files),
                             kb = round(sum(file.size(files)) / 1024),
                             identical = same)
}

if (have_progress) prog_end(pr, "done")

res <- do.call(rbind, rows)
print(res, row.names = FALSE)

if (!all(res$identical)) {
  stop("a replay did not reproduce its run: ",
       paste(res$vignette[!res$identical], collapse = ", "), call. = FALSE)
}

cat("\nevery replay reproduces its run; commit vignettes/results/\n")
