# Time of every chunk in every vignette, and two chain-count variants.
#
# What is being tested. The eleven vignettes knit in 9.2 minutes
# (_dev/vignette-timing.R, 2026-09-24), which is over CRAN's check budget before
# compilation, tests and examples are added. The claim at stake is that the time
# is concentrated in a handful of chunks (diagnostics' `update()` to 2000 + 8000
# draws and comparison's two `kfold()` calls are the expected ones) rather than
# spread thinly across all of them. Two chain-count questions ride along: whether
# diagnostics at 2 chains instead of 4 saves much, and whether causal at 2 chains
# (with the multisession plan it used to have) costs much more than at 1.
#
# The measurement. Wall-clock seconds per chunk from a knitr chunk hook, each
# vignette rendered by rmarkdown::render() in a fresh R process against the
# installed package, with `_R_CHECK_LIMIT_CORES_=TRUE` so that parallel code
# gets CRAN's two cores. The total render time is recorded as well, so that the
# overhead outside the chunks (startup, pandoc) is the difference. The variants
# are copies with `chains = 4` replaced by `chains = 2` (diagnostics), and with
# `chains = 2` added to the four fits and the plan restored (causal). Read
# `seconds` per chunk, and the variants' totals against the current files'.
#
# What each outcome means. Concentrated: cutting or precomputing a few chunks
# brings the vignettes under budget and nothing else needs to change. Spread
# thinly: only a global approach works (skipping evaluation on CRAN, or
# precomputing whole vignettes). Causal at 2 chains within about 20 seconds of
# its 1-chain time: reverting is cheap. Roughly double: it competes with the
# cuts. Diagnostics at 2 chains saving little: chain count is not the lever
# there, and the long refit or the diagnostics themselves are.
#
# Usage: Rscript _dev/vignette-chunk-timing.R <output dir> [vignette names...]
# With no names, every vignette is timed and neither variant is built.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- normalizePath(args[1L], mustWork = FALSE)
only <- args[-1L]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")

pkg <- normalizePath(".")
work <- file.path(out_dir, "vignettes")
unlink(work, recursive = TRUE)
dir.create(work)
file.copy(list.files(file.path(pkg, "vignettes"), full.names = TRUE), work,
          recursive = TRUE)

# The two variants, written beside the originals so that relative paths
# (references.bib, survival-results.rds) still resolve. `sub` edits every line
# containing the pattern; `line` replaces lines that match it exactly, which is
# what keeps one fit from being edited twice.
variant <- function(from, to, sub = list(), line = list()) {
  s <- readLines(file.path(work, from))
  for (e in sub) {
    hit <- grepl(e[1L], s, fixed = TRUE)
    stopifnot(any(hit))
    s[hit] <- gsub(e[1L], e[2L], s[hit], fixed = TRUE)
  }
  for (e in line) {
    hit <- s == e[1L]
    stopifnot(sum(hit) == 1L)
    s[hit] <- e[2L]
  }
  writeLines(s, file.path(work, to))
}

# The variants are built only when named, since they are edits of the files as
# they stood on 2026-09-24 and stop matching once those files move on.
wanted <- function(stem) stem %in% only

if (wanted("diagnostics-2chains")) variant("diagnostics.Rmd", "diagnostics-2chains.Rmd",
        sub = list(c("chains = 4", "chains = 2")))

if (wanted("causal-2chains")) variant("causal.Rmd", "causal-2chains.Rmd", line = list(
  c("run <- TRUE",
    "run <- TRUE; if (rlang::is_installed(\"future\")) future::plan(future::multisession)"),
  c("  data = rhc, family = binomial()",
    "  data = rhc, family = binomial(), chains = 2"),
  c("                data = rhc, family = binomial(), sparsity = FALSE)",
    "                data = rhc, family = binomial(), sparsity = FALSE, chains = 2)"),
  c("               family = binomial())",
    "               family = binomial(), chains = 2)"),
  c("  data = lalonde, family = dpm()",
    "  data = lalonde, family = dpm(), chains = 2")))

rmds <- sort(list.files(work, pattern = "[.]Rmd$"))
if (length(only)) rmds <- rmds[sub("[.]Rmd$", "", rmds) %in% only]

hook <- '
knitr::knit_hooks$set(timeit = function(before, options, envir) {
  if (before) {
    assign(".chunk_t0", proc.time()[["elapsed"]], envir = globalenv())
  } else {
    dt <- proc.time()[["elapsed"]] - get(".chunk_t0", envir = globalenv())
    code <- trimws(options$code)
    code <- code[nzchar(code) & !startsWith(code, "#")]
    row <- data.frame(label = options$label, engine = options$engine,
                      eval = !isFALSE(options$eval), seconds = round(dt, 2),
                      first_line = substr(if (length(code)) code[1L] else "", 1L, 70L))
    f <- Sys.getenv("CHUNK_CSV")
    utils::write.table(row, f, sep = "\\t", row.names = FALSE, quote = FALSE,
                       append = file.exists(f), col.names = !file.exists(f))
  }
  invisible(NULL)
})
knitr::opts_chunk$set(timeit = TRUE)
'

res_file <- file.path(out_dir, "chunk-timing.rds")

pr <- prog_init(total = length(rmds), title = "Vignette chunk times",
                unit = "vignette", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

totals <- list()
chunks <- list()

for (f in rmds) {
  stem <- sub("[.]Rmd$", "", f)
  log <- file.path(out_dir, paste0(stem, ".log"))
  csv <- file.path(out_dir, paste0(stem, ".tsv"))
  unlink(csv)

  script <- file.path(out_dir, paste0(stem, "-render.R"))
  writeLines(c(hook,
               sprintf("rmarkdown::render('%s', output_options = list(keep_md = TRUE), quiet = TRUE)",
                       file.path(work, f))),
             script)

  t0 <- Sys.time()
  status <- system2(file.path(R.home("bin"), "Rscript"), script,
                    stdout = log, stderr = log,
                    env = c("_R_CHECK_LIMIT_CORES_=TRUE",
                            paste0("CHUNK_CSV=", csv)))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- identical(status, 0L)

  ch <- if (file.exists(csv)) {
    read.delim(csv, quote = "", comment.char = "", stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  if (nrow(ch)) {
    ch <- cbind(vignette = stem, order = seq_len(nrow(ch)), ch)
    chunks[[stem]] <- ch
  }

  totals[[stem]] <- data.frame(vignette = stem, total = round(secs, 1),
                               in_chunks = round(sum(ch$seconds), 1),
                               n_chunks = nrow(ch), ok = ok)
  prog_tick(pr, label = sprintf("%s (%.0f s)", f, secs), ok = ok)

  saveRDS(list(totals = do.call(rbind, totals), chunks = do.call(rbind, chunks),
               complete = FALSE, done = length(totals), total = length(rmds)),
          res_file)
}

totals <- do.call(rbind, totals)
chunks <- do.call(rbind, chunks)
saveRDS(list(totals = totals, chunks = chunks, complete = TRUE,
             done = nrow(totals), total = length(rmds)), res_file)

on.exit()
prog_end(pr, "done")

print(totals[order(-totals$total), ], row.names = FALSE)
