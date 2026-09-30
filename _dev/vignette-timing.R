# Vignette build time, against CRAN's check budget.
#
# What is being tested. CRAN expects `R CMD check` to finish in about ten
# minutes, and the check knits every vignette again. Timed on 2026-09-23, seven
# of the eleven vignettes took at least 12 minutes between them (bartisan 238 s,
# causal 139, comparison 127, effects 65-103, families 64, varying 48,
# importance 22), with implementation and diagnostics never timed. The
# non-diagnostics vignettes were then moved from four chains to one. The claim at
# stake: with that change, all eleven vignettes knit in under about 8 minutes in
# total, leaving room for the tests and examples inside the ten.
#
# The measurement. Wall-clock seconds to render each vignette with
# rmarkdown::render(), each in a fresh R process against the installed package
# (never under load_all(); see CLAUDE.md), with `_R_CHECK_LIMIT_CORES_=TRUE` so
# that anything parallel gets the two cores CRAN allows. The column to read is
# `seconds`, and its total; `ok` says whether the knit finished at all, since a
# fast failure is not a fast vignette.
#
# What each outcome means. Under about 8 minutes in total: the vignettes can be
# evaluated in full on CRAN and nothing needs precomputing. Over it: the
# vignettes at the top of the per-vignette table need their heavy chunks
# precomputed or skipped on CRAN, the way vignette("survival") reads saved
# results, and the table says which ones. diagnostics is the expected culprit,
# since its `update()` runs four chains of 2000 + 8000 iterations.
#
# Usage: Rscript _dev/vignette-timing.R <output dir>

args <- commandArgs(trailingOnly = TRUE)
out_dir <- normalizePath(args[1L], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")

pkg <- normalizePath(".")
work <- file.path(out_dir, "vignettes")
unlink(work, recursive = TRUE)
dir.create(work)
file.copy(list.files(file.path(pkg, "vignettes"), full.names = TRUE), work,
          recursive = TRUE)

rmds <- sort(list.files(work, pattern = "[.]Rmd$"))
res_file <- file.path(out_dir, "timing.rds")

pr <- prog_init(total = length(rmds), title = "Vignette build times",
                unit = "vignette", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()

for (f in rmds) {
  log <- file.path(out_dir, sub("[.]Rmd$", ".log", f))
  expr <- sprintf(
    "rmarkdown::render('%s', output_options = list(keep_md = TRUE), quiet = TRUE)",
    file.path(work, f))

  t0 <- Sys.time()
  status <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(expr)),
                    stdout = log, stderr = log,
                    env = "_R_CHECK_LIMIT_CORES_=TRUE")
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  rows[[length(rows) + 1L]] <- data.frame(vignette = f, seconds = round(secs, 1),
                                          ok = identical(status, 0L))
  prog_tick(pr, label = sprintf("%s (%.0f s)", f, secs),
            ok = identical(status, 0L))

  res <- do.call(rbind, rows)
  saveRDS(list(res = res, complete = FALSE, done = nrow(res),
               total = length(rmds)), res_file)
}

res <- do.call(rbind, rows)
res <- res[order(-res$seconds), ]
saveRDS(list(res = res, complete = TRUE, done = nrow(res),
             total = length(rmds)), res_file)

on.exit()
prog_end(pr, "done")

print(res, row.names = FALSE)
cat(sprintf("\nTotal: %.1f minutes\n", sum(res$seconds) / 60))
