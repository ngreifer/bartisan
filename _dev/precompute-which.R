# Which vignettes need precomputing?
#
# What is tested. `_dev/precompute-vignettes.R` renders a vignette once, here,
# so that CRAN does not, and it was first applied to six: diagnostics,
# comparison, effects, varying, bartisan and causal. The convention should
# cover only the vignettes that need it, since a precomputed vignette is
# checked when it is regenerated rather than by CRAN on every submission. The
# two measurements on record disagree about the middle of the ranking. On a
# quiet machine on 2026-09-29: diagnostics 128 s, effects 80, comparison 72,
# varying 53, bartisan 50, causal 46. During the precompute of 2026-10-02:
# diagnostics 230, comparison 142, bartisan 118, varying 112, effects 109,
# causal 107. Only diagnostics and comparison lead in both.
#
# What is measured. Each of the six sources rendered live with
# `rmarkdown::render()`, each in a fresh R process against the installed
# package, one after another on an otherwise idle machine, with
# `_R_CHECK_LIMIT_CORES_=TRUE` as CRAN sets it, two passes. The source is the
# `.Rmd.orig`, copied to a temporary `.Rmd` so that rendering it neither reads
# nor writes the shipped file. Per vignette, wall-clock seconds in each pass and
# their mean.
#
# What each outcome would mean. A vignette whose mean is over about a minute
# stays precomputed and one under it goes back to plain `.Rmd`; the total of
# the ones left live, against the ten minutes CRAN budgets for a whole check
# with compilation, tests and examples on top, says whether that line is in the
# right place or should move.
#
# What happened (2026-10-02). The machine was not idle: two single-threaded
# `_dev/sbc.R` jobs held two of the M4's four performance cores, and every
# vignette came out about 2.3 times its 2026-09-29 figure. With those jobs
# paused, `WHICH_VIGNETTES="varying comparison"` gave 46.8 and 63.4 seconds
# against 53 and 72, so the slowdown was load and the 2026-09-29 numbers stand.
# Diagnostics, effects and comparison stay precomputed; the other three are
# live again. `_dev/TASKS.md` has the entry.
#
# Run with: Rscript _dev/precompute-which.R
# Writes:   _dev/precompute-which.rds, or _dev/precompute-which-<WHICH_TAG>.rds

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

# `WHICH_VIGNETTES` (space separated) and `WHICH_PASSES` re-time a subset, and
# `WHICH_TAG` keeps that run's results apart from the full one's.
VIGNETTES <- strsplit(Sys.getenv("WHICH_VIGNETTES",
  "diagnostics comparison effects varying bartisan causal"), " +")[[1L]]
PASSES <- as.integer(Sys.getenv("WHICH_PASSES", "2"))
TAG <- Sys.getenv("WHICH_TAG")
OUT <- if (nzchar(TAG)) sprintf("_dev/precompute-which-%s.rds", TAG) else "_dev/precompute-which.rds"

vdir <- normalizePath("vignettes")
grid <- expand.grid(vignette = VIGNETTES, pass = seq_len(PASSES),
                    stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "Which vignettes need precomputing",
                unit = "render", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  work <- file.path(tempdir(), sprintf("time-%s-%d", g$vignette, g$pass))
  dir.create(work, showWarnings = FALSE, recursive = TRUE)

  # Everything the vignette reads from its own directory, and the source as a
  # plain `.Rmd`, in a scratch directory of its own.
  for (f in c("references.bib", list.files(vdir, pattern = "\\.rds$"))) {
    if (file.exists(file.path(vdir, f))) file.copy(file.path(vdir, f), work)
  }
  file.copy(file.path(vdir, paste0(g$vignette, ".Rmd.orig")),
            file.path(work, paste0(g$vignette, ".Rmd")))

  script <- sprintf(
    "setwd(%s); t0 <- Sys.time(); rmarkdown::render(%s, quiet = TRUE); cat(as.numeric(difftime(Sys.time(), t0, units = 'secs')))",
    deparse(work), deparse(paste0(g$vignette, ".Rmd")))

  out <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(script)),
                 stdout = TRUE, stderr = FALSE,
                 env = "_R_CHECK_LIMIT_CORES_=TRUE")
  secs <- suppressWarnings(as.numeric(utils::tail(out, 1L)))

  rows[[k]] <- data.frame(vignette = g$vignette, pass = g$pass, seconds = secs,
                          ok = is.finite(secs))
  prog_tick(pr, label = sprintf("%s, pass %d", g$vignette, g$pass),
            ok = is.finite(secs))
  saveRDS(list(rows = do.call(rbind, rows), complete = k == nrow(grid)),
          OUT)
}

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, rows)
out <- do.call(rbind, lapply(split(res, res$vignette), function(x) {
  data.frame(vignette = x$vignette[1L],
             passes = paste(round(x$seconds, 1), collapse = ", "),
             mean = round(mean(x$seconds), 1))
}))
print(out[order(-out$mean), ], row.names = FALSE)
