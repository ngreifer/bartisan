# Where does the test suite's time go under CRAN's settings?
#
# What is tested. The full check of 2026-10-02 (`_dev/SHIP.md`) ran
# `tests/testthat.R` in 11.6 minutes with `NOT_CRAN` unset, more than the ten
# CRAN allows for a whole check. The proposal is to run on CRAN only the tests
# whose assertions hold for any seed: argument checks, structure, and exact
# identities such as a replay reproducing the recorded draws. The claim at
# stake is that the time sits in the other kind, the tests whose outcome
# depends on which draws came out (recovery within a tolerance, coverage,
# mixing thresholds), so that skipping those on CRAN brings the suite well
# under budget. It could be false: an identity test that fits a model costs a
# fit whatever it asserts, and the time could sit there instead.
#
# What is measured. Every test file, run with `testthat::test_file()` against
# the package installed from the same commit into a private library, in one R
# process with `NOT_CRAN` unset and `_R_CHECK_LIMIT_CORES_=TRUE` as CRAN sets
# them, on a machine with nothing else running. Per `test_that()` block:
# elapsed seconds, expectations, skips and failures, and per file the elapsed
# seconds including helpers. The blocks are then joined with a classification,
# read from the test code, of each one as seed-free or seed-dependent. The
# number to read is the total time of the seed-free blocks that run on CRAN
# today, against the total of every block that does.
#
# What each outcome would mean. If the seed-free blocks take a few minutes,
# skipping the rest on CRAN is the fix, and it gives up nothing CRAN's run is
# good for, since a tolerance test that fails there for an unlucky seed is a
# false alarm. If they take most of the 11.6 minutes, skipping the
# seed-dependent ones does not help much, and the time has to come from
# shrinking the fits that the identity tests use, since an identity holds at
# any chain length.
#
# Run with: Rscript _dev/test-timing.R <library> <test directory> <output.rds>

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
lib <- normalizePath(args[1L])
tdir <- normalizePath(args[2L])
out <- args[3L]

Sys.unsetenv("NOT_CRAN")
Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")

# The private library goes first, here and in any worker a test starts.
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
stopifnot(dirname(find.package("bartisan")) == lib)

library(testthat)

files <- sort(list.files(tdir, pattern = "^test-.*\\.R$", full.names = TRUE))

pr <- prog_init(total = length(files), title = "Test time under CRAN's settings",
                unit = "file", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

tests <- list()
per_file <- list()
for (k in seq_along(files)) {
  f <- files[k]
  t0 <- Sys.time()
  res <- testthat::test_file(f, reporter = "silent", package = "bartisan",
                             load_package = "installed",
                             stop_on_failure = FALSE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  df <- as.data.frame(res)
  tests[[k]] <- data.frame(file = basename(f), test = df$test, real = df$real,
                           nb = df$nb, skipped = df$skipped,
                           failed = df$failed, error = df$error)
  per_file[[k]] <- data.frame(file = basename(f), seconds = secs,
                              not_cran = Sys.getenv("NOT_CRAN"))

  prog_tick(pr, label = sprintf("%s, %.0f s", basename(f), secs),
            ok = !any(df$failed > 0 | df$error))
  saveRDS(list(tests = do.call(rbind, tests), files = do.call(rbind, per_file),
               complete = k == length(files)), out)
}

on.exit()
prog_end(pr, "done")

pf <- do.call(rbind, per_file)
tt <- do.call(rbind, tests)
cat(sprintf("\n%d files, %.1f minutes; %d blocks, %d skipped, %d with failures or errors\n",
            nrow(pf), sum(pf$seconds) / 60, nrow(tt), sum(tt$skipped),
            sum(tt$failed > 0 | tt$error)))
print(utils::head(pf[order(-pf$seconds), ], 15), row.names = FALSE)
