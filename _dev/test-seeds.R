# Do the tests that run on CRAN hold for any seed?
#
# What is tested. Since 2026-10-02 the rule in `tests/testthat/helper-bartisan.R`
# is that a test runs on CRAN only if its assertions hold whatever the draws,
# and the fits of many such tests were shrunk to `quick_control()` on the
# strength of it. The claim could be false in two ways: a block classified as
# seed-free from its code may not be, and an assertion that holds only almost
# surely on a long chain (a forest splits at least once, two fits differ, a
# warning does not fire) can stop holding on a short one.
#
# What is measured. The CRAN subset, with `NOT_CRAN` unset, run once per offset
# with every `set.seed(s)` the tests make replaced by `set.seed(s + offset)`, so
# that each run sees different data and different draws from the committed
# ones. The replacement is a helper file written into a temporary copy of the
# test directory; the package's tests are not changed. Per offset, the blocks
# that fail or error. The list to read is every block that fails under any
# offset.
#
# What each outcome would mean. If none fails, the CRAN subset is seed-free in
# practice at the sizes it now uses, and an unlucky draw on a CRAN machine
# cannot fail it. A block that fails under some offset depends on the draws:
# its assertion or its fit size is fixed, or it moves behind skip_on_cran().
#
# Run with: Rscript _dev/test-seeds.R <library> <test directory> <offsets> <output.rds>
#   offsets as "1:8"

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
lib <- normalizePath(args[1L])
src <- normalizePath(args[2L])
offsets <- eval(str2lang(args[3L]))
out <- args[4L]

Sys.unsetenv("NOT_CRAN")
Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
stopifnot(dirname(find.package("bartisan")) == lib)

tdir <- file.path(tempdir(), "seeded-tests")
dir.create(tdir, showWarnings = FALSE)
file.copy(list.files(src, full.names = TRUE), tdir, overwrite = TRUE)

# Sourced after helper-bartisan.R, into the same environment the test files run
# in, so every set.seed() they and the helpers make finds this one first.
writeLines(c(
  '.seed_offset <- as.integer(Sys.getenv("BARTISAN_SEED_OFFSET", "0"))',
  'if (.seed_offset != 0L) {',
  '  set.seed <- function(seed, ...) base::set.seed(seed + .seed_offset, ...)',
  '}'), file.path(tdir, "helper-zz-seed-offset.R"))

files <- sort(list.files(tdir, pattern = "^test-.*\\.R$", full.names = TRUE))
grid <- expand.grid(file = files, offset = offsets, stringsAsFactors = FALSE)

pr <- prog_init(total = nrow(grid), title = "CRAN tests under shifted seeds",
                unit = "file", kind = "measurement")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  Sys.setenv(BARTISAN_SEED_OFFSET = g$offset)
  res <- as.data.frame(testthat::test_file(g$file, reporter = "silent",
                                           package = "bartisan",
                                           load_package = "installed",
                                           stop_on_failure = FALSE))
  bad <- res[res$failed > 0 | res$error, , drop = FALSE]
  rows[[k]] <- data.frame(offset = rep(g$offset, nrow(bad)),
                          file = rep(basename(g$file), nrow(bad)),
                          test = bad$test)
  prog_tick(pr, label = sprintf("offset %d, %s", g$offset, basename(g$file)),
            ok = nrow(bad) == 0L)
  saveRDS(list(failures = do.call(rbind, rows), done = k, total = nrow(grid)),
          out)
}

on.exit()
prog_end(pr, "done")

fails <- do.call(rbind, rows)
cat(sprintf("\n%d offsets x %d files: %d failing blocks\n", length(offsets),
            length(files), nrow(fails)))
if (nrow(fails)) print(fails, row.names = FALSE)
