# A project-level .Rprofile REPLACES the user-level one -- R sources only the
# first it finds. Source the global one explicitly so the usual devtools/usethis
# setup, options and helpers are still available, then set this project's own
# options below so they win on any conflict.
if (file.exists("~/.Rprofile")) {
  source("~/.Rprofile")
}

# devtools and pkgload compile this package's C++ with `-O0` by default, through
# pkgbuild's "extra flags". For a package whose whole point is a sampler, that is
# a five- to twentyfold slowdown that looks like nothing at all: the fits are
# correct, the tests pass, and only the clock is wrong. It has cost real time
# here more than once.
#
# Turning the injection off leaves R's own flags in place, which is `-O2` from
# `R CMD config CXXFLAGS`. Measured on this package, `-O3`, `-mcpu=native` and
# `-flto` are all within one percent of `-O2`, so there is nothing to gain from
# going further -- the only thing that matters is not being at `-O0`.
#
# `bartisan:::.bartisan_optimized()` reports what the loaded library actually is,
# and `bartisan()` warns once per session when the answer is no.
options(pkg.build_extra_flags = FALSE)

if (interactive()) {
  message("bartisan: compiling at -O2 (pkg.build_extra_flags = FALSE)")
}
