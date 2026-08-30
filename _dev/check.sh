#!/bin/sh
#
# Build and check the package without leaving anything in the package
# directory. `R CMD build` always writes its tarball to the working directory
# and has no flag for it, so the build is run *from* the scratch directory with
# the package as its argument; `R CMD check` does have `--output`.
#
# Why bother: the package lives in a Dropbox-synced tree, and a `.Rcheck`
# directory is a few hundred regenerable files that would sync, plus a tarball
# that would sit in the repo root asking to be gitignored.
#
# Usage:  _dev/check.sh [extra R CMD check flags]
# Scratch: $GENBART_CHECK_DIR, or $TMPDIR/bartisan-check, or /tmp/bartisan-check.

set -e

pkg=$(cd "$(dirname "$0")/.." && pwd)
scratch=${GENBART_CHECK_DIR:-${TMPDIR:-/tmp}/bartisan-check}

mkdir -p "$scratch"

printf 'documenting\n'
Rscript -e 'suppressMessages(roxygen2::roxygenise("."))' 2>&1 |
  grep -v 'Operation not permitted' || true

cd "$scratch"

printf 'building in %s\n' "$scratch"
R CMD build "$pkg" > build.log 2>&1 || { tail -30 build.log; exit 1; }

tarball=$(ls -t "$scratch"/*.tar.gz | head -1)

printf 'checking %s\n' "$(basename "$tarball")"
R CMD check --no-manual --output="$scratch" "$@" "$tarball" > check.log 2>&1 || true

# The OpenMP runtime writes a line to stderr about /tmp in this sandbox, and
# R CMD check reports any stderr during a step as a problem. It is not one.
grep -E '^\* checking|^Status:' check.log |
  grep -vE '\.\.\. OK$' || true

printf '\nlogs: %s/check.log, %s/build.log\n' "$scratch" "$scratch"
