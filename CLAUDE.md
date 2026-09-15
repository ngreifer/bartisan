# bartisan

## Standards live in the dotfiles, not here

Before writing anything in this package, read:

- `~/.config/agents/R-PACKAGE-DEV.md` — conventions for code under `R/`, `tests/`,
  `man/`, `DESCRIPTION`, `NAMESPACE`, NEWS and vignettes. **Read this first** for
  package work; it takes precedence over `CODING.md` § R where the two differ.
- `~/.config/agents/CODING.md` — general R style: package loading, pipes,
  function preference, assignment, parallelization.
- `~/.config/agents/WRITING.md` — prose voice for anything a human reads, which
  here means the vignettes, NEWS and roxygen blocks.

These are the authority. This file only records what is specific to this package.

## The check that keeps being missed

`R-PACKAGE-DEV.md` § Markdown files: **a `.md` or `.Rmd` file contains markdown
only, never Rd markup.** Nothing renders `\pkg{}`, `\eqn{}` or `\code{}` in a
markdown file. Pandoc's markdown reader silently deletes the macro *and its
contents*, so the words disappear from the sentence and `R CMD check` stays
silent. It has gone wrong here more than once, most recently across all nine
vignettes and NEWS.md at once.

Run the check from that section before finishing any documentation change:

```sh
grep -n '\\[A-Za-z]\+{' NEWS.md README.md vignettes/*.Rmd
```

Every hit is either a violation or LaTeX inside `$...$`; there is no third case.
Roxygen blocks in `R/` are the only place Rd markup belongs.

## Knit the vignettes against an installed package, not under `load_all()`

`kfold()` sends its K refits to workers with `future.packages = "bartisan"`, and
the thing it sends is the `bartisan` function object itself, since
`kfold_call()` rebuilds the original call. Serializing a closure whose
environment is a namespace writes a *reference* to that namespace rather than
its contents, so the worker resolves it by loading the package by name, which is
the **installed** one. Under `pkgload::load_all()` with an install that predates
the change in hand, any internal function added since is missing there, and
`vignette("comparison")` comes back with

```
#> Error in `prior_record()`: could not find function "prior_record"
```

on every `kfold()` chunk. Nothing is wrong with the code; `R CMD INSTALL .`
first and the same knit is clean. This looks exactly like the real
`diagnosis_block` bug (an unexported function a worker could not find), which
was *not* an install artifact, so the two have to be told apart by installing
rather than by inspection.

## Package specifics

- `_dev/` is development scratch and is gitignored apart from an allowlist in
  `.gitignore`: `TASKS.md`, `SHIP.md`, `benchmark.Rmd`, `check.sh` and the
  `survival-*.R` scripts. `TASKS.md` is the running record of what was tried,
  measured and rejected; `SHIP.md` is the release assessment.
- `src/Makevars` sets `CXX_STD = CXX17`. The code needs it, and without it a
  `-O0` build fails to link where `-O2` silently succeeds.
- The simulations behind `vignette("survival")` are reproducible from
  `_dev/survival-sim.R`, `-bins.R`, `-timing.R` and `-results.R`. The vignette
  reads their saved output so it builds without refitting.
- Timing claims want a quiet machine. `_dev/survival-timing.R` exists because
  timings taken during other work are not comparable.
