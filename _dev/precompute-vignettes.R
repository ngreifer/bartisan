# Render the heavy vignettes once, here, so that CRAN does not render them.
#
# The problem. The eleven vignettes take about 7.5 minutes to build on a quiet
# machine, and six of them hold nearly all of it: diagnostics 217 seconds,
# effects 102, comparison 102, varying 77, bartisan 63 and causal 61, against
# 35 for the other five together (`_dev/vignette-chunk-times.md`). With
# compilation, tests and examples on top, a check is well past the ten minutes
# CRAN budgets for one platform, and the fits are what costs the time rather
# than anything a reader would want shortened.
#
# The shape of the fix. A vignette that must be cheap is written as
# `<name>.Rmd.orig`, the live source, and this script knits it to
# `<name>.Rmd`, the shipped file. Knitting runs every chunk once, here, and
# writes a markdown document in which each chunk has become a fenced block of
# the same code followed by a fenced block of its real output, with figures
# written to `vignettes/figures/` and referenced from the text. The shipped
# `.Rmd` therefore evaluates nothing: rendering it is pandoc's work alone,
# about a second. A reader of the built vignette sees the same code, the same
# printed output and the same figures they would have seen had it run, because
# they are the ones it produced. What they do not see is a note that it was
# precomputed, which is why this script exists rather than a hand-maintained
# pair of files: the two can only agree if one is generated from the other.
#
# Storing the fitted objects instead was the obvious alternative and is not
# available: a four-chain `rhc` fit serializes to 42 megabytes, against the
# five a package may occupy.
#
# What this does not cover. The vignettes left live are checked by CRAN in the
# ordinary way, and the precomputed ones are checked here, by this script, on
# whatever schedule it is run. `_dev/justfile` has a recipe and
# `.github/workflows/precompute-vignettes.yaml` runs it on a push that touches
# a `.orig` file and on request. Nothing stops a `.Rmd` from drifting from its
# `.orig` in between, so the workflow fails if re-knitting changes a file that
# the push did not update.
#
# Run with: Rscript _dev/precompute-vignettes.R [name ...]
#   with no arguments, every `vignettes/*.Rmd.orig`.

A <- path.expand("~/.claude/skills/live-progress/assets")
if (file.exists(file.path(A, "progress.R"))) source(file.path(A, "progress.R"))

stopifnot(requireNamespace("knitr", quietly = TRUE))

root <- getwd()
vdir <- file.path(root, "vignettes")
stopifnot(dir.exists(vdir))

args <- commandArgs(trailingOnly = TRUE)
origs <- {
  if (length(args)) file.path(vdir, paste0(args, ".Rmd.orig"))
  else sort(Sys.glob(file.path(vdir, "*.Rmd.orig")))
}
missing <- origs[!file.exists(origs)]
if (length(missing)) stop("no such source: ", paste(basename(missing), collapse = ", "))
if (!length(origs)) {
  cat("no .Rmd.orig files; nothing to precompute\n")
  quit(save = "no")
}

old <- setwd(vdir)
on.exit(setwd(old), add = TRUE)

have_progress <- exists("prog_init")
pr <- if (have_progress) {
  prog_init(total = length(origs), title = "Precomputing vignettes",
            unit = "vignette", kind = "build")
}

for (orig in origs) {
  name <- sub("\\.Rmd\\.orig$", "", basename(orig))
  t0 <- Sys.time()

  # A fresh environment per vignette, so that one cannot see another's
  # objects and appear to work when it would not on its own.
  env <- new.env(parent = globalenv())

  # The vignette is knitted by `rmarkdown::render()` exactly as a live build
  # would knit it, with the output format its own YAML names, and stopped
  # before pandoc. What it returns is the intermediate markdown of that live
  # build, and that file, with its figures moved, is what ships. An earlier
  # version called `knitr::knit()` directly, which uses knitr's own figure
  # hooks rather than rmarkdown's, and a diff of the two rendered pages found
  # every figure carrying a visible caption, "plot of chunk <label>", that the
  # live page does not have. Taking the live build's own intermediate is what
  # makes the two pages the same rather than similar.
  src <- paste0(name, "__src.Rmd")
  file.copy(basename(orig), src, overwrite = TRUE)
  md <- rmarkdown::render(src, run_pandoc = FALSE, clean = FALSE,
                          envir = env, quiet = TRUE)

  # Figures go to a tracked directory, prefixed by the vignette's name so that
  # two vignettes cannot collide on a chunk label. Stale figures of this
  # vignette go first, so a chunk that was removed does not leave one behind.
  # `html_vignette` embeds them when the shipped file is rendered, so they are
  # read at build time and not installed.
  dir.create("figures", showWarnings = FALSE)
  unlink(Sys.glob(file.path("figures", paste0(name, "-*"))))
  figsrc <- file.path(paste0(name, "__src_files"), "figure-html")
  for (png in list.files(figsrc, full.names = TRUE)) {
    file.copy(png, file.path("figures", paste0(name, "-", basename(png))),
              overwrite = TRUE)
  }

  lines <- readLines(md, warn = FALSE)
  lines <- gsub(paste0(name, "__src_files/figure-html/"),
                paste0("figures/", name, "-"), lines, fixed = TRUE)
  writeLines(lines, paste0(name, ".Rmd"))


  # Knitting evaluates the inline code in the YAML too, which would freeze the
  # date to the day this ran and make a precomputed vignette visibly older than
  # a live one on the same site. The line is put back as it was written, so
  # every vignette still reports the day it was built.
  orig_lines <- readLines(basename(orig), warn = FALSE)
  out_file <- paste0(name, ".Rmd")
  out_lines <- readLines(out_file, warn = FALSE)
  at_orig <- grep("^date:", orig_lines)[1L]
  at_out <- grep("^date:", out_lines)[1L]

  if (!is.na(at_orig) && !is.na(at_out) &&
      grepl("`r ", orig_lines[at_orig], fixed = TRUE)) {
    out_lines[at_out] <- orig_lines[at_orig]
    writeLines(out_lines, out_file)
  }

  # The check that makes "indistinguishable" a property rather than a hope. A
  # live page is pandoc applied to the knitted intermediate, which is `md`;
  # the shipped page is a render of the shipped file, which has no chunks left
  # to knit and so goes straight to pandoc too. If the two pages are not the
  # same bytes, the shipped file is not what a reader would have seen, and
  # that is a failure here rather than something found later on the site.
  page <- function(input) {
    f <- tempfile(fileext = ".html")
    rmarkdown::render(input, output_file = f, envir = new.env(), quiet = TRUE)
    readLines(f, warn = FALSE)
  }
  same <- identical(page(md), page(out_file))

  unlink(c(src, md, sub("\\.knit\\.md$", ".utf8.md", md),
           paste0(name, "__src_files")), recursive = TRUE)

  if (!same) {
    stop(sprintf("%s.Rmd does not reproduce the live render of %s.Rmd.orig",
                 name, name), call. = FALSE)
  }

  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("%-16s %6.1f s  matches the live render\n", name, secs))
  if (have_progress) prog_tick(pr, label = name, secs = secs)
}

if (have_progress) prog_end(pr, "done")

figs <- Sys.glob(file.path(vdir, "figures", "*"))
cat(sprintf("\n%d figure%s, %.1f MB in vignettes/figures/\n", length(figs),
            if (length(figs) == 1L) "" else "s",
            sum(file.size(figs)) / 1024^2))
