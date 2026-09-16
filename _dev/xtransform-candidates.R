# Candidate replacements for `x_transform`, measured on the same eight
# data-generating processes as `xtransform-sim.R`.
#
# `range` failed only where the bulk of x is compressed into a sliver of the
# range AND the truth needs fine structure inside it, because a cutpoint drawn
# uniformly on a node's live range then never lands there. Everything below
# attacks that by compressing the tail instead of the bulk, while staying
# monotone and (unlike the ECDF) differentiable, so a soft-rule fit keeps a
# derivative.
#
# Each is a monotone map to a bounded interval, which is exactly what
# `make_unit_map()` returns, so a prototype that pre-transforms and then fits
# with `x_transform = "range"` is the same model: `range` on an already-monotone
# image is an affine rescale, and a cutpoint uniform on the live range is
# affine-equivariant.
#
#   robust     logistic squash on a median/MAD scale. Smooth everywhere, and
#              outliers saturate instead of stretching the coordinate.
#   smoothcdf  kernel-smoothed ECDF, bandwidth at Azzalini's (1981) n^(-1/3)
#              rate for the distribution function. This is the closest thing
#              reachable without touching the engine to "cutpoints on the
#              quantile scale, coordinate smooth".
#   winsor     cap at the 1st and 99th percentiles, then range. Crude, and the
#              baseline the other two have to beat.
source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages(library(bartisan))

SCRATCH <- "/private/tmp/claude-501/-Users-NoahGreifer-Dropbox-Research-R-bartisan/42a89c8f-bd1a-4a13-9717-48692000a78c/scratchpad/xt"
OUT <- file.path(SCRATCH, "candidates.rds")
N_REP <- 12L; N_TEST <- 1000L; SIZES <- c(500L, 2000L)

src <- readLines(file.path(SCRATCH, "xtransform-sim.R"))
eval(parse(text = paste(src[grep("^scenarios <- list", src):
                            (grep("^ctl <- function", src) - 1L)], collapse = "\n")))

# Fitted on the training column, then applied to both, so nothing leaks.
fit_map <- function(kind, xtrain) {
  if (identical(kind, "robust")) {
    med <- stats::median(xtrain); s <- stats::mad(xtrain)
    if (!is.finite(s) || s <= 0) s <- stats::IQR(xtrain) / 1.349
    if (!is.finite(s) || s <= 0) s <- stats::sd(xtrain)
    if (!is.finite(s) || s <= 0) return(identity)
    return(function(z) stats::plogis((z - med) / s))
  }
  if (identical(kind, "smoothcdf")) {
    n <- length(xtrain)
    s <- min(stats::sd(xtrain), stats::IQR(xtrain) / 1.349)
    if (!is.finite(s) || s <= 0) s <- stats::sd(xtrain)
    if (!is.finite(s) || s <= 0) return(identity)
    h <- s * n^(-1 / 3)
    xs <- sort(xtrain)
    return(function(z) vapply(z, function(v) mean(stats::pnorm((v - xs) / h)), 0))
  }
  if (identical(kind, "winsor")) {
    q <- stats::quantile(xtrain, c(.01, .99), names = FALSE)
    return(function(z) pmin(pmax(z, q[1L]), q[2L]))
  }
  identity
}

apply_maps <- function(x, maps) {
  for (j in names(maps)) x[[j]] <- maps[[j]](x[[j]])
  x
}

ARMS <- c("quantile", "range", "robust", "smoothcdf", "winsor")

one_cell <- function(scn, n, rep) {
  set.seed(1000L * rep + n)
  tr <- scenarios[[scn]]$gen(n)
  te <- scenarios[[scn]]$gen(N_TEST)
  sdf <- stats::sd(tr$f)
  train_y <- tr$f + rnorm(n, sd = 0.3 * sdf)

  out <- NULL
  for (arm in ARMS) {
    native <- if (arm %in% c("quantile", "range")) arm else "range"
    maps <- if (arm %in% c("quantile", "range")) list() else
      setNames(lapply(names(tr$x), function(j) fit_map(arm, tr$x[[j]])), names(tr$x))

    xtr <- apply_maps(tr$x, maps)
    xte <- apply_maps(te$x, maps)

    fit <- bartisan(y ~ ., cbind(xtr, y = train_y), family = gaussian(),
                    control = bartisan_control(num_trees = 50L, num_burn = 200L,
                                               num_draws = 500L, verbose = FALSE,
                                               x_transform = native))
    dr <- predict(fit, newdata = xte, type = "response", draws = TRUE)
    fhat <- colMeans(dr)
    lo <- apply(dr, 2L, quantile, 0.025); hi <- apply(dr, 2L, quantile, 0.975)
    out <- rbind(out, data.frame(
      scenario = scn, n = n, rep = rep, arm = arm,
      rmse = sqrt(mean((fhat - te$f)^2)) / sdf,
      cover = mean(te$f >= lo & te$f <= hi),
      width = mean(hi - lo) / sdf))
  }
  out
}

cells <- expand.grid(scn = names(scenarios), n = SIZES, stringsAsFactors = FALSE)
pr <- prog_init(total = nrow(cells) * N_REP, title = "x_transform candidates",
                unit = "cell", workers = 1L, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- list()
for (k in seq_len(nrow(cells))) {
  for (r in seq_len(N_REP)) {
    res <- tryCatch(one_cell(cells$scn[k], cells$n[k], r),
                    error = function(e) data.frame(scenario = cells$scn[k],
                      n = cells$n[k], rep = r, arm = NA_character_,
                      rmse = NA_real_, cover = NA_real_, width = NA_real_))
    rows[[length(rows) + 1L]] <- res
    prog_tick(pr, label = sprintf("%s n=%d rep %d", cells$scn[k], cells$n[k], r))
  }
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE), OUT)
}
saveRDS(list(res = do.call(rbind, rows), complete = TRUE), OUT)
on.exit(); prog_end(pr, "done")
