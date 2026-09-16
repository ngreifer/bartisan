# Two kernels and three bandwidth rules, in the setting that matters: BART's
# held-out accuracy, not the MISE of the distribution function they are derived
# from. Those are different objectives, and the transform exists to serve the
# first.
source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")
suppressMessages(library(bartisan))
SCRATCH <- "/private/tmp/claude-501/-Users-NoahGreifer-Dropbox-Research-R-bartisan/42a89c8f-bd1a-4a13-9717-48692000a78c/scratchpad/xt"
OUT <- file.path(SCRATCH, "kernel-bandwidth.rds")
N_REP <- 10L; N_TEST <- 1000L; SIZES <- c(500L, 2000L)

src <- readLines(file.path(SCRATCH, "xtransform-sim.R"))
eval(parse(text = paste(src[grep("^scenarios <- list", src):
                            (grep("^ctl <- function", src) - 1L)], collapse = "\n")))

scale_of <- function(x) {
  s <- min(stats::sd(x), stats::IQR(x) / 1.349)
  if (!is.finite(s) || s <= 0) s <- stats::sd(x)
  s
}

bz_bandwidth <- function(x, kernel, variant, n_grid = 512L, n_h = 40L) {
  x <- sort(x); n <- length(x); s <- scale_of(x)
  hs <- exp(seq(log(s * n^(-0.7)), log(s * 1.5), length.out = n_h))
  g <- seq(x[1L] - 5 * max(hs), x[n] + 5 * max(hs), length.out = n_grid)
  Fn <- findInterval(g, x) / n
  w <- diff(g)
  obj <- vapply(hs, function(h) {
    p <- bartisan:::.bartisan_smooth_cdf(x, g, h, kernel)
    Fh <- p[, 1L]
    y <- if (identical(variant, "Emp")) {
      p[, 2L] / n^2 - Fh^2 / n + (Fh - Fn)^2
    } else {
      p[, 2L] / (n * (n - 1)) - Fh^2 / (n - 1) +
        (n / (n - 1)) * (Fh - Fn)^2 - p[, 3L] / (n * (n - 1))
    }
    sum(w * (y[-length(y)] + y[-1L]) / 2)
  }, 0)
  hs[which.min(obj)]
}

# Returns a closure mapping new values, fitted on the training column.
fit_map <- function(xtrain, kernel, rule) {
  x <- sort(xtrain[!is.na(xtrain)]); n <- length(x)
  if (length(unique(x)) < 3L) return(function(z) z)
  # A bandwidth is not comparable across kernels: the Epanechnikov reaches
  # exactly h either side where the Gaussian reaches several times that, so the
  # same number smooths less. Calibrated empirically rather than derived: the
  # ratio of the two kernels' Bergmann-Zaehle optima was 2.0 to 2.2 across
  # normal, lognormal, outlier and bimodal samples at two sizes, so the plug-in
  # rule is scaled by 2.1 for the Epanechnikov. The selectors need no such
  # correction, each finding the right bandwidth for its own kernel.
  h <- if (identical(rule, "azzalini")) {
    scale_of(x) * n^(-1 / 3) * (if (kernel == 0L) 2.1 else 1)
  }
  else bz_bandwidth(x, kernel, rule)
  if (!is.finite(h) || h <= 0) return(function(z) z)
  g <- seq(x[1L], x[n], length.out = 1024L)
  at <- bartisan:::.bartisan_smooth_cdf(x, g, h, kernel)[, 1L]
  lo <- at[1L]; hi <- at[length(at)]
  if (!(hi > lo)) return(function(z) z)
  f <- stats::approxfun(g, (at - lo) / (hi - lo), rule = 2L)
  function(z) pmin(pmax(f(z), 0), 1)
}

ARMS <- list(
  gauss_azz = list(kernel = 1L, rule = "azzalini"),
  epan_azz  = list(kernel = 0L, rule = "azzalini"),
  gauss_Emp = list(kernel = 1L, rule = "Emp"),
  epan_Emp  = list(kernel = 0L, rule = "Emp"),
  epan_EmpC = list(kernel = 0L, rule = "EmpC"))

one_cell <- function(scn, n, rep) {
  set.seed(1000L * rep + n)
  tr <- scenarios[[scn]]$gen(n); te <- scenarios[[scn]]$gen(N_TEST)
  sdf <- stats::sd(tr$f); ytr <- tr$f + rnorm(n, sd = 0.3 * sdf)
  out <- NULL
  for (a in names(ARMS)) {
    spec <- ARMS[[a]]
    t0 <- Sys.time()
    maps <- lapply(names(tr$x), function(j) fit_map(tr$x[[j]], spec$kernel, spec$rule))
    names(maps) <- names(tr$x)
    build <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    xtr <- tr$x; xte <- te$x
    for (j in names(maps)) { xtr[[j]] <- maps[[j]](xtr[[j]]); xte[[j]] <- maps[[j]](xte[[j]]) }
    fit <- bartisan(y ~ ., cbind(xtr, y = ytr), family = gaussian(),
                    control = bartisan_control(num_trees = 50L, num_burn = 200L,
                      num_draws = 500L, verbose = FALSE, x_transform = "range"))
    dr <- predict(fit, newdata = xte, type = "response", draws = TRUE)
    fhat <- colMeans(dr)
    lo <- apply(dr, 2L, quantile, .025); hi <- apply(dr, 2L, quantile, .975)
    out <- rbind(out, data.frame(scenario = scn, n = n, rep = rep, arm = a,
      rmse = sqrt(mean((fhat - te$f)^2)) / sdf,
      cover = mean(te$f >= lo & te$f <= hi), build = build))
  }
  out
}

cells <- expand.grid(scn = names(scenarios), n = SIZES, stringsAsFactors = FALSE)
pr <- prog_init(total = nrow(cells) * N_REP, title = "kernel and bandwidth",
                unit = "cell", workers = 1L, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)
rows <- list()
for (k in seq_len(nrow(cells))) {
  for (r in seq_len(N_REP)) {
    rows[[length(rows) + 1L]] <- tryCatch(one_cell(cells$scn[k], cells$n[k], r),
      error = function(e) data.frame(scenario = cells$scn[k], n = cells$n[k],
        rep = r, arm = NA_character_, rmse = NA_real_, cover = NA_real_, build = NA_real_))
    prog_tick(pr, label = sprintf("%s n=%d rep %d", cells$scn[k], cells$n[k], r))
  }
  saveRDS(list(res = do.call(rbind, rows), complete = FALSE), OUT)
}
saveRDS(list(res = do.call(rbind, rows), complete = TRUE), OUT)
on.exit(); prog_end(pr, "done")
