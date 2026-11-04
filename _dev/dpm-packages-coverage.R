# Do other R implementations of Dirichlet process mixtures of normals cover
# where bartisan's does not, and does drawing the baseline make the difference?
#
# What is tested. bartisan's Dirichlet process families fix the baseline's
# variance scale lambda and mean-scale ratio k_0 at values calibrated from a
# linear fit (George et al. 2019), and their intervals under-cover on bimodal
# and skewed errors (`_dev/survival-skew-coverage.R`, `_dev/dpm-k0-*.R`).
# Other packages differ in exactly this: *dirichletprocess* fixes its baseline
# at constants; *BNPmix*'s `PYdensity()` draws k0 and the variance scale b0 by
# default (`hyper = TRUE`) and fixes them at data-based values otherwise; and
# *spBayesSurv*'s `anovaDDP()`, a survival model on log T with a linear mean,
# fixes its components' variance prior at a scale from a log-normal `survreg()`
# fit. Chattopadhyay, Chakraborty and Dunson (2023, JMLR, sec. "Simulation
# Experiments") report "severe under-coverage" of pointwise density intervals
# with *dirichletprocess* at its fixed baseline. Two claims, either of which
# could be false:
#   (a) For the error density with the error itself observed, the
#       implementations that fix the baseline (*dirichletprocess*; *BNPmix*
#       with `hyper = FALSE`) under-cover on the bimodal and skewed errors, and
#       the one that draws it (*BNPmix* with `hyper = TRUE`) covers clearly
#       better. *bayesm*'s `rDPGibbs()`, which draws every baseline parameter
#       on grids, was meant to be a second; it cannot be, because it returns
#       one draw of the next observation's component per iteration
#       (`thetaNp1draw`, with `probdraw` a column of ones; bayesm 3.1-7,
#       `R/rDPGibbs_rcpp.R` line 163) rather than the mixture, so it gives the
#       predictive density and no pointwise interval for it.
#   (b) For S(t | x) under a linear mean with censoring, `anovaDDP()` with its
#       fixed variance prior under-covers on the bimodal error as `dpm_aft()`
#       does at its defaults, and `dpm_aft()` with k_0 and lambda drawn covers
#       better on the same data.
# They fail if the drawing implementations under-cover as badly as the fixed
# ones (then the shortfall belongs to Dirichlet process mixtures of normals at
# these sample sizes, not to a fixed baseline), or if the fixed ones cover
# (then bartisan's calibration is the outlier, not the practice of fixing).
#
# The measurement. Part A, density: the four error shapes of
# `_dev/dpm-k0-numeric.R` (normal, skewed gamma(1.5), heavy t_3, bimodal
# +-0.9 + N(0, .19), each standardized and times 0.6), n = 500, twenty
# replicates, each package at its documented defaults on the observed errors
# (*dirichletprocess* on standardized data, as its vignette does, 3000 burn-in
# and 2000 kept iterations as in Chattopadhyay et al.; the others 2500 + 2500).
# Per fit: `cover`, the share of 101 points over the middle 96% of the true
# error distribution whose 90% pointwise interval holds the true density
# (READ THIS); `width`; `mae`, the mean absolute error of the posterior mean
# density. The packages see the errors themselves, an advantage over
# bartisan's families, which also estimate the mean; the comparable bartisan
# numbers are `dens_cover` in `_dev/dpm-k0-numeric.R` and `_dev/dpm-k0-oracle.R`.
# Part B, survival: log T = 1 + 0.8 x1 - 0.6 x2 + W with three uniform
# predictors, W log-normal (N(0, .7^2)), bimodal (+-1.3 + N(0, .35^2)) or
# skewed (0.7 times a standardized gamma(1.5)), log censoring times
# N(mean log T + 1.2, 1.5^2) (about 30% censored), n = 700, twenty replicates;
# `anovaDDP()` at its defaults (3000 + 2000), `dpm_aft()` at its defaults and
# with k_0 and lambda drawn (gamma priors of shape 1.5 centered on the fixed
# values), 50 trees, 500 + 500. Per fit, on 100 test subjects by 20 times
# between the 5% and 95% quantiles of log T: `cover` of 90% intervals for
# S(t | x) (READ THIS), `width`, `rmse`, and `cover_inner`, the coverage over
# the cells whose true S(t | x) is below one, added before this run: below the
# skewed error's lower edge the truth is exactly one, which nothing built from
# normals can cover, so on that truth `cover_inner` is the column to read. Methods are compared within data set.
#
# What each outcome would mean.
#   - (a) and (b) hold: drawing the baseline is what separates the
#     implementations that cover from those that do not, and bartisan's default should
#     follow, subject to `_dev/dpm-k0-numeric.R` finding no cost to the mean.
#   - The drawing implementations under-cover too: the shortfall is a property
#     of these mixtures, a default change will not bring nominal coverage, and
#     the documentation has to say what the intervals do and do not cover.
#   - The fixed implementations cover: bartisan's calibration (a linear fit's
#     residual spread, nu = 10) is the problem rather than fixing as such, and
#     the remedy may be a better fixed calibration.
#
# Run with: BARTISAN_LIB=<library with the working tree> \
#           PKG_LIB=<library with BNPmix, dirichletprocess, spBayesSurv> \
#             Rscript _dev/dpm-packages-coverage.R [reps] [workers]
# Writes:   _dev/dpm-packages-coverage.rds, and one file per data set in
#           _dev/dpm-packages-coverage/ as it finishes.

A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

libs <- c(Sys.getenv("BARTISAN_LIB"), Sys.getenv("PKG_LIB"))
.libPaths(c(libs[nzchar(libs)], .libPaths()))

args <- commandArgs(TRUE)
reps <- if (length(args) >= 1L) as.integer(args[1L]) else 20L
workers <- if (length(args) >= 2L) as.integer(args[2L]) else 4L

study <- local({
  units_dir <- "_dev/dpm-packages-coverage"
  level <- 0.90
  lo_q <- (1 - level) / 2
  hi_q <- 1 - lo_q

  # Part A: the shapes of W, each with mean 0 and variance 1.
  shapes <- list(
    normal = list(draw = function(n) stats::rnorm(n),
                  dens = function(w) stats::dnorm(w)),
    skewed = list(
      draw = function(n) (stats::rgamma(n, shape = 1.5) - 1.5) / sqrt(1.5),
      dens = function(w) stats::dgamma(1.5 + sqrt(1.5) * w, 1.5) * sqrt(1.5)),
    heavy = list(draw = function(n) stats::rt(n, df = 3) / sqrt(3),
                 dens = function(w) stats::dt(w * sqrt(3), 3) * sqrt(3)),
    bimodal = list(
      draw = function(n) {
        ifelse(stats::rbinom(n, 1L, 0.5) == 1L, 0.9, -0.9) +
          stats::rnorm(n, sd = sqrt(0.19))
      },
      dens = function(w) {
        0.5 * stats::dnorm(w, 0.9, sqrt(0.19)) +
          0.5 * stats::dnorm(w, -0.9, sqrt(0.19))
      })
  )

  ranges <- lapply(shapes, function(sh) {
    set.seed(1)
    stats::quantile(sh$draw(1e6), c(0.02, 0.98), names = FALSE)
  })

  score_density <- function(D, truth) {
    lo <- apply(D, 2L, stats::quantile, lo_q)
    hi <- apply(D, 2L, stats::quantile, hi_q)
    c(cover = mean(truth >= lo & truth <= hi), width = mean(hi - lo),
      mae = mean(abs(colMeans(D) - truth)))
  }

  density_unit <- function(unit) {
    set.seed(unit$seed)
    n <- 500L
    e <- 0.6 * shapes[[unit$shape]]$draw(n)
    grid <- 0.6 * seq(ranges[[unit$shape]][1L], ranges[[unit$shape]][2L],
                      length.out = 101L)
    truth <- shapes[[unit$shape]]$dens(grid / 0.6) / 0.6

    fits <- list(
      bnpmix_drawn = function() {
        f <- BNPmix::PYdensity(e, mcmc = list(niter = 5000, nburn = 2500,
                                              model = "LS", hyper = TRUE,
                                              print_message = FALSE),
                               prior = list(strength = 1, discount = 0),
                               output = list(grid = grid, out_type = "FULL"))
        f$density
      },
      bnpmix_fixed = function() {
        f <- BNPmix::PYdensity(e, mcmc = list(niter = 5000, nburn = 2500,
                                              model = "LS", hyper = FALSE,
                                              print_message = FALSE),
                               prior = list(strength = 1, discount = 0),
                               output = list(grid = grid, out_type = "FULL"))
        f$density
      },
      dirichletprocess = function() {
        m <- mean(e)
        s <- stats::sd(e)
        dp <- dirichletprocess::Fit(
          dirichletprocess::DirichletProcessGaussian((e - m) / s), 5000,
          progressBar = FALSE)
        t(vapply(3001:5000, function(i) {
          dirichletprocess::PosteriorFunction(dp, i)((grid - m) / s) / s
        }, numeric(length(grid))))
      })

    rows <- lapply(names(fits), function(method) {
      started <- Sys.time()
      D <- fits[[method]]()
      seconds <- as.numeric(Sys.time() - started, units = "secs")
      sc <- score_density(D, truth)
      data.frame(part = "density", truth = unit$shape, rep = unit$rep,
                 method = method, cover = sc[["cover"]], cover_inner = NA_real_,
                 width = sc[["width"]], mae = sc[["mae"]], seconds = seconds)
    })

    do.call(rbind, rows)
  }

  # Part B: survival with a linear mean.
  errors <- list(
    lognormal = list(draw = function(n) stats::rnorm(n, sd = 0.7),
                     surv = function(z) stats::pnorm(z, sd = 0.7,
                                                     lower.tail = FALSE)),
    bimodal = list(
      draw = function(n) {
        ifelse(stats::rbinom(n, 1L, 0.5) == 1L, 1.3, -1.3) +
          stats::rnorm(n, sd = 0.35)
      },
      surv = function(z) {
        0.5 * stats::pnorm(z, 1.3, 0.35, lower.tail = FALSE) +
          0.5 * stats::pnorm(z, -1.3, 0.35, lower.tail = FALSE)
      }),
    skewed = list(
      draw = function(n) 0.7 * (stats::rgamma(n, 1.5) - 1.5) / sqrt(1.5),
      surv = function(z) {
        stats::pgamma(1.5 + sqrt(1.5) * z / 0.7, 1.5, lower.tail = FALSE)
      })
  )

  survival_unit <- function(unit) {
    set.seed(unit$seed)
    n <- 700L
    mk <- function(n) {
      data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                 x3 = stats::runif(n))
    }
    m <- function(d) 1 + 0.8 * d$x1 - 0.6 * d$x2

    train <- mk(n)
    test <- mk(100L)
    log_t <- m(train) + errors[[unit$shape]]$draw(n)
    log_c <- stats::rnorm(n, mean(log_t) + 1.2, 1.5)
    train$time <- exp(pmin(log_t, log_c))
    train$status <- as.numeric(log_t <= log_c)

    times <- exp(stats::quantile(log_t, seq(0.05, 0.95, length.out = 20L),
                                 names = FALSE))
    truth <- t(vapply(log(times), function(l) {
      errors[[unit$shape]]$surv(l - m(test))
    }, numeric(nrow(test))))

    # Cells below the skewed error's lower edge have a true S(t | x) of
    # exactly one, which nothing built from normals can cover.
    inner <- truth < 1 - 1e-9

    score <- function(est, lo, hi) {
      hit <- truth >= lo & truth <= hi
      c(cover = mean(hit), cover_inner = mean(hit[inner]),
        width = mean(hi - lo), rmse = sqrt(mean((est - truth)^2)))
    }

    ctrl <- bartisan::bartisan_control(num_trees = 50, num_burn = 500,
                                       num_draws = 500)
    formula <- survival::Surv(time, status) ~ x1 + x2 + x3

    bart_fit <- function(arm) {
      fam <- bartisan::dpm_aft()
      fam[names(arm)] <- arm
      fit <- bartisan::bartisan(formula, data = train, family = fam,
                                control = ctrl, verbose = FALSE)
      drawn <- predict(fit, newdata = test, type = "survival", times = times,
                       draws = TRUE)
      score(t(apply(drawn, c(2L, 3L), mean)),
            t(apply(drawn, c(2L, 3L), stats::quantile, lo_q)),
            t(apply(drawn, c(2L, 3L), stats::quantile, hi_q)))
    }

    fits <- list(
      anovaDDP = function() {
        out <- utils::capture.output(
          f <- spBayesSurv::anovaDDP(formula, data = train,
                                     mcmc = list(nburn = 3000, nsave = 2000,
                                                 nskip = 0, ndisplay = 1e6)))
        g <- spBayesSurv::GetCurves(f, xnewdata = test, tgrid = times,
                                    CI = level, PLOT = FALSE)
        score(g$Shat, g$Shatlow, g$Shatup)
      },
      # The fixed baseline was the default when this was run; since
      # 2026-11-03 it has to be asked for.
      dpm_aft_fixed = function() bart_fit(list(update_baseline = FALSE)),
      dpm_aft_drawn = function() {
        bart_fit(list(k0_shape = 1.5, lambda_shape = 1.5))
      })

    rows <- lapply(names(fits), function(method) {
      started <- Sys.time()
      sc <- fits[[method]]()
      seconds <- as.numeric(Sys.time() - started, units = "secs")
      data.frame(part = "survival", truth = unit$shape, rep = unit$rep,
                 method = method, cover = sc[["cover"]],
                 cover_inner = sc[["cover_inner"]], width = sc[["width"]],
                 mae = sc[["rmse"]], seconds = seconds)
    })

    do.call(rbind, rows)
  }

  one_unit <- function(unit) {
    suppressMessages(library(survival))
    out <- if (unit$part == "density") density_unit(unit) else survival_unit(unit)
    saveRDS(out, file.path(units_dir, sprintf("%s-%s-%02d.rds", unit$part,
                                              unit$shape, unit$rep)))
    out
  }

  environment()
})

dir.create(study$units_dir, showWarnings = FALSE)

units_a <- expand.grid(rep = seq_len(reps), shape = names(study$shapes),
                       stringsAsFactors = FALSE)
units_a$part <- "density"
units_b <- expand.grid(rep = seq_len(reps), shape = names(study$errors),
                       stringsAsFactors = FALSE)
units_b$part <- "survival"
units <- rbind(units_a, units_b)
units$seed <- 20261104L + seq_len(nrow(units))
units <- split(units, seq_len(nrow(units)))

pr <- prog_init(total = length(units), title = "Other packages' coverage",
                unit = "data set", workers = workers, kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

future::plan(future::multisession, workers = workers)

results <- prog_future_lapply(pr, units, study$one_unit,
                              label = function(u) {
                                paste(u$part, u$shape, u$rep)
                              })

results <- do.call(rbind, Filter(Negate(is.null), results))
saveRDS(list(main = results, complete = TRUE, when = Sys.time()),
        "_dev/dpm-packages-coverage.rds")
on.exit()
prog_end(pr, "done", sprintf("%d data sets", length(units)))
cat("wrote _dev/dpm-packages-coverage.rds\n")
