# Simulation behind vignettes/survival.Rmd.
#
# Six data-generating truths crossed with the five survival families and the
# discrete-time route, scored on held-out data against the true survival
# function. Run with
#
#   Rscript _dev/survival-sim.R
#
# which writes _dev/survival-sim.rds. The vignette embeds the numbers rather
# than sourcing this, so that it builds without refitting anything.

suppressPackageStartupMessages({
  library(bartisan)
  library(survival)
})

N_TRAIN <- 700L
N_TEST  <- 700L
N_REP   <- 5L
N_GRID  <- 20L
CENSOR  <- 0.3

ctrl <- bartisan_control(num_trees = 50, num_burn = 500, num_save = 500)

# ---- covariates and the signal ------------------------------------------------

# One fixed nonlinear function of five covariates, standardized so that its
# spread is the same in every scenario and the scenarios differ only in how the
# time distribution is built around it.
make_x <- function(n) {
  data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n),
             x4 = runif(n), x5 = runif(n))
}

signal <- function(d) {
  raw <- 0.9 * sin(3 * d$x1) + 0.8 * (d$x2 - 0.5)^2 * 4 - 0.7 * d$x3 +
    0.6 * d$x1 * d$x2
  (raw - 0.30) / 0.55
}

XNAMES <- paste0("x", 1:5)
FORM <- Surv(time, status) ~ x1 + x2 + x3 + x4 + x5

# ---- the six truths ----------------------------------------------------------
#
# Each returns the event times and a function giving the true S(t | x) for the
# same rows, so that every family is scored against the same target.

# Proportional hazards, Weibull baseline: the one truth that is both a
# proportional hazards and an accelerated failure time model.
truth_weibull <- function(d, k = 1.6, b = 1) {
  g <- signal(d)
  t_of <- function(u) b * (-log(u) / exp(g))^(1 / k)
  list(time = t_of(runif(nrow(d))),
       surv = function(tt) outer(exp(g), (tt / b)^k, function(e, L) exp(-L * e)))
}

# Proportional hazards with a baseline hazard that rises then falls, which no
# parametric family in the package can represent.
truth_turnover <- function(d, k = 2.5, b = 1) {
  g <- signal(d)
  u <- runif(nrow(d))
  list(time = b * (exp(-log(u) / exp(g)) - 1)^(1 / k),
       surv = function(tt) {
         outer(exp(g), log1p((tt / b)^k), function(e, L) exp(-L * e))
       })
}

# Accelerated failure time with normal errors: lognormal_aft() is exactly right
# and dpm_aft() nests it.
truth_lognormal <- function(d, sigma = 0.7) {
  m <- signal(d)
  list(time = exp(m + sigma * rnorm(nrow(d))),
       surv = function(tt) {
         1 - outer(m, log(tt), function(mm, lt) pnorm((lt - mm) / sigma))
       })
}

# Accelerated failure time with a two-component error, centered so the
# predictor is still the conditional mean of log T.
truth_bimodal <- function(d, mu = 1.3, s = 0.35) {
  m <- signal(d)
  n <- nrow(d)
  comp <- rbinom(n, 1, 0.5)
  w <- ifelse(comp == 1, mu, -mu) + s * rnorm(n)
  list(time = exp(m + w),
       surv = function(tt) {
         1 - outer(m, log(tt), function(mm, lt) {
           0.5 * pnorm((lt - mm - mu) / s) + 0.5 * pnorm((lt - mm + mu) / s)
         })
       })
}

# Accelerated failure time with heavy-tailed errors.
truth_heavy <- function(d, sigma = 0.6, df = 3) {
  m <- signal(d)
  list(time = exp(m + sigma * rt(nrow(d), df)),
       surv = function(tt) {
         1 - outer(m, log(tt), function(mm, lt) pt((lt - mm) / sigma, df))
       })
}

# Non-proportional hazards: the covariate effect reverses sign over time, so
# survival curves cross. Times drawn by inverting a numerically integrated
# cumulative hazard.
truth_crossing <- function(d, tau = 1.2, k = 1.6, b = 1) {
  g <- signal(d)
  n <- nrow(d)
  grid <- seq(0, 6, length.out = 2001)
  # baseline hazard k t^(k-1) / b^k, multiplied by exp(g * (1 - 2 min(t/tau, 1)))
  base_h <- k * pmax(grid, 1e-8)^(k - 1) / b^k
  tilt <- 1 - 2 * pmin(grid / tau, 1)
  cumhaz <- function(gi) {
    h <- base_h * exp(gi * tilt)
    c(0, cumsum((h[-1] + h[-length(h)]) / 2 * diff(grid)))
  }
  H <- vapply(g, cumhaz, numeric(length(grid)))     # grid x n
  target <- -log(runif(n))
  time <- vapply(seq_len(n), function(i) {
    j <- findInterval(target[i], H[, i])
    if (j >= length(grid)) grid[length(grid)] else grid[j + 1]
  }, numeric(1))
  list(time = time,
       surv = function(tt) {
         idx <- pmax(findInterval(tt, grid), 1L)
         exp(-t(H[idx, , drop = FALSE]))
       })
}

TRUTHS <- list(
  "Weibull PH"        = truth_weibull,
  "hazard turns over" = truth_turnover,
  "log-normal errors" = truth_lognormal,
  "bimodal errors"    = truth_bimodal,
  "heavy-tailed errors" = truth_heavy,
  "crossing hazards"  = truth_crossing
)

# ---- censoring ---------------------------------------------------------------

# Independent exponential censoring, its rate solved for the target censoring
# fraction so that every scenario is compared at the same amount of censoring.
apply_censoring <- function(time, target) {
  if (target <= 0) {
    return(list(time = time, status = rep(1, length(time))))
  }
  frac <- function(rate) mean(rexp(length(time), rate) < time)
  rate <- tryCatch(
    uniroot(function(lr) frac(exp(lr)) - target,
            c(log(1e-4), log(1e3)), tol = 1e-3)$root,
    error = function(e) log(1))
  cens <- rexp(length(time), exp(rate))
  list(time = pmin(time, cens), status = as.numeric(time <= cens))
}

# ---- scoring ----------------------------------------------------------------

# Held-out log predictive score, on the density of T for every family. The
# accelerated failure time families report the density of log T for an event and
# the survival function for a censored observation; only the former carries a
# Jacobian, so the correction is -sum(delta * log t).
log_score <- function(fit, test) {
  ld <- predict(fit, newdata = test, type = "density", log = TRUE)
  jac <- if (identical(fit$family$family, "ph")) 0 else
    sum(test$status * log(test$time))
  sum(ld) - jac
}

# RMSE of S(t | x) over the evaluation grid, the worst single time, and how
# well the ordering of subjects by survival is recovered -- averaged over the
# grid rather than read at one time, because under crossing hazards the
# ordering at the median time is uninformative by construction.
surv_error <- function(shat, strue) {
  per_t <- sqrt(colMeans((shat - strue)^2))
  rank <- mean(vapply(seq_len(ncol(shat)), function(j) {
    suppressWarnings(cor(shat[, j], strue[, j], method = "spearman"))
  }, numeric(1)), na.rm = TRUE)
  c(rmse = sqrt(mean((shat - strue)^2)), worst = max(per_t), rank = rank)
}

FAMILIES <- list(
  "weibull_aft()"     = function() weibull_aft(),
  "loglogistic_aft()" = function() loglogistic_aft(),
  "lognormal_aft()"   = function() lognormal_aft(),
  "ph()"              = function() ph(),
  "dpm_aft()"         = function() dpm_aft()
)

# ---- the discrete-time route -------------------------------------------------

expand_dt <- function(d, edges, xnames) {
  reached <- pmin(pmax(findInterval(d$time, edges, rightmost.closed = TRUE), 1L),
                  length(edges))
  long <- d[rep(seq_len(nrow(d)), reached), xnames, drop = FALSE]
  long$tbin <- unlist(lapply(reached, seq_len))
  long$ev <- 0
  long$ev[cumsum(reached)] <- d$status
  long
}

fit_dt <- function(train, test, edges, ctrl) {
  long <- expand_dt(train, edges, XNAMES)
  fit <- bartisan(ev ~ ., data = long, family = binomial("probit"),
                  control = ctrl, verbose = FALSE)
  # hazard at every (test row, bin), then S as the running product
  K <- length(edges)
  grid_dat <- test[rep(seq_len(nrow(test)), each = K), XNAMES, drop = FALSE]
  grid_dat$tbin <- rep(seq_len(K), nrow(test))
  h <- matrix(predict(fit, newdata = grid_dat, type = "response"),
              nrow = K, ncol = nrow(test))
  list(surv = t(apply(1 - h, 2, cumprod)), rows = nrow(long))
}

# ---- one replicate ----------------------------------------------------------

one_rep <- function(truth_name, rep, censor = CENSOR, families = FAMILIES) {
  set.seed(20260828L + 1000L * match(truth_name, names(TRUTHS)) + rep)

  gen <- TRUTHS[[truth_name]]
  train <- make_x(N_TRAIN)
  test  <- make_x(N_TEST)
  tr <- gen(train)
  te <- gen(test)

  ctr <- apply_censoring(tr$time, censor)
  train$time <- ctr$time; train$status <- ctr$status
  cte <- apply_censoring(te$time, censor)
  test$time <- cte$time; test$status <- cte$status

  # A common evaluation grid: quantiles of the uncensored event times, which is
  # also the grid the discrete-time model is given.
  edges <- unique(quantile(tr$time, seq(0.05, 0.95, length.out = N_GRID)))
  strue <- te$surv(edges)

  # Three subjects for the survival-curve figures, picked by their true
  # survival at the middle of the grid so that the labels mean the same thing
  # whether the truth is built on the hazard scale or the time scale.
  mid <- strue[, ncol(strue) %/% 2]
  who <- c("poor prognosis" = which.min(abs(mid - quantile(mid, 0.1))),
           "typical"        = which.min(abs(mid - median(mid))),
           "good prognosis" = which.min(abs(mid - quantile(mid, 0.9))))

  out <- list()
  for (nm in names(families)) {
    tick <- proc.time()[["elapsed"]]
    fit <- bartisan(FORM, data = train, family = families[[nm]](),
                    control = ctrl, verbose = FALSE)
    secs <- proc.time()[["elapsed"]] - tick
    shat <- predict(fit, newdata = test, type = "survival", times = edges)
    err <- surv_error(shat, strue)
    out[[nm]] <- data.frame(
      truth = truth_name, family = nm, rep = rep, censor = censor,
      s_rmse = err[["rmse"]], s_worst = err[["worst"]],
      logscore = log_score(fit, test), rank = err[["rank"]], secs = secs)

    # Fitted curves for a few representative subjects, kept from the first
    # replicate only, for the figures.
    if (rep == 1L && censor == CENSOR) {
      curves[[nm]] <<- data.frame(
        truth = truth_name, family = nm,
        subject = rep(names(who), each = length(edges)),
        time = rep(edges, length(who)),
        surv = as.vector(t(shat[who, , drop = FALSE])))
    }
    if (rep == 1L && censor == CENSOR && identical(nm, "dpm_aft()")) {
      ed <- error_density(fit)
      dens[[truth_name]] <<- data.frame(truth = truth_name, at = ed$at,
                                        mean = ed$mean, lower = ed$lower,
                                        upper = ed$upper)
    }
  }

  tick <- proc.time()[["elapsed"]]
  dt <- fit_dt(train, test, edges, ctrl)
  err <- surv_error(dt$surv, strue)
  out[["discrete time"]] <- data.frame(
    truth = truth_name, family = "discrete-time probit", rep = rep,
    censor = censor, s_rmse = err[["rmse"]], s_worst = err[["worst"]],
    logscore = NA_real_, rank = err[["rank"]],
    secs = proc.time()[["elapsed"]] - tick)

  if (rep == 1L && censor == CENSOR) {
    curves[["discrete time"]] <<- data.frame(
      truth = truth_name, family = "discrete-time probit",
      subject = rep(names(who), each = length(edges)),
      time = rep(edges, length(who)),
      surv = as.vector(t(dt$surv[who, , drop = FALSE])))
    curves[["truth"]] <<- data.frame(
      truth = truth_name, family = "truth",
      subject = rep(names(who), each = length(edges)),
      time = rep(edges, length(who)),
      surv = as.vector(t(strue[who, , drop = FALSE])))
  }

  do.call(rbind, out)
}

curves <- list()
dens <- list()
all_curves <- list()
all_dens <- list()

# ---- run --------------------------------------------------------------------

main <- list()
for (nm in names(TRUTHS)) {
  for (r in seq_len(N_REP)) {
    curves <- list(); dens <- list()
    main[[length(main) + 1L]] <- one_rep(nm, r)
    if (r == 1L) {
      all_curves[[nm]] <- do.call(rbind, curves)
      all_dens[[nm]] <- do.call(rbind, dens)
    }
    cat(sprintf("[main] %-20s rep %d done\n", nm, r)); flush(stdout())
  }
}
main <- do.call(rbind, main)
all_curves <- do.call(rbind, all_curves)
all_dens <- do.call(rbind, all_dens)

# Censoring sweep on the turning-over baseline, where the families disagree
# most, to see whether the ordering survives heavy censoring.
sweep <- list()
for (cz in c(0, 0.25, 0.5, 0.7)) {
  for (r in 1:3) {
    sweep[[length(sweep) + 1L]] <- one_rep("hazard turns over", r, censor = cz)
    cat(sprintf("[sweep] censor %.2f rep %d done\n", cz, r)); flush(stdout())
  }
}
sweep <- do.call(rbind, sweep)

saveRDS(list(main = main, sweep = sweep, curves = all_curves,
             dens = all_dens, when = Sys.time()),
        "_dev/survival-sim.rds")
cat("wrote _dev/survival-sim.rds\n")
