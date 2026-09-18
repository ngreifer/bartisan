# A zero-augmented gamma with a single additive predictor and a cutpoint: the
# direct [0, Inf) analogue of ordbeta(). One eta drives both the probability of
# a positive value and its conditional mean.
#
#   P(y = 0 | x)     = 1 - expit(eta - cut)
#   y | y > 0, x     ~ Gamma(shape = k, mean = exp(eta))
#
# Prototyped through custom_family() so the model can be checked against a known
# truth before any of it is written in C++.
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

suppressMessages(library(bartisan))

OUT <- "_dev/zero-augmented-prototype.rds"

logdens <- function(y, eta, aux) {
  e <- eta[, 1L]
  cut <- aux[[1L]]
  k <- exp(aux[[2L]])
  lp <- e - cut
  out <- numeric(length(y))
  z <- y > 0
  out[!z] <- -log1p(exp(lp[!z]))                       # log(1 - expit(lp))
  out[z] <- lp[z] - log1p(exp(lp[z])) +
    stats::dgamma(y[z], shape = k, rate = k / exp(e[z]), log = TRUE)
  out
}

set.seed(20260906)
n <- 1500L
p <- 5L
x <- matrix(stats::runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
truth <- 1.5 + 1.2 * sin(pi * x[, 1L]) + 1.0 * x[, 2L] - 0.8 * x[, 3L]
CUT <- 1.0
K <- 2.0
pos <- stats::rbinom(n, 1L, stats::plogis(truth - CUT))
y <- pos * stats::rgamma(n, shape = K, rate = K / exp(truth))
d <- data.frame(y = y, x)
cat(sprintf("simulated: %.0f%% zeros; positive part mean %.1f, sd %.1f\n",
            100 * mean(y == 0), mean(y[y > 0]), stats::sd(y[y > 0])))

fam <- custom_family(logdens, num_predictors = 1L, start = mean(log(y[y > 0])),
                     aux_names = c("cut", "log_shape"),
                     aux_start = c(0, 0), name = "zagamma")

pr <- prog_init(total = 1L, title = "A zero-augmented family in R", unit = "fit",
                kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the fit finished"), add = TRUE)

t0 <- Sys.time()
set.seed(1)
fit <- bartisan(y ~ ., data = d, family = fam, chains = 2,
                num_burn = 400, num_draws = 400, verbose = FALSE)
secs <- as.numeric(Sys.time() - t0, units = "secs")
prog_tick(pr, secs = secs, label = "zero-augmented gamma")

eta_hat <- colMeans(fit[["eta"]][[1L]])
aux <- colMeans(fit[["aux"]])
cat(sprintf("\n%.0fs for 2 chains x 800 sweeps (logdens is evaluated in R)\n", secs))
cat(sprintf("eta:        cor with truth %.3f, RMSE %.3f, bias %+.3f\n",
            stats::cor(eta_hat, truth), sqrt(mean((eta_hat - truth)^2)),
            mean(eta_hat - truth)))
cat(sprintf("cut:        %.2f   (truth %.2f)\n", aux[["cut"]], CUT))
cat(sprintf("shape:      %.2f   (truth %.2f)\n", exp(aux[["log_shape"]]), K))

# The estimand a causal analysis wants, E[y | x], which is closed form here.
mu_hat <- stats::plogis(eta_hat - aux[["cut"]]) * exp(eta_hat)
mu_true <- stats::plogis(truth - CUT) * exp(truth)
cat(sprintf("E[y|x]:     cor with truth %.3f, RMSE %.1f (sd of truth %.1f)\n",
            stats::cor(mu_hat, mu_true), sqrt(mean((mu_hat - mu_true)^2)),
            stats::sd(mu_true)))
saveRDS(list(secs = secs, eta_hat = eta_hat, truth = truth, aux = aux,
             cut_truth = CUT, shape_truth = K, mu_hat = mu_hat,
             mu_true = mu_true, complete = TRUE),
        OUT)

on.exit()
prog_end(pr, "done")
cat("wrote", OUT, "\n")
