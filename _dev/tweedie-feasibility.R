# Tweedie ties the probability of a zero to the mean deterministically:
#   P(y = 0 | x) = exp(-mu(x)^(2-p) / (phi (2-p)))
# so there is no freedom in the level once (phi, p) are set. The cutpoint model
# has a free level. Which of those the data want is testable: estimate mu(x),
# bin on it, and compare the observed share of zeros per bin against the best
# curve each model can draw.
A <- path.expand("~/.claude/skills/live-progress/assets")
source(file.path(A, "progress.R"))

suppressMessages(library(bartisan))
data("lalonde", package = "cobalt")

OUT <- "_dev/tweedie-feasibility.rds"
y <- lalonde[["re78"]]
cat(sprintf("re78: n = %d, %.1f%% zeros, mean %.0f, sd %.0f\n",
            length(y), 100*mean(y == 0), mean(y), sd(y)))

# Marginal solve first: can any (phi, p) match the zero share and the variance?
solve_marg <- function(mu, v, p0) {
  f <- function(p) { phi <- v / mu^p; exp(-mu^(2-p)/(phi*(2-p))) - p0 }
  r <- tryCatch(stats::uniroot(f, c(1.001, 1.999)), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  list(p = r$root, phi = v / mu^r$root)
}
m <- solve_marg(mean(y), var(y), mean(y == 0))
cat(sprintf("\nmarginally, matching the mean, variance and zero share needs p = %.3f, phi = %.1f\n",
            m$p, m$phi))
cat(sprintf("  implied sd %.0f (observed %.0f), implied P(0) %.3f (observed %.3f)\n",
            sqrt(m$phi*mean(y)^m$p), sd(y),
            exp(-mean(y)^(2-m$p)/(m$phi*(2-m$p))), mean(y == 0)))

# Now the conditional constraint. mu(x) from a flexible fit of E[y|x].
pr <- prog_init(total = 1L, title = "Tweedie feasibility on re78", unit = "fit",
                kind = "simulation")
on.exit(prog_end(pr, "failed", "aborted before the fit finished"), add = TRUE)

set.seed(4)
t0 <- Sys.time()
fit <- bartisan(re78 ~ ., data = lalonde, family = gaussian(), chains = 4,
                num_burn = 400, num_draws = 800, verbose = FALSE)
prog_tick(pr, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
          label = "E[y | x]")
mu <- pmax(fitted(fit), 50)
b <- cut(mu, stats::quantile(mu, seq(0, 1, length.out = 9L)), include.lowest = TRUE)
obs <- tapply(y == 0, b, mean)
mub <- tapply(mu, b, mean)
nb <- tapply(y, b, length)

# Best Tweedie curve through those bins, and the best cutpoint curve.
tw_ss <- function(par) {
  p <- 1 + 0.998/(1+exp(-par[1])); phi <- exp(par[2])
  sum(nb * (obs - exp(-mub^(2-p)/(phi*(2-p))))^2)
}
cut_ss <- function(par) sum(nb * (obs - (1 - stats::plogis(par[1]*log(mub) - par[2])))^2)
o1 <- stats::optim(c(0, 5), tw_ss); o2 <- stats::optim(c(1, 5), cut_ss)
p_hat <- 1 + 0.998/(1+exp(-o1$par[1])); phi_hat <- exp(o1$par[2])

cat(sprintf("\nweighted SS of the zero share across 8 bins of mu(x):\n"))
cat(sprintf("  Tweedie   (p = %.3f, phi = %.0f):  %.5f\n", p_hat, phi_hat, o1$value))
cat(sprintf("  cutpoint  (slope %.2f, cut %.2f):  %.5f\n", o2$par[1], o2$par[2], o2$value))
cat(sprintf("\n%-14s %8s %8s %10s %10s\n", "bin mean mu", "n", "obs P(0)", "Tweedie", "cutpoint"))
for (i in seq_along(obs))
  cat(sprintf("%14.0f %8d %8.3f %10.3f %10.3f\n", mub[i], nb[i], obs[i],
              exp(-mub[i]^(2-p_hat)/(phi_hat*(2-p_hat))),
              1 - stats::plogis(o2$par[1]*log(mub[i]) - o2$par[2])))
saveRDS(list(marginal = m,
             bins = data.frame(
               mu = as.numeric(mub), n = as.numeric(nb),
               obs = as.numeric(obs),
               tweedie = exp(-as.numeric(mub)^(2 - p_hat) /
                               (phi_hat * (2 - p_hat))),
               cutpoint = 1 - stats::plogis(o2$par[1] * log(as.numeric(mub)) -
                                              o2$par[2])),
             tweedie = list(p = p_hat, phi = phi_hat, ss = o1$value),
             cutpoint = list(slope = o2$par[1], cut = o2$par[2],
                             ss = o2$value),
             complete = TRUE),
        OUT)

on.exit()
prog_end(pr, "done")
cat("wrote", OUT, "\n")

# --- Why the intractable density is not a blocker -----------------------------
#
# For 1 < p < 2 the Tweedie density has no closed form, but it factors as
#
#   log f(y; mu, phi, p) = (1/phi)[y mu^(1-p)/(1-p) - mu^(2-p)/(2-p)]
#                          + log W(y, phi, p)
#
# and `log W` contains no mu at all. Verified numerically against the compound
# Poisson-gamma series: log f minus the bracket agrees to 6 decimal places at
# mu = 0.8, 2.0 and 5.0. The bracket is the whole eta-dependent part, and
# `Family::compute_eta_free()` in src/family.h exists for exactly this: terms
# that "cancel from every acceptance ratio", refreshed when a nuisance parameter
# changes. So the series is evaluated once per sweep, not once per leaf.
#
# The series length depends on y, phi and p and not on mu, so it does not grow
# as the forest moves. On re78 with the marginally implied p = 1.433 and
# phi = 180.7 it needs at most 20 terms, the peak sitting at j = 1 to 5; at
# phi = 10 the worst observation needs 143. A saddlepoint fallback is wanted
# only for extreme y/phi.
#
# P(y = 0) is closed form, exp(-mu^(2-p)/(phi(2-p))), so the zero part of the
# likelihood needs no series at all.
#
# And the score and information are closed form from the exponential-dispersion
# structure, with mu = exp(eta):
#
#   dlogf/deta   = mu^(1-p) (y - mu) / phi
#   -d2logf/deta2 = mu^(2-p) / phi          (expected information; positive)
