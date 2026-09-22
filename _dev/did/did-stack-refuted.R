# The discriminating test. The paper's DGPs give units no persistent
# heterogeneity at all -- only a common 0.75 shift for the ever-treated -- and
# their Eq. 3 has no unit effect to absorb one. Real panels do have them. Add
# a unit effect to their DGP 2 / Setting 1 and see which form survives.
suppressPackageStartupMessages(library(bartisan))
setwd("/Users/NoahGreifer/Dropbox/Research/R/bartisan"); source("_dev/did/did-bcf-dgp.R")
source("_dev/did/did-stack-fun.R")

run <- function(sigma_u, reps = 5L) {
  cat(sprintf("\n== unit-effect sd = %.1f (paper's DGP is 0) ==\n", sigma_u))
  for (r in seq_len(reps)) {
    d <- dgp_A(seed = 9000L + r)
    set.seed(9000L + r + 500L)
    u <- rnorm(length(unique(d$id)), 0, sigma_u)
    d$y <- d$y + u[d$id]
    d$G2 <- ifelse(is.finite(d$G), d$G, Inf); d$k_fin <- ifelse(is.finite(d$k), d$k, -99)
    tr <- d$Dit == 1
    fl <- suppressMessages(suppressWarnings(bartisan(
      y ~ cohort + t + x1+x2+x3+x4+x5+x6+x7 +
        vc(Dit, ~ x1+x2+x3+x4+x5+x6+x7 + k_fin), data = d, family = gaussian(),
      control = bartisan_control(num_trees = c(150L,25L), num_burn = 400L,
                                 num_draws = 600L, chains = 2L, verbose = FALSE))))
    al <- rowMeans(coef(fl, draws = TRUE)[[1L]][, tr, drop = FALSE])

    s <- stack_cells(d, "id", "t", "y", "G2", paste0("x", 1:7))
    fs <- suppressMessages(suppressWarnings(bartisan(
      dY ~ cell + x1+x2+x3+x4+x5+x6+x7 +
        vc(treat, ~ x1+x2+x3+x4+x5+x6+x7 + cohort + since_f), data = s,
      family = gaussian(),
      control = bartisan_control(num_trees = c(50L,25L), num_burn = 400L,
                                 num_draws = 600L, chains = 2L, verbose = FALSE))))
    trs <- s$treat == 1
    as_ <- rowMeans(coef(fs, draws = TRUE)[[1L]][, trs, drop = FALSE])
    cat(sprintf("  rep %d  level %+.3f [%+.3f,%+.3f] cov %-5s | stacked %+.3f [%+.3f,%+.3f] cov %s\n",
        r, mean(al), quantile(al,.025), quantile(al,.975),
        quantile(al,.025) <= 3 && 3 <= quantile(al,.975),
        mean(as_), quantile(as_,.025), quantile(as_,.975),
        quantile(as_,.025) <= 3 && 3 <= quantile(as_,.975)))
    utils::flush.console()
  }
}
run(0); run(2)
