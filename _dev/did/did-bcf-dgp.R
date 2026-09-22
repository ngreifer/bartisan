# The two DGPs from Souto & Louzada Neto (2025), arXiv:2505.09706.
# Cell A = their DGP 2 / Setting 1; Cell B = their DGP 5 / Setting 3.

make_panel <- function(N = 200L, Tt = 8L) {
  d <- expand.grid(t = seq_len(Tt), id = seq_len(N))
  d[order(d$id, d$t), c("id", "t")]
}

# DGP 2, Setting 1: staggered, homogeneous tau = 3, fully linear.
dgp_A <- function(seed, N = 200L, Tt = 8L, tau_mult = 1) {
  set.seed(seed)
  d <- make_panel(N, Tt); n <- nrow(d)
  d$x1 <- rbinom(n, 1L, 0.66)
  for (j in 2:6) d[[paste0("x", j)]] <- rnorm(n)
  d$x7 <- sample(1:4, n, TRUE, prob = c(0.3, 0.1, 0.2, 0.4))
  beta <- c(-0.75, 0.5, -0.5, -1.30, 1.8, 2.5, -1.0)
  X <- as.matrix(d[, paste0("x", 1:7)])

  g <- sample(rep(0:3, each = N / 4L))            # 50 units per group
  G <- c(Inf, 4, 5, 6)[g + 1L][d$id]
  d$G <- G
  d$cohort <- factor(ifelse(is.finite(G), as.character(G), "never"),
                     levels = c("never", "4", "5", "6"))
  d$ever <- as.integer(is.finite(G))
  d$k <- d$t - G
  d$Dit <- as.integer(is.finite(G) & d$t >= G)

  d$tau_true <- 3 * tau_mult
  d$y <- -0.5 + 0.75 * d$ever + 0.2 * d$t + as.vector(X %*% beta) +
    d$tau_true * d$Dit + rnorm(n)
  d$x7 <- factor(d$x7)
  d
}

# DGP 5, Setting 3: staggered + selection on statics, conditional heterogeneity,
# every covariate non-linear, quadratic time trend, tau_base = 5.
dgp_B <- function(seed, N = 200L, Tt = 8L, tau_mult = 1) {
  set.seed(seed)
  d <- make_panel(N, Tt); n <- nrow(d)

  x1_u <- rbinom(N, 1L, 0.66)                      # static, drives assignment
  x8_u <- rnorm(N)                                 # static, drives assignment
  d$x1 <- x1_u[d$id]
  d$x8 <- x8_u[d$id]
  d$x2 <- rbinom(n, 1L, 0.45)
  d$x3 <- sample(1:4, n, TRUE, prob = c(0.3, 0.1, 0.2, 0.4))   # drives CHTE
  d$x4 <- rnorm(n)                                             # drives CHTE
  for (j in 5:7) d[[paste0("x", j)]] <- rnorm(n)

  # Assignment: utility maximization over {never, 4, 5, 6} on the statics.
  V <- cbind(0,
             0.1 + 0.8 * x1_u + 0.6 * x8_u,
             -0.5 * x1_u - 0.7 * x8_u,
             0.1 + 0.3 * x1_u + 0.4 * x8_u)
  U <- V + matrix(rnorm(N * 4L, 0, 0.5), N, 4L)
  g <- max.col(U) - 1L
  G <- c(Inf, 4, 5, 6)[g + 1L][d$id]
  d$G <- G
  d$cohort <- factor(ifelse(is.finite(G), as.character(G), "never"),
                     levels = c("never", "4", "5", "6"))
  d$ever <- as.integer(is.finite(G))
  d$k <- d$t - G
  d$Dit <- as.integer(is.finite(G) & d$t >= G)

  # Non-linear baseline: exp on x1,x2; square on x3,x4; abs on x5,x6;
  # sqrt|.| on x7,x8; quadratic time.
  beta <- c(-0.75, 0.5, -0.5, -1.30, 1.8, 2.5, -1.0, 0.3)
  base <- exp(d$x1) * beta[1] + exp(d$x2) * beta[2] +
    d$x3^2 * beta[3] + d$x4^2 * beta[4] +
    abs(d$x5) * beta[5] + abs(d$x6) * beta[6] +
    sqrt(abs(d$x7)) * beta[7] + sqrt(abs(d$x8)) * beta[8]

  # tau_mult = 0 is the null cell: the paper notes that with no effect DGP 5 is
  # structurally DGP 3, so the heterogeneity vanishes with the effect.
  tau_base <- 5
  d$tau_true <- tau_mult *
    ifelse(d$x3 %in% c(1, 3), tau_base + 1.5 * sqrt(abs(d$x4)),
    ifelse(d$x3 == 2,         tau_base,
                              tau_base - 0.5 * sqrt(abs(d$x4))))

  d$y <- -0.5 + 0.75 * d$ever + 0.2 * d$t^2 + base +
    d$tau_true * d$Dit + rnorm(n)
  d$x3 <- factor(d$x3)
  d
}

covs_A <- paste0("x", 1:7)
covs_B <- paste0("x", 1:8)
