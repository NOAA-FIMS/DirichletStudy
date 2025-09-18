## file = DM_saturated_vs_linear_Ng_simplex.R
## DM curvature/accuracy comparison over a K=3 simplex mesh
## Outputs:
##  - CSV: "dm_rmse_simplex_mesh.csv"
##  - Plots: "rmse_ii_ternary.png", "rmse_iii_ternary.png"

## ---- packages --------------------------------------------------------------
if (!requireNamespace("ggtern", quietly = TRUE)) install.packages("ggtern")
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtern)
  library(patchwork)
})

## ---- simplex mesh (mk_simplex.r style) ------------------------------------
mk_simplex <- function(h = 0.05, interior = TRUE) {
  M <- round(1 / h)
  out <- vector("list", (M + 1) * (M + 2) / 2)
  k <- 0L
  for (i in 0:M) for (j in 0:(M - i)) {
    k3 <- M - i - j
    if (interior && (i == 0 || j == 0 || k3 == 0)) next
    k <- k + 1L
    out[[k]] <- c(i, j, k3) / M
  }
  do.call(rbind, out[seq_len(k)])
}

## ---- utilities -------------------------------------------------------------
softmax_from_z <- function(z) { v <- c(z, 0); ev <- exp(v - max(v)); ev / sum(ev) }

dm_loglik_group <- function(x, p, a0) {
  N <- sum(x)
  lgamma(a0) - lgamma(N + a0) + sum(lgamma(x + a0 * p) - lgamma(a0 * p))
}

## scores (ii): p shared; a0_g free (opt over z (length K-1) and t_g = log a0_g)
score_param_ii <- function(X, p, a0_vec) {
  G <- nrow(X); K <- length(p)
  g_p <- numeric(K); g_t <- numeric(G)
  for (g in 1:G) {
    xg <- X[g, ]; a0 <- a0_vec[g]; Ng <- sum(xg)
    g_p <- g_p + a0 * (digamma(xg + a0 * p) - digamma(a0 * p))
    g_t[g] <- a0 * ( digamma(a0) - digamma(Ng + a0)
                     + sum(p * (digamma(xg + a0 * p) - digamma(a0 * p))) )
  }
  J <- diag(p) - tcrossprod(p)                     # projection to simplex tangent
  g_z <- as.vector(crossprod(J[, 1:(K-1), drop=FALSE], g_p))
  list(g_z = g_z, g_t = g_t)
}

## scores (iii): p shared; theta shared (opt over z and t = log theta)
score_param_iii <- function(X, p, theta) {
  G <- nrow(X); K <- length(p)
  g_p <- numeric(K); g_theta <- 0
  for (g in 1:G) {
    xg <- X[g, ]; Ng <- sum(xg); a0 <- theta * Ng
    g_p <- g_p + a0 * (digamma(xg + a0 * p) - digamma(a0 * p))
    g_theta <- g_theta + Ng * ( digamma(a0) - digamma(Ng + a0)
                                + sum(p * (digamma(xg + a0 * p) - digamma(a0 * p))) )
  }
  J <- diag(p) - tcrossprod(p)
  g_z <- as.vector(crossprod(J[, 1:(K-1), drop=FALSE], g_p))
  g_t <- theta * g_theta
  list(g_z = g_z, g_t = g_t)
}

rmse <- function(a, b) sqrt(mean((a - b)^2))

## ---- simulation settings ---------------------------------------------------
K <- 3
G <- 10L
theta_true <- 0.2
Nmin  <- 25
Nmax  <- 250

set.seed(26267)

## ---- mesh of p_true --------------------------------------------------------
mesh <- mk_simplex(h = 0.05, interior = TRUE)  # exclude boundaries to keep α0 p_k > 0
nmesh <- nrow(mesh)

## Vary N_g over orders of magnitude (kept fixed across mesh points for comparability)
N_vec <- matrix(
  round(exp(runif(nmesh * G, log(Nmin), log(Nmax)))),
  nrow = nmesh, ncol = G, byrow = TRUE
)

## ---- loop over mesh --------------------------------------------------------
run_one <- function(p_true,N_row) {
  ## simulate counts for this p_true
  X <- matrix(0L, nrow = G, ncol = K)
  for (g in 1:G) {
    a0g <- theta_true * N_row[g]
    alpha <- a0g * p_true
    phi <- rgamma(K, shape = alpha); phi <- phi / sum(phi)
    X[g, ] <- as.vector(rmultinom(1, size = N_row[g], prob = phi))
  }
  
  ## MLE (ii): (p, a0_g)
  pooled <- colSums(X); p_init <- pmax(pooled, 1) / sum(pmax(pooled, 1))
  z0 <- log(p_init[1:(K-1)]) - log(p_init[K])
  t0 <- rep(0, G)  # log a0_g = 0 => a0_g = 1
  
  nll_ii <- function(par) {
    z <- par[1:(K-1)]; t <- par[K:(K-1+G)]
    p <- softmax_from_z(z); a0 <- exp(t)
    -sum(vapply(1:G, function(g)
      dm_loglik_group(X[g,], p, a0[g]), 0.0))
  }
  grad_ii <- function(par) {
    z <- par[1:(K-1)]; t <- par[K:(K-1+G)]
    p <- softmax_from_z(z); a0 <- exp(t)
    sc <- score_param_ii(X, p, a0)
    -c(sc$g_z, sc$g_t)
  }
  fit_ii <- try(optim(c(z0, t0), nll_ii, grad_ii, method = "BFGS",
                      control = list(maxit = 500, reltol = 1e-10)), silent = TRUE)
  
  if (inherits(fit_ii, "try-error")) {
    p_hat_ii <- rep(NA_real_, K)
  } else {
    zhat <- fit_ii$par[1:(K-1)]
    p_hat_ii <- softmax_from_z(zhat)
  }
  
  ## MLE (iii): (p, theta), α0_g = θ N_g
  nll_iii <- function(par) {
    z <- par[1:(K-1)]; t <- par[K]
    p <- softmax_from_z(z); theta <- exp(t)
    -sum(vapply(1:G, function(g)
      dm_loglik_group(X[g,], p, theta * N_vec[g]), 0.0))
  }
  grad_iii <- function(par) {
    z <- par[1:(K-1)]; t <- par[K]
    p <- softmax_from_z(z); theta <- exp(t)
    sc <- score_param_iii(X, p, theta)
    -c(sc$g_z, sc$g_t)
  }
  fit_iii <- try(optim(c(z0, 0), nll_iii, grad_iii, method = "BFGS",
                       control = list(maxit = 500, reltol = 1e-10)), silent = TRUE)
  
  if (inherits(fit_iii, "try-error")) {
    p_hat_iii <- rep(NA_real_, K)
  } else {
    zhat <- fit_iii$par[1:(K-1)]
    p_hat_iii <- softmax_from_z(zhat)
  }
  rmse_ii  = rmse(p_hat_ii,  p_true)
  rmse_iii = rmse(p_hat_iii, p_true)
  c(p_true,rmse_ii,rmse_iii)
}

# res_mat <- t(apply(mesh, 1, run_one))
res_mat <- do.call(
  rbind,
  lapply(seq_len(nmesh), function(i) {
    run_one(p_true = mesh[i, ], N_row = N_vec[i, ])
  })
)
out <- as.data.frame(res_mat)
# N_str <- paste(N_vec, collapse = "-")
# out$N_vec <- N_str
N_str <- apply(N_vec, 1, function(v) paste(v, collapse = "-"))
out$N_sizes <- N_str
names(out) <- c("p1","p2","p3","rmse_ii","rmse_iii","Sample size by group")

## ---- write CSV -------------------------------------------------------------
csv_path <- "dm_rmse_simplex_mesh.csv"
write.csv(out, csv_path, row.names = FALSE)
message("Wrote: ", csv_path)

## ---- ternary plots ---------------------------------------------------------

rmse_min <- min(out$rmse_ii, out$rmse_iii, na.rm = TRUE)
rmse_max <- max(out$rmse_ii, out$rmse_iii, na.rm = TRUE)

p_ii <- ggtern::ggtern(out, aes(x = p1, y = p2, z = p3, colour = rmse_ii)) +
  geom_point(shape = 16, size = 2, alpha = 0.9) +
  labs(title = "(ii) RMSE over simplex (K=3, G=10, θ=0.2, α0=beta)", colour = "RMSE") +
  scale_colour_viridis_c(limits = c(rmse_min, rmse_max), oob = scales::squish) +
  theme_bw()

p_iii <- ggtern::ggtern(out, aes(x = p1, y = p2, z = p3, colour = rmse_iii)) +
  geom_point(shape = 16, size = 2, alpha = 0.9) +
  labs(title = "(iii) RMSE over simplex (K=3, G=10, θ=0.2, α0=θ*N)", colour = "RMSE") +
  scale_colour_viridis_c(limits = c(rmse_min, rmse_max), oob = scales::squish) +
  theme_bw()

ggsave("rmse_ii_ternary.png", p_ii, width = 6, height = 5, dpi = 300)
ggsave("rmse_iii_ternary.png", p_iii, width = 6, height = 5, dpi = 300)

print(p_ii)
print(p_iii)
