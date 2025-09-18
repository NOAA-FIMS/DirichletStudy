## file = DM_rmse_All_3_simplex.R
## Compare accuracy of 3 DM forms over a K=3 simplex mesh
## Accuracy measure is the root mean square error of the
## estimated DM proportion based on a set of G independent
## samples from a single population with observed sample sizes
## drawn from a loguniform(nmin, nmax) distribution
## 
## The three Dirichlet multinomial forms are
## (i)   unconstrained concentration vector alpha
## (ii)  constrained total concentration alpha0
## (iii) linear constraint on input sample size theta
##
## Output files:
##  - CSV: "dm_rmse_simplex_mesh.csv"
##  - Plots: 
##    "rmse_i_ternary.png","rmse_ii_ternary.png", "rmse_iii_ternary.png"
##
## Output to console tests of whether significant differences exist
## Friedman test:** Tests whether the three RMSE distributions differ 
## overall across matched simplex points using a nonparametric 
## repeated-measures ANOVA.
## Kendall’s W:** Quantifies the effect size from the Friedman test, 
## indicating how consistently one method outperforms others across the grid (0–1 scale).
## Pairwise Wilcoxon signed-rank (Holm-adjusted):** Pinpoints which 
## method pairs have significant RMSE differences while respecting the 
## within-point pairing and controlling familywise error.
## Hodges–Lehmann median difference with 95% CI:** Estimates the 
## typical magnitude and direction of RMSE change between each pair 
## of methods with robust confidence intervals.
## Percent-change summaries (median, IQR):** Expresses practical 
## impact as relative RMSE improvement/worsening between methods 
## across simplex points.
## Mixed-effects model (optional parametric check):** Confirms results 
## by modeling method as a fixed effect and simplex point as a random 
## effect, with estimated marginal means for pairwise contrasts.

## ---- packages --------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(
  library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
))
suppressWarnings(suppressPackageStartupMessages(
  library(ggtern, quietly = TRUE, warn.conflicts = FALSE)
))
suppressWarnings(suppressPackageStartupMessages(
  library(patchwork, quietly = TRUE, warn.conflicts = FALSE)
))
suppressWarnings(suppressPackageStartupMessages(
  library(lme4, quietly = TRUE, warn.conflicts = FALSE)
))
suppressWarnings(suppressPackageStartupMessages(
  library(lmerTest, quietly = TRUE, warn.conflicts = FALSE)
))
suppressWarnings(suppressPackageStartupMessages(
  library(emmeans, quietly = TRUE, warn.conflicts = FALSE)
))

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

# scores (i)
score_param_i <- function(t, X) {
  alpha <- exp(t)
  a0    <- sum(alpha)
  g_alpha <- rep(0.0, length(alpha))
  for (g in 1:G) {
    xg <- X[g, ]
    Ng <- sum(xg)
    common <- digamma(a0) - digamma(Ng + a0)         # scalar
    g_alpha <- g_alpha +
      common + (digamma(xg + alpha) - digamma(alpha))# vector add
  }
  g_t <- alpha * g_alpha
  -g_t   # gradient of nll
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

logfile <- "DM_rmse_All_3_simplex.lst"
## sink(logfile, split = TRUE)
## sink(logfile, split = TRUE)
sink(logfile, split = TRUE, type = "output")


## ---- simulation settings ---------------------------------------------------
K <- 3
G <- 4L
theta_true <- 0.2
Nmin  <- 25
Nmax  <- 250

set.seed(644)

## ---- mesh of p_true --------------------------------------------------------
mesh <- mk_simplex(h = 0.02, interior = TRUE)  # exclude boundaries to keep α0 p_k > 0
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
 
## MLE (i): unconstrained alpha (shared across groups) -------------------
## Start at alpha = alpha0 * p_init with a small positive floor
alpha0_init <- 1.0
pooled <- colSums(X); p_init <- pmax(pooled, 1) / sum(pmax(pooled, 1))
t0_i <- log(pmax(alpha0_init * p_init, 1e-4))   # t = log(alpha), length K

## Negative joint log-likelihood for (i)
nll_i <- function(t,X) {
  alpha <- exp(t)                 # length K
  a0    <- sum(alpha)
  # Sum over groups
  -sum(vapply(1:G, function(g) {
    xg <- X[g, ]
    Ng <- sum(xg)
    lgamma(a0) - lgamma(Ng + a0) +
      sum(lgamma(xg + alpha) - lgamma(alpha))
  }, 0.0))
}
 
fit_i <- optim(t0_i, nll_i, score_param_i, method = "BFGS",
               control = list(maxit = 500, reltol = 1e-10), 
               X = X)

alpha_hat_i <- exp(fit_i$par)
p_hat_i     <- alpha_hat_i / sum(alpha_hat_i)

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
  rmse_i   = rmse(p_hat_i,   p_true)
  rmse_ii  = rmse(p_hat_ii,  p_true)
  rmse_iii = rmse(p_hat_iii, p_true)
  c(p_true,rmse_i,rmse_ii,rmse_iii)
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
names(out) <- c("p1","p2","p3","rmse_i","rmse_ii","rmse_iii","Sample size by group")

## ---- write CSV -------------------------------------------------------------
csv_path <- "DM_rmse_All_3_simplex.csv"
write.csv(out, csv_path, row.names = FALSE)
message("Wrote: ", csv_path)

## ---- ternary plots ---------------------------------------------------------

rmse_min <- min(out$rmse_i,out$rmse_ii, out$rmse_iii, na.rm = TRUE)
rmse_max <- max(out$rmse_i,out$rmse_ii, out$rmse_iii, na.rm = TRUE)

p_i <- ggtern::ggtern(out, aes(x = p1, y = p2, z = p3, colour = rmse_i)) +
  geom_point(shape = 16, size = 2, alpha = 0.9) +
  labs(title = sprintf("(i) RMSE over simplex (K=3, G=%d, θ=%.1f, unconstrained α)", G, theta_true), colour = "RMSE") +
  scale_colour_viridis_c(limits = c(rmse_min, rmse_max), oob = scales::squish) +
  theme_bw()

p_ii <- ggtern::ggtern(out, aes(x = p1, y = p2, z = p3, colour = rmse_ii)) +
  geom_point(shape = 16, size = 2, alpha = 0.9) +
  labs(title = sprintf("(i) RMSE over simplex (K=3, G=%d, θ=%.1f, α0=beta)", G, theta_true), colour = "RMSE") +
  scale_colour_viridis_c(limits = c(rmse_min, rmse_max), oob = scales::squish) +
  theme_bw()

p_iii <- ggtern::ggtern(out, aes(x = p1, y = p2, z = p3, colour = rmse_iii)) +
  geom_point(shape = 16, size = 2, alpha = 0.9) +
  labs(title = sprintf("(i) RMSE over simplex (K=3, G=%d, θ=%.1f, α0=θ*N)", G, theta_true), colour = "RMSE") +
    scale_colour_viridis_c(limits = c(rmse_min, rmse_max), oob = scales::squish) +
  theme_bw()

ggsave("rmse_i_ternary.png", p_i, width = 6, height = 5, dpi = 300)
ggsave("rmse_ii_ternary.png", p_ii, width = 6, height = 5, dpi = 300)
ggsave("rmse_iii_ternary.png", p_iii, width = 6, height = 5, dpi = 300)

print(p_i)
print(p_ii)
print(p_iii)

n <- length(out$rmse_i)
stopifnot(length(out$rmse_ii) == n, length(out$rmse_iii) == n)

df <- data.frame(
  id     = rep(seq_len(n), times = 3),                 # block (simplex point)
  method = factor(rep(c("i","ii","iii"), each = n),
                  levels = c("i","ii","iii")),
  rmse   = c(out$rmse_i, out$rmse_ii, out$rmse_iii)
)

## Overall nonparametric repeated-measures test
ft <- friedman.test(rmse ~ method | id, data = df)
print(ft)

## Kendall’s W effect size (approx. from Friedman chi-square)
k <- nlevels(df$method)
W <- as.numeric(ft$statistic) / (n * (k - 1))
cat(sprintf("Kendall's W ≈ %.3f\n", W))

## Post-hoc paired comparisons with multiplicity control
pw <- pairwise.wilcox.test(df$rmse, df$method,
                           paired = TRUE, p.adjust.method = "holm",
                           exact = FALSE)
print(pw)

## Assuming df as before (id, method ∈ {i,ii,iii}, rmse)
wide <- reshape(df, idvar = "id", timevar = "method", direction = "wide")

## 1) Paired, robust effect: Hodges–Lehmann (median) difference + 95% CI
W2_1 <- wilcox.test(wide$rmse.ii, wide$rmse.i,  paired = TRUE, conf.int = TRUE)   # ii - i
W3_1 <- wilcox.test(wide$rmse.iii, wide$rmse.i, paired = TRUE, conf.int = TRUE)   # iii - i
W3_2 <- wilcox.test(wide$rmse.iii, wide$rmse.ii,paired = TRUE, conf.int = TRUE)   # iii - ii
print(W2_1)
print(W3_1)
print(W3_2)

## 2) Practical impact: percent change in RMSE (median & IQR)
pct <- transform(wide,
                 pct_ii_i   = 100*(rmse.ii - rmse.i)/rmse.i,
                 pct_iii_i  = 100*(rmse.iii - rmse.i)/rmse.i,
                 pct_iii_ii = 100*(rmse.iii - rmse.ii)/rmse.ii
)

tmp <- sapply(pct[c("pct_ii_i","pct_iii_i","pct_iii_ii")],
       \(x) c(median = median(x), IQR = IQR(x), q25 = quantile(x, .25), q75 = quantile(x, .75)))
print(tmp)

## 3) Standardized paired effect r from Wilcoxon Z
r_from_wilcox <- function(x, y){
  wt <- wilcox.test(x, y, paired=TRUE, exact=FALSE)
  # Approximate Z from p-value and sign of statistic
  z  <- qnorm(wt$p.value/2, lower.tail = FALSE) * sign(median(y - x))
  r  <- as.numeric(z)/sqrt(sum(!is.na(x) & !is.na(y)))
  return(r)
}
Rho_1_2 <- r_from_wilcox(wide$rmse.i, wide$rmse.ii)
Rho_1_3 <- r_from_wilcox(wide$rmse.i, wide$rmse.iii)
Rho_2_3 <- r_from_wilcox(wide$rmse.ii, wide$rmse.iii)
print(Rho_1_2)
print(Rho_1_3)
print(Rho_2_3)

m <- lmer(rmse ~ method + (1 | id), data = df)
print(anova(m))                               # overall method effect)
print(emmeans(m, pairwise ~ method, adjust = "holm"))  # post-hoc

on.exit(sink(), add = TRUE) 