## file = DM_rmse_All_3_simplex.R
## C:\Users\Jon.Brodziak\Documents\GitHub\DirichletStudy\R\scripts
## Compare accuracy of 3 DM forms over a K=3 simplex mesh
## Accuracy measure is the root mean square error of the
## estimated DM proportion based on a set of G independent
## samples from a single population with observed sample sizes
## drawn from a log-uniform(nmin, nmax) distribution
## 
## The three Dirichlet multinomial forms are
## (i)   unconstrained concentration vector alpha
## (ii)  constrained total concentration alpha0
## (iii) linear constraint on input sample size theta
##
## Input file:
##  - txt: "DM_rmse_All_3_simplex.inp"
## Output files:
##  - CSV: "DM_rmse_All_3_simplex.csv"
##  - Plots: "rmse_i_ternary.png","rmse_ii_ternary.png", "rmse_iii_ternary.png"
##  - List: "DM_rmse_All_3_simplex.lst"
##
## Output to console tests of whether significant differences exist
## Friedman test: Tests whether the three RMSE distributions differ 
## overall across matched simplex points using a nonparametric 
## repeated-measures ANOVA.
## Kendall’s W: Quantifies the effect size from the Friedman test, 
## indicating how consistently one method outperforms others across the grid (0–1 scale).
## Pairwise Wilcoxon signed-rank (Holm-adjusted): Pinpoints which 
## method pairs have significant RMSE differences while respecting the 
## within-point pairing and controlling familywise error.
## Hodges–Lehmann median difference with 95% CI: Estimates the 
## typical magnitude and direction of RMSE change between each pair 
## of methods with robust confidence intervals.
## Percent-change summaries (median, IQR): Expresses practical 
## impact as relative RMSE improvement/worsening between methods 
## across simplex points.
## Mixed-effects model (optional parametric check): Confirms results 
## by modeling method as a fixed effect and simplex point as a random 
## effect, with estimated marginal means for pairwise contrasts.

## ---- packages --------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(
  library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)))
suppressWarnings(suppressPackageStartupMessages(
  library(lme4, quietly = TRUE, warn.conflicts = FALSE)))
suppressWarnings(suppressPackageStartupMessages(
  library(lmerTest, quietly = TRUE, warn.conflicts = FALSE)))
suppressWarnings(suppressPackageStartupMessages(
  library(emmeans, quietly = TRUE, warn.conflicts = FALSE)))

## ---- simplex functions ----------------------------------------------------  

## ---- generate K compositions of N------------------------------------------
# nexcom(N, K) outputs a list containing the set of compositions
#              of the integer N into K parts. The composition list
#              is ordered lexicographically from the first composition
#              of P[first] = (N, 0, 0, ..., 0) to the last composition
#              of P[last] = (0, 0, 0, ..., N)
#
# function nexcom.step
# nexcom.step(N, K, P, MTC, I, J) returns the next K-dimensional composition 
#                                 vector of N, P[next], given the current 
#                                 composition vector of N, P[current].
# function arguments are:
#	N - The positive integer to compose into K parts.
#	K - The positive  integer number of parts, or categories in a composition vector. 
#	P - The current composition vector in the lexicographic set of all vectors.
#	MTC	- The logical flag indicating if the set of compositions of N is not complete
#       as in an acronym for "More To Come".
#   If MTC = TRUE, then the current composition is not the last composition
#                  in the lexicographic set of all vectors.
#		If MTC = FALSE, then the current composition is the last composition.
#	I - An index variable.
#	J - An index variable.
#
# returns - A list containing (P, MTC, I, J).
#
nexcom.step <- function (N, K, P, MTC, I, J) {
  
  if (MTC == FALSE) {
    P[1] <- N
    I <- N
    J <- 0	
    if (K != 1) {
      for (ii in 2:K) {
        P[ii] <- 0
      }
      MTC <- (P[K] != N)
      return(list(P = P, MTC = MTC, I = I, J = J))
    }
  }
  if (I > 1) J <- 0
  J <- J + 1
  I <- P[J]
  P[J] <- 0
  P[1] <- I - 1
  P[J + 1] <- P[J + 1] + 1
  MTC <- (P[K] != N)
  return(list(P = P, MTC = MTC, I = I, J = J))
}

# function nexcom
# nexcom() is the wrapper function that calls nexcom.step to step through
# the set of lexicographically-ordered K-compositions of N.
#
# function arguments are:
#	N - The positive integer to compose into K parts.
#	K - The number of parts of N, or categories in the composition. 
#
# returns - A data frame containing all possible K-part compositions of N.
#
nexcom <- function (N, K) {
  
  rn.comp <- nexcom.step(N, K, P = integer(K), MTC = FALSE, I = 0, J = 0)
  df.comp <- data.frame(P = rbind(rn.comp$P)) 
  
  ii <- 0
  while(rn.comp$MTC == TRUE) {
    rn.comp <- nexcom.step(N, K, P = rn.comp$P, MTC = rn.comp$MTC, I = rn.comp$I, J = rn.comp$J)
    df.comp <- rbind(df.comp, data.frame(P = rbind(rn.comp$P)))
    ii <- ii + 1
  }
  return(df.comp)
}

# --- mk_simplex() using nexcom(N, K) ----------------------------------------
# Inputs:
#   K : integer >= 2
#   h : numeric with 0 < h < 1/K
# Output:
#   Matrix with K columns (p1..pK), rows sum to 1.
#   attr(., "mesh_size") = 1/N, where N = ceiling(1/h).
mk_simplex <- function(K, h) {
  # validate
  if (!is.numeric(K) || length(K) != 1 || K != as.integer(K) || K < 2)
    stop("K must be a single integer >= 2.")
  if (!is.numeric(h) || length(h) != 1 || !(h > 0) || !(h < 1 / K))
    stop("h must be a single numeric with 0 < h < 1/K.")
  
  # choose grid so effective step <= h
  N <- as.integer(ceiling(1 / h))
  if (N < 1L) N <- 1L
  h_eff <- 1 / N
  
  # all K-part compositions of N (nonnegative integers summing to N)
  # nexcom(N, K) returns an N x K integer data frame
  ints <- nexcom(N, K)
  
  # scale to the open simplex and convert df to matrix
  simplex <- (ints + 1) / (N + K)
  simplex <- as.matrix(simplex)
}

## ---- utilities -------------------------------------------------------------
softmax_from_z <- function(z) { v <- c(z, 0); ev <- exp(v - max(v)); ev / sum(ev) }

logfile <- "DM_rmse_All_3_simplex.lst"
sink(logfile, split = TRUE, type = "output")

dm_loglik_group <- function(x, p, a0) {
  N <- sum(x)
  lgamma(a0) - lgamma(N + a0) + sum(lgamma(x + a0 * p) - lgamma(a0 * p))
}

# scores (i)
score_param_i <- function(t, X) {
  alpha_vec <- exp(t)
  a0    <- sum(alpha_vec)
  g_alpha <- rep(0.0, length(alpha_vec))
  for (g in 1:G) {
    xg <- X[g, ]
    Ng <- sum(xg)
    common <- digamma(a0) - digamma(Ng + a0)         # scalar
    g_alpha <- g_alpha +
      common + (digamma(xg + alpha_vec) - digamma(alpha_vec))# vector add
  }
  g_t <- alpha_vec * g_alpha
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

## ---- simulation settings ---------------------------------------------------
# ---- Read parameters from DM_rmse_All_3_simplex.inp ----
inp_path <- "DM_rmse_All_3_simplex.inp"

if (!file.exists(inp_path)) {
  stop(sprintf("Input file '%s' not found. Create it first.", inp_path))
}

# Parse key=value pairs, ignore comments and blank lines
tbl <- read.table(
  file = inp_path,
  sep = "=",
  comment.char = "#",
  strip.white = TRUE,
  blank.lines.skip = TRUE,
  col.names = c("key", "value"),
  colClasses = c("character", "character")
)

kv <- setNames(tbl$value, tbl$key)

# Set typed parameters
K          <- as.integer(kv[["K"]])
G          <- as.integer(kv[["G"]])
h          <- as.numeric(kv[["h"]])
theta_true <- as.numeric(kv[["theta_true"]])
Nmin       <- as.integer(kv[["Nmin"]])
Nmax       <- as.integer(kv[["Nmax"]])
nsamples       <- as.integer(kv[["nsamples"]])
random.seed       <- as.integer(kv[["random.seed"]])

# basic validation
if (any(is.na(c(K, G, h, theta_true, Nmin, Nmax, nsamples, random.seed)))) {
  stop("One or more required parameters are missing or not numeric in the .inp file.")
}
if (Nmin > Nmax) stop("Nmin must be <= Nmax.")

params <- list(K = K, G = G, h = h, theta_true = theta_true, Nmin = Nmin, 
               Nmax = Nmax, nsamples = nsamples, random.seed = random.seed)

set.seed(random.seed)

## ---- mesh of p_true -------------------------------------------------------
mesh <- mk_simplex(K, h)
nmesh <- nrow(mesh)

## Vary N_g over orders of magnitude
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
    alpha_vec <- a0g * p_true
    phi <- rgamma(K, shape = alpha_vec, rate = 1); phi <- phi / sum(phi)
    X[g, ] <- as.vector(rmultinom(1, size = N_row[g], prob = phi))
  }
 
## MLE (i): unconstrained alpha_vec (shared across groups) -------------------
## Start at alpha_vec = alpha0 * p_init with a small positive floor
alpha0_init <- 1.0
pooled <- colSums(X); p_init <- pmax(pooled, 1) / sum(pmax(pooled, 1))
t0_i <- log(pmax(alpha0_init * p_init, 1e-4))   # t = log(alpha_vec), length K

## Negative joint log-likelihood for (i)
nll_i <- function(t,X) {
  alpha_vec <- exp(t)                 # length K
  a0    <- sum(alpha_vec)
  # Sum over groups
  -sum(vapply(1:G, function(g) {
    xg <- X[g, ]
    Ng <- sum(xg)
    lgamma(a0) - lgamma(Ng + a0) +
      sum(lgamma(xg + alpha_vec) - lgamma(alpha_vec))
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
# Replace your block with this:

out <- as.data.frame(res_mat)

# K = number of true proportions assumed to be the first K columns of res_mat
K <- ncol(out) - 3L  # (p1..pK, then rmse_i, rmse_ii, rmse_iii)

# Ensure N_vec is a matrix with one row per simulation and G columns (groups)
N_mat <- if (is.null(dim(N_vec))) {
  matrix(N_vec, nrow = nrow(out), byrow = TRUE)
} else {
  N_vec
}

# Coerce sample sizes to integers (safe rounding, then integer cast)
N_int <- matrix(as.integer(round(N_mat)),
                nrow = nrow(N_mat),
                ncol = ncol(N_mat))
colnames(N_int) <- paste0("N", seq_len(ncol(N_int)))  # N1..NG

# Bind G sample-size columns to out
out <- cbind(out, N_int)

# Name the existing columns for proportions and RMSEs
colnames(out)[seq_len(K)] <- paste0("p", seq_len(K))
colnames(out)[K + seq_len(3)] <- c("rmse_i","rmse_ii","rmse_iii")

# Write CSV (no row names)
# write.csv(out, file = "DM_rmse_by_simulation.csv", row.names = FALSE)

## ---- write CSV -------------------------------------------------------------
csv_path <- "DM_rmse_All_3_simplex.csv"
write.csv(out, csv_path, row.names = FALSE)
message("Wrote: ", csv_path)

## ---- ternary plots ---------------------------------------------------------
if (K == 3) {
  
suppressWarnings(suppressPackageStartupMessages(
  library(ggtern, quietly = TRUE, warn.conflicts = FALSE)))
  
rmse_min <- min(out$rmse_i,out$rmse_ii, out$rmse_iii, na.rm = TRUE)
rmse_max <- max(out$rmse_i,out$rmse_ii, out$rmse_iii, na.rm = TRUE)

p_i <- ggtern::ggtern(out, aes(x = p1, y = p2, z = p3, colour = rmse_i)) +
  geom_point(shape = 16, size = 2, alpha = 0.9) +
  labs(title = sprintf("(i) RMSE over simplex (K=3, G=%d, θ=%.1f, unconstrained α)", G, theta_true), colour = "RMSE") +
  scale_colour_viridis_c(limits = c(rmse_min, rmse_max), oob = scales::squish) +
  theme_bw()

p_ii <- ggtern::ggtern(out, aes(x = p1, y = p2, z = p3, colour = rmse_ii)) +
  geom_point(shape = 16, size = 2, alpha = 0.9) +
  labs(title = sprintf("(ii) RMSE over simplex (K=3, G=%d, θ=%.1f, α0=beta)", G, theta_true), colour = "RMSE") +
  scale_colour_viridis_c(limits = c(rmse_min, rmse_max), oob = scales::squish) +
  theme_bw()

p_iii <- ggtern::ggtern(out, aes(x = p1, y = p2, z = p3, colour = rmse_iii)) +
  geom_point(shape = 16, size = 2, alpha = 0.9) +
  labs(title = sprintf("(iii) RMSE over simplex (K=3, G=%d, θ=%.1f, α0=θ*N)", G, theta_true), colour = "RMSE") +
    scale_colour_viridis_c(limits = c(rmse_min, rmse_max), oob = scales::squish) +
  theme_bw()

ggsave("rmse_i_ternary.png", p_i, width = 6, height = 5, dpi = 300)
ggsave("rmse_ii_ternary.png", p_ii, width = 6, height = 5, dpi = 300)
ggsave("rmse_iii_ternary.png", p_iii, width = 6, height = 5, dpi = 300)

print(p_i)
print(p_ii)
print(p_iii)

}

print(params)

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

on.exit(sink()) 
# on.exit(sink(), add = TRUE) 