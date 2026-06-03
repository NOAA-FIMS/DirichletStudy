## file = hake_multinomial_om.R
## Runs nsims iterations per simplex point

## ---- packages --------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(
  library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)))
suppressWarnings(suppressPackageStartupMessages(
  library(lme4, quietly = TRUE, warn.conflicts = FALSE)))
suppressWarnings(suppressPackageStartupMessages(
  library(lmerTest, quietly = TRUE, warn.conflicts = FALSE)))
suppressWarnings(suppressPackageStartupMessages(
  library(emmeans, quietly = TRUE, warn.conflicts = FALSE)))
suppressWarnings(suppressPackageStartupMessages(
  library(future, quietly = TRUE, warn.conflicts = FALSE)))
suppressWarnings(suppressPackageStartupMessages(
  library(future.apply, quietly = TRUE, warn.conflicts = FALSE)))

## ---- simplex_output_flag = 1 creates two simplex files --------------------
## file1 = outprefix.csv stores results by simplex point
## file2 = outprefix_simplex_samples.csv stores the sampled simplex points ---
simplex_output_flag <- 1

## ---- ggplot_output_flag = 1 creates all ggplot2 *.png files ---------------
ggplot_output_flag <- 1

## ---- simplex functions ----------------------------------------------------  
nexcom.step <- function (N, K, P, MTC, I, J) {
  if (MTC == FALSE) {
    P[1] <- N
    I <- N
    J <- 0	
    if (K != 1) {
      for (ii in 2:K) { P[ii] <- 0 }
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

mk_simplex <- function(K, h) {
  if (!is.numeric(K) || length(K) != 1 || K != as.integer(K) || K < 2)
    stop("K must be a single integer >= 2.")
  if (!is.numeric(h) || length(h) != 1 || !(h > 0) || !(h < 1 / K))
    stop("h must be a single numeric with 0 < h < 1/K.")
  N <- as.integer(ceiling(1 / h))
  if (N < 1L) N <- 1L
  h_eff <- 1 / N
  ints <- nexcom(N, K)
  simplex <- (ints + 1) / (N + K)
  simplex <- as.matrix(simplex)
}

## ---- utilities -------------------------------------------------------------
softmax_from_z <- function(z) { v <- c(z, 0); ev <- exp(v - max(v)); ev / sum(ev) }

dm_loglik_group <- function(x, p, a0) {
  N <- sum(x)
  lgamma(a0) - lgamma(N + a0) + sum(lgamma(x + a0 * p) - lgamma(a0 * p))
}

score_param_i <- function(t, X) {
  alpha_vec <- exp(t)
  a0    <- sum(alpha_vec)
  g_alpha <- rep(0.0, length(alpha_vec))
  for (g in 1:nrow(X)) {
    xg <- X[g, ]
    Ng <- sum(xg)
    common <- digamma(a0) - digamma(Ng + a0)
    g_alpha <- g_alpha + common + (digamma(xg + alpha_vec) - digamma(alpha_vec))
  }
  g_t <- alpha_vec * g_alpha
  -g_t
}

score_param_ii <- function(X, p, a0_vec) {
  G <- nrow(X); K <- length(p)
  g_p <- numeric(K); g_t <- numeric(G)
  for (g in 1:G) {
    xg <- X[g, ]; a0 <- a0_vec[g]; Ng <- sum(xg)
    g_p <- g_p + a0 * (digamma(xg + a0 * p) - digamma(a0 * p))
    g_t[g] <- a0 * ( digamma(a0) - digamma(Ng + a0)
                     + sum(p * (digamma(xg + a0 * p) - digamma(a0 * p))) )
  }
  J <- diag(p) - tcrossprod(p)
  g_z <- as.vector(crossprod(J[, 1:(K-1), drop=FALSE], g_p))
  list(g_z = g_z, g_t = g_t)
}

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
L1_norm <- function(a, b) sum(abs(a - b))
Linf_norm <- function(a, b) max(abs(a - b))

## ---- simulation settings ---------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
inp_path <- if (length(args) >= 1) args[1] else "hake_multinomial_om.inp"
if (!file.exists(inp_path)) stop(sprintf("Input file '%s' not found.", inp_path))

tbl <- read.table(
  file = inp_path, sep = "=", comment.char = "#", strip.white = TRUE,
  blank.lines.skip = TRUE, col.names = c("key", "value"), colClasses = c("character", "character")
)
kv <- setNames(tbl$value, tbl$key)

K          <- as.integer(kv[["K"]])
G          <- as.integer(kv[["G"]])
h          <- as.numeric(kv[["h"]])
theta_true <- as.numeric(kv[["theta_true"]])
theta_CV   <- as.numeric(kv[["theta_CV"]])
Nmin       <- as.integer(kv[["Nmin"]])
Nmax       <- as.integer(kv[["Nmax"]])
nsims      <- as.integer(kv[["nsims"]])
random.seed<- as.integer(kv[["random.seed"]])
od_mult    <- as.numeric(kv[["od_mult"]])
sigma      <- as.numeric(kv[["sigma"]])
dist_code  <- as.integer(kv[["dist_code"]])
mean_nsamp <- as.numeric(kv[["mean_nsamp"]])
ln_sd      <- as.numeric(kv[["ln_sd"]])
nb_size    <- as.numeric(kv[["nb_size"]])
multinomial_om_flag_raw <- kv[["multinomial_om_flag"]]
multinomial_om_flag <- if (
  is.null(multinomial_om_flag_raw) ||
  is.na(multinomial_om_flag_raw) ||
  !nzchar(trimws(multinomial_om_flag_raw))
) 0L else as.integer(multinomial_om_flag_raw)
if (!multinomial_om_flag %in% c(0L, 1L)) stop("multinomial_om_flag must be 0 or 1.")

p_ubound_raw <- kv[["p_ubound"]]
if (is.null(p_ubound_raw) || is.na(p_ubound_raw) || !nzchar(trimws(p_ubound_raw))) {
  p_ubound <- rep(1, K)
} else {
  p_ubound <- as.numeric(unlist(strsplit(p_ubound_raw, "[,[:space:]]+")))[nzchar(unlist(strsplit(p_ubound_raw, "[,[:space:]]+")))]
}

p_lbound_raw <- kv[["p_lbound"]]
if (is.null(p_lbound_raw) || is.na(p_lbound_raw) || !nzchar(trimws(p_lbound_raw))) {
  p_lbound <- rep(0, K)
} else {
  p_lbound <- as.numeric(unlist(strsplit(p_lbound_raw, "[,[:space:]]+")))[nzchar(unlist(strsplit(p_lbound_raw, "[,[:space:]]+")))]
}

params <- list(K=K, G=G, h=h, theta_true=theta_true, theta_CV=theta_CV, dist_code=dist_code,
               mean_nsamp=mean_nsamp, ln_sd=ln_sd, nb_size=nb_size, Nmin=Nmin, Nmax=Nmax,
               nsims=nsims, random.seed=random.seed, od_mult=od_mult, sigma=sigma,
               multinomial_om_flag=multinomial_om_flag,
               p_lbound=p_lbound, p_ubound=p_ubound)

set.seed(random.seed)

dist_levels <- c("Poisson", "Lognormal", "Negative binomial", "Log-uniform")
dist_tag    <- dist_levels[[as.integer(dist_code)]]
base_stem   <- tools::file_path_sans_ext(basename(inp_path))
out_prefix  <- "hake_multinomial_om"

local({
  logfile <- paste0(out_prefix, ".lst")
  base_output_sink <- sink.number(type = "output")

  close_log_sink <- function() {
    while (sink.number(type = "output") > base_output_sink) {
      sink(type = "output")
    }
    invisible(NULL)
  }

  sink(logfile, split = TRUE, type = "output")
  on.exit(close_log_sink(), add = TRUE)

  message("Logging to: ", logfile)
  start_time <- Sys.time()

  cat("Run Time Start:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

## ---- mesh of p_true -------------------------------------------------------
mesh <- mk_simplex(K, h)

## Filter mesh based on bounds FIRST
keep_idx <- apply(sweep(mesh, 2, p_ubound, `<=`), 1, all) & apply(sweep(mesh, 2, p_lbound, `>=`), 1, all)
mesh <- mesh[keep_idx, , drop = FALSE]
nmesh <- nrow(mesh)

cat("Sampled simplex dimensions:", dim(mesh), "\n\n")
if (simplex_output_flag == 1)
  write.csv(mesh, file = "simplex_samples.csv", row.names = FALSE)

## Calculate total runs using nsims
total_runs <- nmesh * nsims

## ---- sample sizes (N_vec) ------------------------------------------------
draw_nsamp <- function(n, dist_code, mean_nsamp, ln_sd, nb_size, Nmin, Nmax) {
  if (dist_code == 1L) return(pmax(1L, as.integer(rpois(n, lambda = mean_nsamp))))
  if (dist_code == 2L) return(pmax(1L, as.integer(round(rlnorm(n, meanlog = log(mean_nsamp) - 0.5 * ln_sd^2, sdlog = ln_sd)))))
  if (dist_code == 3L) return(pmax(1L, as.integer(rnbinom(n, size = nb_size, mu = mean_nsamp))))
  if (dist_code == 4L) {
    mu0 <- (Nmax - Nmin) / (log(Nmax) - log(Nmin))
    s <- mean_nsamp / mu0
    return(pmax(1L, as.integer(round(exp(runif(n, log(max(1L, as.integer(round(Nmin * s)))), log(max(1L + max(1L, as.integer(round(Nmin * s))), as.integer(round(Nmax * s))))))))))
  }
}

# Expand mesh and tracking IDs
mesh_expanded <- mesh[rep(seq_len(nmesh), each = nsims), , drop = FALSE]
mesh_id_vec <- rep(seq_len(nmesh), each = nsims)
sim_id_vec  <- rep(seq_len(nsims), times = nmesh)

# Draw sample sizes for the fully expanded simulation grid
N_draws <- draw_nsamp(total_runs * G, dist_code, mean_nsamp, ln_sd, nb_size, Nmin, Nmax)
N_vec <- matrix(N_draws, nrow = total_runs, ncol = G, byrow = TRUE)

cat("--- Sample-size distribution diagnostics ---\n")
cat(sprintf("Sample Size Distribution = %d (%s)\nTarget mean = %.3f\nRealized mean = %.3f\nRealized SD = %.3f\n", 
            dist_code, dist_tag, mean_nsamp, mean(as.vector(N_vec)), sd(as.vector(N_vec))))
print(quantile(as.vector(N_vec), probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)))
cat("--- End diagnostics ---\n\n")

theta_sigma <- sqrt(log(theta_CV^2 + 1))

## ---- loop over total runs --------------------------------------------------
run_one <- function(p_true, N_row) {
  X <- matrix(0L, nrow = G, ncol = K)
  for (g in 1:G) {
    if (multinomial_om_flag == 1L) {
      z <- log(p_true) + rnorm(K, 0, sigma)
      p_g <- exp(z)
      p_g <- p_g / sum(p_g)
      X[g, ] <- as.vector(rmultinom(1, size = N_row[g], prob = p_g))
    } else {
      theta_g <- theta_true * exp(rnorm(1, mean = 0, sd = theta_sigma) - (theta_sigma^2) / 2)
      a0g <- (theta_g * od_mult) * N_row[g]  
      z <- log(p_true) + rnorm(K, 0, sigma) 
      p_g <- exp(z)
      p_g <- p_g / sum(p_g)
      alpha_vec <- a0g * p_g
      phi <- rgamma(K, shape = alpha_vec, rate = 1); phi <- phi / sum(phi)
      X[g, ] <- as.vector(rmultinom(1, size = N_row[g], prob = phi))
    }
  }
  X_out <- as.vector(t(cbind(N_row, X)))
  
  ## MLE (i)
  alpha0_init <- 1.0; pooled <- colSums(X); p_init <- pmax(pooled, 1) / sum(pmax(pooled, 1))
  t0_i <- log(pmax(alpha0_init * p_init, 1e-4))
  nll_i <- function(t, X) {
    alpha_vec <- exp(t); a0 <- sum(alpha_vec)
    -sum(vapply(1:G, function(g) {
      lgamma(a0) - lgamma(sum(X[g, ]) + a0) + sum(lgamma(X[g, ] + alpha_vec) - lgamma(alpha_vec))
    }, 0.0))
  }
  fit_i <- optim(t0_i, nll_i, score_param_i, method = "BFGS", control = list(maxit = 500, reltol = 1e-10), X = X)
  alpha_hat_i <- exp(fit_i$par); p_hat_i <- alpha_hat_i / sum(alpha_hat_i)
  
  ## MLE (ii)
  z0 <- log(p_init[1:(K-1)]) - log(p_init[K]); t0 <- rep(0, G)
  nll_ii <- function(par) {
    z <- par[1:(K-1)]; t <- par[K:(K-1+G)]; p <- softmax_from_z(z); a0 <- exp(t)
    -sum(vapply(1:G, function(g) dm_loglik_group(X[g,], p, a0[g]), 0.0))
  }
  grad_ii <- function(par) {
    z <- par[1:(K-1)]; t <- par[K:(K-1+G)]; p <- softmax_from_z(z); a0 <- exp(t)
    sc <- score_param_ii(X, p, a0); -c(sc$g_z, sc$g_t)
  }
  fit_ii <- try(optim(c(z0, t0), nll_ii, grad_ii, method = "BFGS", control = list(maxit = 500, reltol = 1e-10)), silent = TRUE)
  p_hat_ii <- if(inherits(fit_ii, "try-error")) rep(NA_real_, K) else softmax_from_z(fit_ii$par[1:(K-1)])
  
  ## MLE (iii)
  theta_fixed <- theta_true
  nll_iii <- function(z) { p <- softmax_from_z(z); -sum(vapply(1:G, function(g) dm_loglik_group(X[g,], p, theta_fixed * N_row[g]), 0.0)) }
  grad_iii <- function(z) { p <- softmax_from_z(z); sc <- score_param_iii(X, p, theta_fixed); -sc$g_z }
  fit_iii <- try(optim(z0, nll_iii, grad_iii, method = "BFGS", control = list(maxit = 500, reltol = 1e-10)), silent = TRUE)
  p_hat_iii <- if(inherits(fit_iii, "try-error")) rep(NA_real_, K) else softmax_from_z(fit_iii$par)
  
  ## MLE (iv)
  p_hat_iv <- colSums(X) / sum(colSums(X))
   
  ## MLE (v)
  nll_v <- function(par) {
    z <- par[1:(K-1)]; t <- par[K]; p <- softmax_from_z(z); theta <- exp(t)
    -sum(vapply(1:G, function(g) dm_loglik_group(X[g,], p, theta * N_row[g]), 0.0))
  }
  grad_v <- function(par) {
    z <- par[1:(K-1)]; t <- par[K]; p <- softmax_from_z(z); theta <- exp(t)
    sc <- score_param_iii(X, p, theta); -c(sc$g_z, sc$g_t)
  }
  fit_v <- try(optim(c(z0, 0), nll_v, grad_v, method = "BFGS", control = list(maxit = 500, reltol = 1e-10)), silent = TRUE)
  p_hat_v <- if(inherits(fit_v, "try-error")) rep(NA_real_, K) else softmax_from_z(fit_v$par[1:(K-1)])
  
  rmse_i = rmse(p_hat_i, p_true);   L1_norm_i = L1_norm(p_hat_i, p_true);   Linf_norm_i = Linf_norm(p_hat_i, p_true)
  rmse_ii = rmse(p_hat_ii, p_true); L1_norm_ii = L1_norm(p_hat_ii, p_true); Linf_norm_ii = Linf_norm(p_hat_ii, p_true)
  rmse_iii = rmse(p_hat_iii, p_true); L1_norm_iii = L1_norm(p_hat_iii, p_true); Linf_norm_iii = Linf_norm(p_hat_iii, p_true)
  rmse_iv = rmse(p_hat_iv, p_true); L1_norm_iv = L1_norm(p_hat_iv, p_true); Linf_norm_iv = Linf_norm(p_hat_iv, p_true)
  rmse_v = rmse(p_hat_v, p_true);   L1_norm_v = L1_norm(p_hat_v, p_true);   Linf_norm_v = Linf_norm(p_hat_v, p_true)

  return(c(
    p_true,
    rmse_i, L1_norm_i, Linf_norm_i, p_hat_i,
    rmse_ii, L1_norm_ii, Linf_norm_ii, p_hat_ii,
    rmse_iii, L1_norm_iii, Linf_norm_iii, p_hat_iii,
    rmse_iv, L1_norm_iv, Linf_norm_iv, p_hat_iv,
    rmse_v, L1_norm_v, Linf_norm_v, p_hat_v,
    X_out
  ))
}

# Run the simulation over all expanded combinations using parallel processing
message("Starting parallel simulation...")

## Safe future parallel setup for Windows batch runs
options(
  parallelly.makeNodePSOCK.timeout = 300,
  parallelly.makeNodePSOCK.connectTimeout = 300
)

n_workers <- as.integer(Sys.getenv("HAKE_N_WORKERS", unset = "4"))
n_workers <- max(1L, min(n_workers, future::availableCores()))

future::plan(future::multisession, workers = n_workers)

res_mat <- do.call(rbind, future_lapply(seq_len(total_runs), function(i) {
  c(mesh_id_vec[i], sim_id_vec[i], run_one(p_true = mesh_expanded[i, ], N_row = N_vec[i, ]))
}, future.seed = TRUE))

out <- as.data.frame(res_mat)
K_true <- K
true_names <- paste0("p", seq_len(K_true))
model_block_names <- c(
  "rmse_i",   "L1_norm_i",   "Linf_norm_i",   paste0("p_hat_i_",   seq_len(K_true)),
  "rmse_ii",  "L1_norm_ii",  "Linf_norm_ii",  paste0("p_hat_ii_",  seq_len(K_true)),
  "rmse_iii", "L1_norm_iii", "Linf_norm_iii", paste0("p_hat_iii_", seq_len(K_true)),
  "rmse_iv",  "L1_norm_iv",  "Linf_norm_iv",  paste0("p_hat_iv_",  seq_len(K_true)),
  "rmse_v",   "L1_norm_v",   "Linf_norm_v",   paste0("p_hat_v_",   seq_len(K_true))
)

X_block_names <- as.vector(sapply(seq_len(G), function(g) c(paste0("N", g), paste0("X", g, "_", seq_len(K)))))
colnames(out) <- c("mesh_id", "sim_id", true_names, model_block_names, X_block_names)

## ---- write CSV -------------------------------------------------------------
if (simplex_output_flag == 1) {
  csv_path <- paste0(out_prefix, ".csv")
  write.csv(out, csv_path, row.names = FALSE)
  message("Wrote: ", csv_path)
}

## ---- ternary plots (using averages for each mesh point to prevent overplotting) --------
if (K == 3) {
  suppressWarnings(suppressPackageStartupMessages(library(ggtern, quietly = TRUE, warn.conflicts = FALSE)))
  
  # Aggregate mean RMSE by mesh_id for cleaner ternary plots
  out_agg <- aggregate(out[, c("rmse_i", "rmse_ii", "rmse_iii", "rmse_iv", "rmse_v", "p1", "p2", "p3")], 
                       by = list(mesh_id = out$mesh_id), FUN = mean, na.rm = TRUE)
  
  rmse_min <- min(out_agg$rmse_i, out_agg$rmse_ii, out_agg$rmse_iii, out_agg$rmse_iv, out_agg$rmse_v, na.rm = TRUE)
  rmse_max <- max(out_agg$rmse_i, out_agg$rmse_ii, out_agg$rmse_iii, out_agg$rmse_iv, out_agg$rmse_v, na.rm = TRUE)
  
  plot_tern <- function(data, rmse_col, title_str) {
    ggtern::ggtern(data, aes(x = p1, y = p2, z = p3, colour = .data[[rmse_col]])) +
      geom_point(shape = 16, size = 2, alpha = 0.9) +
      labs(title = title_str, colour = "Mean RMSE") +
      scale_colour_viridis_c(limits = c(rmse_min, rmse_max), oob = scales::squish) +
      theme_bw()
  }
  
  p_i <- plot_tern(out_agg, "rmse_i", sprintf("(i) RMSE over simplex (K=3, G=%d, theta=%.1f)", G, theta_true))
  p_ii <- plot_tern(out_agg, "rmse_ii", sprintf("(ii) RMSE over simplex (K=3, G=%d, alpha0=beta)", G))
  p_iii <- plot_tern(out_agg, "rmse_iii", sprintf("(iii) RMSE over simplex (K=3, G=%d, alpha0=theta*N)", G))
  p_iv <- plot_tern(out_agg, "rmse_iv", sprintf("(iv) RMSE over simplex (K=3, G=%d, Multinomial)", G))
  p_v <- plot_tern(out_agg, "rmse_v", sprintf("(v) RMSE over simplex (K=3, G=%d, DML theta estimated)", G))
  
  print(p_i)
  print(p_ii)
  print(p_iii)
  print(p_iv)
  print(p_v)
  
  if (ggplot_output_flag == 1) {
    ggsave(paste0("rmse_i_", out_prefix, ".png"), p_i, width = 6, height = 5, dpi = 300)
    ggsave(paste0("rmse_ii_", out_prefix, ".png"), p_ii, width = 6, height = 5, dpi = 300)
    ggsave(paste0("rmse_iii_", out_prefix, ".png"), p_iii, width = 6, height = 5, dpi = 300)
    ggsave(paste0("rmse_iv_", out_prefix, ".png"), p_iv, width = 6, height = 5, dpi = 300)
    ggsave(paste0("rmse_v_", out_prefix, ".png"), p_v, width = 6, height = 5, dpi = 300)
  }
  
  # Aggregate mean L1 norm by mesh_id for cleaner ternary plots
  out_agg_L1 <- aggregate(out[, c("L1_norm_i", "L1_norm_ii", "L1_norm_iii", "L1_norm_iv", "L1_norm_v", "p1", "p2", "p3")],
                          by = list(mesh_id = out$mesh_id), FUN = mean, na.rm = TRUE)

  L1_norm_min <- min(out_agg_L1$L1_norm_i, out_agg_L1$L1_norm_ii, out_agg_L1$L1_norm_iii, out_agg_L1$L1_norm_iv, out_agg_L1$L1_norm_v, na.rm = TRUE)
  L1_norm_max <- max(out_agg_L1$L1_norm_i, out_agg_L1$L1_norm_ii, out_agg_L1$L1_norm_iii, out_agg_L1$L1_norm_iv, out_agg_L1$L1_norm_v, na.rm = TRUE)

  plot_tern_L1 <- function(data, L1_norm_col, title_str) {
    ggtern::ggtern(data, aes(x = p1, y = p2, z = p3, colour = .data[[L1_norm_col]])) +
      geom_point(shape = 16, size = 2, alpha = 0.9) +
      labs(title = title_str, colour = "Mean L1 norm") +
      scale_colour_viridis_c(limits = c(L1_norm_min, L1_norm_max), oob = scales::squish) +
      theme_bw()
  }

  p_L1_i <- plot_tern_L1(out_agg_L1, "L1_norm_i", sprintf("(i) L1 norm over simplex (K=3, G=%d, theta=%.1f)", G, theta_true))
  p_L1_ii <- plot_tern_L1(out_agg_L1, "L1_norm_ii", sprintf("(ii) L1 norm over simplex (K=3, G=%d, alpha0=beta)", G))
  p_L1_iii <- plot_tern_L1(out_agg_L1, "L1_norm_iii", sprintf("(iii) L1 norm over simplex (K=3, G=%d, alpha0=theta*N)", G))
  p_L1_iv <- plot_tern_L1(out_agg_L1, "L1_norm_iv", sprintf("(iv) L1 norm over simplex (K=3, G=%d, Multinomial)", G))
  p_L1_v <- plot_tern_L1(out_agg_L1, "L1_norm_v", sprintf("(v) L1 norm over simplex (K=3, G=%d, DML theta estimated)", G))

  print(p_L1_i)
  print(p_L1_ii)
  print(p_L1_iii)
  print(p_L1_iv)
  print(p_L1_v)

  if (ggplot_output_flag == 1) {
    ggsave(paste0("L1_norm_i_", out_prefix, ".png"), p_L1_i, width = 6, height = 5, dpi = 300)
    ggsave(paste0("L1_norm_ii_", out_prefix, ".png"), p_L1_ii, width = 6, height = 5, dpi = 300)
    ggsave(paste0("L1_norm_iii_", out_prefix, ".png"), p_L1_iii, width = 6, height = 5, dpi = 300)
    ggsave(paste0("L1_norm_iv_", out_prefix, ".png"), p_L1_iv, width = 6, height = 5, dpi = 300)
    ggsave(paste0("L1_norm_v_", out_prefix, ".png"), p_L1_v, width = 6, height = 5, dpi = 300)
  }
  
  # Aggregate mean Linf norm by mesh_id for cleaner ternary plots
  out_agg_Linf <- aggregate(out[, c("Linf_norm_i", "Linf_norm_ii", "Linf_norm_iii", "Linf_norm_iv", "Linf_norm_v", "p1", "p2", "p3")],
                            by = list(mesh_id = out$mesh_id), FUN = mean, na.rm = TRUE)

  Linf_norm_min <- min(out_agg_Linf$Linf_norm_i, out_agg_Linf$Linf_norm_ii, out_agg_Linf$Linf_norm_iii, out_agg_Linf$Linf_norm_iv, out_agg_Linf$Linf_norm_v, na.rm = TRUE)
  Linf_norm_max <- max(out_agg_Linf$Linf_norm_i, out_agg_Linf$Linf_norm_ii, out_agg_Linf$Linf_norm_iii, out_agg_Linf$Linf_norm_iv, out_agg_Linf$Linf_norm_v, na.rm = TRUE)

  plot_tern_Linf <- function(data, Linf_norm_col, title_str) {
    ggtern::ggtern(data, aes(x = p1, y = p2, z = p3, colour = .data[[Linf_norm_col]])) +
      geom_point(shape = 16, size = 2, alpha = 0.9) +
      labs(title = title_str, colour = "Mean Linf norm") +
      scale_colour_viridis_c(limits = c(Linf_norm_min, Linf_norm_max), oob = scales::squish) +
      theme_bw()
  }

  p_Linf_i <- plot_tern_Linf(out_agg_Linf, "Linf_norm_i", sprintf("(i) Linf norm over simplex (K=3, G=%d, theta=%.1f)", G, theta_true))
  p_Linf_ii <- plot_tern_Linf(out_agg_Linf, "Linf_norm_ii", sprintf("(ii) Linf norm over simplex (K=3, G=%d, alpha0=beta)", G))
  p_Linf_iii <- plot_tern_Linf(out_agg_Linf, "Linf_norm_iii", sprintf("(iii) Linf norm over simplex (K=3, G=%d, alpha0=theta*N)", G))
  p_Linf_iv <- plot_tern_Linf(out_agg_Linf, "Linf_norm_iv", sprintf("(iv) Linf norm over simplex (K=3, G=%d, Multinomial)", G))
  p_Linf_v <- plot_tern_Linf(out_agg_Linf, "Linf_norm_v", sprintf("(v) Linf norm over simplex (K=3, G=%d, DML theta estimated)", G))

  print(p_Linf_i)
  print(p_Linf_ii)
  print(p_Linf_iii)
  print(p_Linf_iv)
  print(p_Linf_v)

  if (ggplot_output_flag == 1) {
    ggsave(paste0("Linf_norm_i_", out_prefix, ".png"), p_Linf_i, width = 6, height = 5, dpi = 300)
    ggsave(paste0("Linf_norm_ii_", out_prefix, ".png"), p_Linf_ii, width = 6, height = 5, dpi = 300)
    ggsave(paste0("Linf_norm_iii_", out_prefix, ".png"), p_Linf_iii, width = 6, height = 5, dpi = 300)
    ggsave(paste0("Linf_norm_iv_", out_prefix, ".png"), p_Linf_iv, width = 6, height = 5, dpi = 300)
    ggsave(paste0("Linf_norm_v_", out_prefix, ".png"), p_Linf_v, width = 6, height = 5, dpi = 300)
  }
}

print(params)

## Statistics
n_obs <- nrow(out)
df <- data.frame(
  id     = rep(seq_len(n_obs), times = 5), # id now correctly identifies the specific simulation dataset run
  method = factor(rep(c("i","ii","iii","iv","v"), each = n_obs), levels = c("i","ii","iii","iv","v")),
  rmse   = c(out$rmse_i, out$rmse_ii, out$rmse_iii, out$rmse_iv, out$rmse_v)
)

ft <- friedman.test(rmse ~ method | id, data = df)
print(ft)

k <- nlevels(df$method)
W <- as.numeric(ft$statistic) / (n_obs * (k - 1))
cat(sprintf("Kendall's W ~= %.3f\n", W))

pw <- pairwise.wilcox.test(df$rmse, df$method, paired = TRUE, p.adjust.method = "holm", exact = FALSE)
print(pw)

wide <- reshape(df, idvar = "id", timevar = "method", direction = "wide")
methods <- c("i","ii","iii","iv","v")
pair_list <- combn(methods, 2, simplify = FALSE)

wilcox_ci <- lapply(pair_list, function(pr) {
  a <- wide[[paste0("rmse.", pr[1])]]; b <- wide[[paste0("rmse.", pr[2])]]
  wt <- wilcox.test(b, a, paired = TRUE, conf.int = TRUE, exact = FALSE)
  list(comp = paste0(pr[2], " - ", pr[1]), test = wt)
})
for (obj in wilcox_ci) { cat("\nPaired Wilcoxon (", obj$comp, "):\n", sep = ""); print(obj$test) }

pct <- wide
for (ab in pair_list) {
  nm <- paste0("pct_", ab[2], "_", ab[1])
  pct[[nm]] <- 100 * (pct[[paste0("rmse.", ab[2])]] - pct[[paste0("rmse.", ab[1])]]) / pct[[paste0("rmse.", ab[1])]]
}
pct_cols <- paste0("pct_", sapply(pair_list, `[`, 2), "_", sapply(pair_list, `[`, 1))
tmp <- sapply(pct[pct_cols], \(x) c(median = median(x, na.rm=TRUE), IQR = IQR(x, na.rm=TRUE), q25 = quantile(x, .25, na.rm=TRUE), q75 = quantile(x, .75, na.rm=TRUE)))
cat("Percent change in RMSE for method X relative to method Y", "\n\n")
print(tmp)

m <- lmer(rmse ~ method + (1 | id), data = df)
if(any(grepl("singular", m@optinfo$conv$lme4$messages))) {
  warning("Linear Mixed Model fit is singular. Variance estimates may be unreliable.")
}
print(anova(m))
print(emmeans(m, pairwise ~ method, adjust = "holm"))

df_L1 <- data.frame(
  id      = rep(seq_len(n_obs), times = 5),
  method  = factor(rep(c("i","ii","iii","iv","v"), each = n_obs), levels = c("i","ii","iii","iv","v")),
  L1_norm = c(out$L1_norm_i, out$L1_norm_ii, out$L1_norm_iii, out$L1_norm_iv, out$L1_norm_v)
)

ft_L1 <- friedman.test(L1_norm ~ method | id, data = df_L1)
print(ft_L1)

k_L1 <- nlevels(df_L1$method)
W_L1 <- as.numeric(ft_L1$statistic) / (n_obs * (k_L1 - 1))
cat(sprintf("Kendall's W for L1 norm ~= %.3f\n", W_L1))

pw_L1 <- pairwise.wilcox.test(df_L1$L1_norm, df_L1$method, paired = TRUE, p.adjust.method = "holm", exact = FALSE)
print(pw_L1)

wide_L1 <- reshape(df_L1, idvar = "id", timevar = "method", direction = "wide")
methods_L1 <- c("i","ii","iii","iv","v")
pair_list_L1 <- combn(methods_L1, 2, simplify = FALSE)

wilcox_ci_L1 <- lapply(pair_list_L1, function(pr) {
  a <- wide_L1[[paste0("L1_norm.", pr[1])]]; b <- wide_L1[[paste0("L1_norm.", pr[2])]]
  wt <- wilcox.test(b, a, paired = TRUE, conf.int = TRUE, exact = FALSE)
  list(comp = paste0(pr[2], " - ", pr[1]), test = wt)
})
for (obj in wilcox_ci_L1) { cat("\nPaired Wilcoxon for L1 norm (", obj$comp, "):\n", sep = ""); print(obj$test) }

pct_L1 <- wide_L1
for (ab in pair_list_L1) {
  nm <- paste0("pct_L1_norm_", ab[2], "_", ab[1])
  pct_L1[[nm]] <- 100 * (pct_L1[[paste0("L1_norm.", ab[2])]] - pct_L1[[paste0("L1_norm.", ab[1])]]) / pct_L1[[paste0("L1_norm.", ab[1])]]
}
pct_cols_L1 <- paste0("pct_L1_norm_", sapply(pair_list_L1, `[`, 2), "_", sapply(pair_list_L1, `[`, 1))
tmp_L1 <- sapply(pct_L1[pct_cols_L1], \(x) c(median = median(x, na.rm=TRUE), IQR = IQR(x, na.rm=TRUE), q25 = quantile(x, .25, na.rm=TRUE), q75 = quantile(x, .75, na.rm=TRUE)))
cat("Percent change in L1 norm for method X relative to method Y", "\n\n")
print(tmp_L1)

m_L1 <- lmer(L1_norm ~ method + (1 | id), data = df_L1)
if(any(grepl("singular", m_L1@optinfo$conv$lme4$messages))) {
  warning("Linear Mixed Model fit for L1 norm is singular. Variance estimates may be unreliable.")
}
print(anova(m_L1))
print(emmeans(m_L1, pairwise ~ method, adjust = "holm"))

df_Linf <- data.frame(
  id        = rep(seq_len(n_obs), times = 5),
  method    = factor(rep(c("i","ii","iii","iv","v"), each = n_obs), levels = c("i","ii","iii","iv","v")),
  Linf_norm = c(out$Linf_norm_i, out$Linf_norm_ii, out$Linf_norm_iii, out$Linf_norm_iv, out$Linf_norm_v)
)

ft_Linf <- friedman.test(Linf_norm ~ method | id, data = df_Linf)
print(ft_Linf)

k_Linf <- nlevels(df_Linf$method)
W_Linf <- as.numeric(ft_Linf$statistic) / (n_obs * (k_Linf - 1))
cat(sprintf("Kendall's W for Linf norm ~= %.3f\n", W_Linf))

pw_Linf <- pairwise.wilcox.test(df_Linf$Linf_norm, df_Linf$method, paired = TRUE, p.adjust.method = "holm", exact = FALSE)
print(pw_Linf)

wide_Linf <- reshape(df_Linf, idvar = "id", timevar = "method", direction = "wide")
methods_Linf <- c("i","ii","iii","iv","v")
pair_list_Linf <- combn(methods_Linf, 2, simplify = FALSE)

wilcox_ci_Linf <- lapply(pair_list_Linf, function(pr) {
  a <- wide_Linf[[paste0("Linf_norm.", pr[1])]]; b <- wide_Linf[[paste0("Linf_norm.", pr[2])]]
  wt <- wilcox.test(b, a, paired = TRUE, conf.int = TRUE, exact = FALSE)
  list(comp = paste0(pr[2], " - ", pr[1]), test = wt)
})
for (obj in wilcox_ci_Linf) { cat("\nPaired Wilcoxon for Linf norm (", obj$comp, "):\n", sep = ""); print(obj$test) }

pct_Linf <- wide_Linf
for (ab in pair_list_Linf) {
  nm <- paste0("pct_Linf_norm_", ab[2], "_", ab[1])
  pct_Linf[[nm]] <- 100 * (pct_Linf[[paste0("Linf_norm.", ab[2])]] - pct_Linf[[paste0("Linf_norm.", ab[1])]]) / pct_Linf[[paste0("Linf_norm.", ab[1])]]
}
pct_cols_Linf <- paste0("pct_Linf_norm_", sapply(pair_list_Linf, `[`, 2), "_", sapply(pair_list_Linf, `[`, 1))
tmp_Linf <- sapply(pct_Linf[pct_cols_Linf], \(x) c(median = median(x, na.rm=TRUE), IQR = IQR(x, na.rm=TRUE), q25 = quantile(x, .25, na.rm=TRUE), q75 = quantile(x, .75, na.rm=TRUE)))
cat("Percent change in Linf norm for method X relative to method Y", "\n\n")
print(tmp_Linf)

m_Linf <- lmer(Linf_norm ~ method + (1 | id), data = df_Linf)
if(any(grepl("singular", m_Linf@optinfo$conv$lme4$messages))) {
  warning("Linear Mixed Model fit for Linf norm is singular. Variance estimates may be unreliable.")
}
print(anova(m_Linf))
print(emmeans(m_Linf, pairwise ~ method, adjust = "holm"))

standardize_emmeans_ci_names <- function(x) {
  nms <- names(x)
  
  if ("lower.CL" %in% nms && "upper.CL" %in% nms) {
    x$ci_lower <- x$lower.CL
    x$ci_upper <- x$upper.CL
    return(x)
  }
  
  if ("asymp.LCL" %in% nms && "asymp.UCL" %in% nms) {
    x$ci_lower <- x$asymp.LCL
    x$ci_upper <- x$asymp.UCL
    return(x)
  }
  
  stop(
    "Could not find confidence interval columns in emmeans output. ",
    "Available columns are: ", paste(nms, collapse = ", ")
  )
}

## ---- Extract and Visualize Statistical Summaries ---------------------------
# 1. Extract the EMMeans into a clean data frame
emm_results <- emmeans(m, pairwise ~ method, adjust = "holm")

emm_summary_df <- standardize_emmeans_ci_names(
  as.data.frame(summary(emm_results$emmeans, infer = TRUE))
)

# Save the summary table to a CSV
write.csv(emm_summary_df, paste0(out_prefix, "_emmeans_summary.csv"), row.names = FALSE)
message("Wrote: ", paste0(out_prefix, "_emmeans_summary.csv"))

# 2. Plot 1: Boxplot of the raw RMSE data
p_box <- ggplot(df, aes(x = method, y = rmse, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.4) +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Distribution of Raw RMSE by Estimation Method",
    subtitle = sprintf("Based on %d simulations (K=%d, G=%d)", n_obs, K, G),
    x = "Estimation Method",
    y = "Root Mean Square Error (RMSE)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# 3. Plot 2: Point-Range plot of the Estimated Marginal Means
p_emm <- ggplot(emm_summary_df, aes(x = method, y = emmean, color = method)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 1) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Estimated Marginal Mean RMSE (with 95% CIs)",
    subtitle = "Parametric estimates from Linear Mixed-Effects Model",
    x = "Estimation Method",
    y = "Estimated Mean RMSE"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Print and save the plots
print(p_box)
print(p_emm)
if (ggplot_output_flag == 1) {
  ggsave(paste0("rmse_boxplot_", out_prefix, ".png"), p_box, width = 6, height = 5, dpi = 300)
  ggsave(paste0("rmse_emmeans_", out_prefix, ".png"), p_emm, width = 6, height = 5, dpi = 300)
}
# 4. Extract the EMMeans for L1 norm into a clean data frame
emm_results_L1 <- emmeans(m_L1, pairwise ~ method, adjust = "holm")

emm_summary_df_L1 <- standardize_emmeans_ci_names(
  as.data.frame(summary(emm_results_L1$emmeans, infer = TRUE))
)

# Save the L1 norm summary table to a CSV
write.csv(emm_summary_df_L1, paste0(out_prefix, "_L1_norm_emmeans_summary.csv"), row.names = FALSE)
message("Wrote: ", paste0(out_prefix, "_L1_norm_emmeans_summary.csv"))

# 5. Plot 3: Boxplot of the raw L1 norm data
p_box_L1 <- ggplot(df_L1, aes(x = method, y = L1_norm, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.4) +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Distribution of Raw L1 Norm by Estimation Method",
    subtitle = sprintf("Based on %d simulations (K=%d, G=%d)", n_obs, K, G),
    x = "Estimation Method",
    y = "L1 Norm"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# 6. Plot 4: Point-Range plot of the Estimated Marginal Means for L1 norm
p_emm_L1 <- ggplot(emm_summary_df_L1, aes(x = method, y = emmean, color = method)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 1) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Estimated Marginal Mean L1 Norm (with 95% CIs)",
    subtitle = "Parametric estimates from Linear Mixed-Effects Model",
    x = "Estimation Method",
    y = "Estimated Mean L1 Norm"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Print and save the L1 norm plots
print(p_box_L1)
print(p_emm_L1)
if (ggplot_output_flag == 1) {
  ggsave(paste0("L1_norm_boxplot_", out_prefix, ".png"), p_box_L1, width = 6, height = 5, dpi = 300)
  ggsave(paste0("L1_norm_emmeans_", out_prefix, ".png"), p_emm_L1, width = 6, height = 5, dpi = 300)
}

# 7. Extract the EMMeans for Linf norm into a clean data frame
emm_results_Linf <- emmeans(m_Linf, pairwise ~ method, adjust = "holm")

emm_summary_df_Linf <- standardize_emmeans_ci_names(
  as.data.frame(summary(emm_results_Linf$emmeans, infer = TRUE))
)

# Save the Linf norm summary table to a CSV
write.csv(emm_summary_df_Linf, paste0(out_prefix, "_Linf_norm_emmeans_summary.csv"), row.names = FALSE)
message("Wrote: ", paste0(out_prefix, "_Linf_norm_emmeans_summary.csv"))

# 8. Plot 5: Boxplot of the raw Linf norm data
p_box_Linf <- ggplot(df_Linf, aes(x = method, y = Linf_norm, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.4) +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Distribution of Raw Linf Norm by Estimation Method",
    subtitle = sprintf("Based on %d simulations (K=%d, G=%d)", n_obs, K, G),
    x = "Estimation Method",
    y = "Linf Norm"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# 9. Plot 6: Point-Range plot of the Estimated Marginal Means for Linf norm
p_emm_Linf <- ggplot(emm_summary_df_Linf, aes(x = method, y = emmean, color = method)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 1) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Estimated Marginal Mean Linf Norm (with 95% CIs)",
    subtitle = "Parametric estimates from Linear Mixed-Effects Model",
    x = "Estimation Method",
    y = "Estimated Mean Linf Norm"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Print and save the Linf norm plots
print(p_box_Linf)
print(p_emm_Linf)
if (ggplot_output_flag == 1) {
  ggsave(paste0("Linf_norm_boxplot_", out_prefix, ".png"), p_box_Linf, width = 6, height = 5, dpi = 300)
  ggsave(paste0("Linf_norm_emmeans_", out_prefix, ".png"), p_emm_Linf, width = 6, height = 5, dpi = 300)
}

gc() # Clear memory
cat("\n")

end_time <- Sys.time()
cat("Run Time End:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n\n")
run_time <- end_time - start_time
cat("Total Elapsed Time:", round(run_time, 2), attr(run_time, "units"), "\n\n")

})