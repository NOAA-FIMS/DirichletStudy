## file = hake_selftest.R
## Purpose: Fast self-test comparing Multinomial (MN) vs Dirichlet-multinomial (DM)
## Design (fast high-value run):
##  - Fix one operating-model “truth” using values in hake_selftest.inp
##  - GEN ∈ {MN, DMstrong}
##  - Nsamp ∈ {25, 125, 250}
##  - FIT ∈ {MN, DM}
##  - R = 200 reps per (GEN, Nsamp)
##
## Run in RStudio:
##   source("hake_selftest.R", echo = FALSE, print.eval = TRUE)
##
## Outputs (prefix set in .inp):
##  - CSV (raw):    <out_prefix>_raw.csv
##  - CSV (summary):<out_prefix>_summary.csv
##  - Plot:         <out_prefix>_rmse_boxplot.png
##  - List:         <out_prefix>.lst

suppressWarnings(suppressPackageStartupMessages({
  library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
}))

# ----------------------------- utilities ----------------------------------

softmax_from_z <- function(z) {
  v <- c(z, 0)
  ev <- exp(v - max(v))
  ev / sum(ev)
}

rmse <- function(a, b) sqrt(mean((a - b)^2))

dm_loglik_group <- function(x, p, a0) {
  # Dirichlet-multinomial log-likelihood for one group with total concentration a0 and mean p
  N <- sum(x)
  lgamma(a0) - lgamma(N + a0) + sum(lgamma(x + a0 * p) - lgamma(a0 * p))
}

# Score for DML: p shared; alpha0_g = theta * N_g  (opt over z and t=log theta)
score_param_dml <- function(X, p, theta, N_row) {
  G <- nrow(X); K <- length(p)
  g_p <- numeric(K)
  g_theta <- 0
  for (g in 1:G) {
    xg <- X[g, ]
    Ng <- sum(xg)
    a0 <- theta * N_row[g]
    g_p <- g_p + a0 * (digamma(xg + a0 * p) - digamma(a0 * p))
    g_theta <- g_theta + N_row[g] * (
      digamma(a0) - digamma(Ng + a0) + sum(p * (digamma(xg + a0 * p) - digamma(a0 * p)))
    )
  }
  J <- diag(p) - tcrossprod(p)
  g_z <- as.vector(crossprod(J[, 1:(K-1), drop=FALSE], g_p))
  g_t <- theta * g_theta
  list(g_z = g_z, g_t = g_t)
}

fit_mn <- function(X) {
  pooled <- colSums(X)
  p_hat <- pooled / sum(pooled)
  list(p_hat = p_hat, theta_hat = NA_real_, conv = TRUE, nll = NA_real_)
}

fit_dml <- function(X, N_row, z0, t0 = 0, maxit = 500) {
  # MLE of (p, theta) under alpha0_g = theta * N_g
  K <- ncol(X)

  nll <- function(par) {
    z <- par[1:(K-1)]
    t <- par[K]
    p <- softmax_from_z(z)
    theta <- exp(t)
    -sum(vapply(1:nrow(X), function(g) dm_loglik_group(X[g,], p, theta * N_row[g]), 0.0))
  }

  grad <- function(par) {
    z <- par[1:(K-1)]
    t <- par[K]
    p <- softmax_from_z(z)
    theta <- exp(t)
    sc <- score_param_dml(X, p, theta, N_row)
    -c(sc$g_z, sc$g_t)
  }

  fit <- try(optim(c(z0, t0), nll, grad, method = "BFGS",
                   control = list(maxit = maxit, reltol = 1e-10)),
             silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(list(p_hat = rep(NA_real_, K), theta_hat = NA_real_, conv = FALSE, nll = NA_real_))
  }
  conv <- isTRUE(fit$convergence == 0)
  zhat <- fit$par[1:(K-1)]
  p_hat <- softmax_from_z(zhat)
  theta_hat <- exp(fit$par[K])
  list(p_hat = p_hat, theta_hat = theta_hat, conv = conv, nll = fit$value)
}

simulate_counts <- function(gen, p_true, N_row, theta_true, sigma, od_mult) {
  # Returns X: G x K integer counts
  G <- length(N_row)
  K <- length(p_true)
  X <- matrix(0L, nrow = G, ncol = K)

  if (gen == "MN") {
    for (g in 1:G) {
      X[g, ] <- as.vector(rmultinom(1, size = N_row[g], prob = p_true))
    }
    return(X)
  }

  if (gen == "DMstrong") {
    for (g in 1:G) {
      a0g <- (theta_true * od_mult) * N_row[g]  # smaller => more overdispersion
      z <- log(p_true) + rnorm(K, 0, sigma)     # induces between-sample variation in p_g
      p_g <- exp(z); p_g <- p_g / sum(p_g)
      alpha_vec <- a0g * p_g
      # DM draw via normalized gammas (Dirichlet), then multinomial
      phi <- rgamma(K, shape = alpha_vec, rate = 1); phi <- phi / sum(phi)
      X[g, ] <- as.vector(rmultinom(1, size = N_row[g], prob = phi))
    }
    return(X)
  }

  stop("Unknown gen. Use GEN in {MN, DMstrong}.")
}

# ----------------------------- read input ----------------------------------

args <- commandArgs(trailingOnly = TRUE)
inp_path <- if (length(args) >= 1) args[1] else "hake_selftest.inp"
if (!file.exists(inp_path)) stop(sprintf("Input file '%s' not found.", inp_path))

tbl <- read.table(
  file = inp_path,
  sep = "=",
  comment.char = "#",
  strip.white = TRUE,
  blank.lines.skip = TRUE,
  col.names = c("key", "value"),
  colClasses = c("character", "character")
)
kv <- setNames(trimws(tbl$value), trimws(tbl$key))

get_num <- function(key, default = NA_real_) {
  if (!key %in% names(kv)) return(default)
  as.numeric(kv[[key]])
}
get_int <- function(key, default = NA_integer_) {
  if (!key %in% names(kv)) return(default)
  as.integer(kv[[key]])
}
get_chr <- function(key, default = NA_character_) {
  if (!key %in% names(kv)) return(default)
  as.character(kv[[key]])
}
get_vec_num <- function(key, K_expected = NULL, default = NULL) {
  if (!key %in% names(kv)) return(default)
  raw <- kv[[key]]
  parts <- unlist(strsplit(raw, "[,[:space:]]+"))
  parts <- parts[nzchar(parts)]
  v <- as.numeric(parts)
  if (any(is.na(v))) stop(sprintf("%s must be numeric.", key))
  if (!is.null(K_expected) && length(v) != K_expected) {
    stop(sprintf("%s must have length %d.", key, K_expected))
  }
  v
}

K          <- get_int("K")
G          <- get_int("G")
theta_true <- get_num("theta_true")
sigma      <- get_num("sigma")
od_mult    <- get_num("od_mult")
R          <- get_int("R", 200L)
random_seed<- get_int("random.seed", 12345L)

out_prefix <- get_chr("out_prefix", "hake_selftest")
p_true     <- get_vec_num("p_true", K_expected = K, default = rep(1/K, K))
nsamp_levels <- get_vec_num("nsamp_levels", default = c(25, 125, 250))
nsamp_levels <- as.integer(nsamp_levels)

# validation
stopifnot(is.finite(K), K >= 2, is.finite(G), G >= 1)
stopifnot(all(p_true > 0), abs(sum(p_true) - 1) < 1e-8)
stopifnot(all(nsamp_levels >= 1))
stopifnot(is.finite(theta_true) && theta_true > 0)
stopifnot(is.finite(sigma) && sigma > 0)
stopifnot(is.finite(od_mult) && od_mult > 0)
stopifnot(R >= 1)

set.seed(random_seed)

params <- list(
  K = K, G = G, p_true = p_true,
  theta_true = theta_true, sigma = sigma, od_mult = od_mult,
  GEN = c("MN", "DMstrong"), Nsamp = nsamp_levels, R = R,
  random.seed = random_seed, out_prefix = out_prefix
)

# ----------------------------- logging -------------------------------------

logfile <- paste0(out_prefix, ".lst")
sink(logfile, split = TRUE, type = "output")

cat("Logging to: ", logfile, "\n", sep = "")
cat("Parameters:\n")
print(params)

# ----------------------------- scenario grid --------------------------------

scen <- expand.grid(
  GEN = c("MN", "DMstrong"),
  Nsamp = nsamp_levels,
  rep = seq_len(R),
  stringsAsFactors = FALSE
)

# preallocate results list (faster than rbind)
res <- vector("list", nrow(scen) * 2L)  # *2 fits per replicate
idx <- 0L

# starting values for p (in z-space)
z0 <- log(p_true[1:(K-1)]) - log(p_true[K])

# ----------------------------- main loop ------------------------------------

for (ii in seq_len(nrow(scen))) {
  
  print(ii)
  
  gen  <- scen$GEN[ii]
  Ns   <- scen$Nsamp[ii]
  rr   <- scen$rep[ii]

  # fixed N_row for this scenario (same for each group)
  N_row <- rep(as.integer(Ns), G)

  # simulate one replicate dataset
  X <- simulate_counts(gen, p_true, N_row, theta_true, sigma, od_mult)

  # FIT = MN
  mn <- fit_mn(X)
  idx <- idx + 1L
  res[[idx]] <- data.frame(
    GEN = gen, Nsamp = Ns, rep = rr, FIT = "MN",
    conv = mn$conv,
    rmse_p = rmse(mn$p_hat, p_true),
    theta_hat = mn$theta_hat,
    stringsAsFactors = FALSE
  )
  for (k in 1:K) res[[idx]][[paste0("p", k, "_hat")]] <- mn$p_hat[k]

  # FIT = DM (DML with theta estimated)
  dm <- fit_dml(X, N_row, z0 = z0, t0 = log(theta_true), maxit = 500)
  idx <- idx + 1L
  res[[idx]] <- data.frame(
    GEN = gen, Nsamp = Ns, rep = rr, FIT = "DM",
    conv = dm$conv,
    rmse_p = rmse(dm$p_hat, p_true),
    theta_hat = dm$theta_hat,
    stringsAsFactors = FALSE
  )
  for (k in 1:K) res[[idx]][[paste0("p", k, "_hat")]] <- dm$p_hat[k]
}

raw <- do.call(rbind, res)

# ----------------------------- summaries ------------------------------------

# convergence rate + RMSE summaries
summ <- aggregate(
  cbind(conv = as.numeric(raw$conv), rmse_p = raw$rmse_p, theta_hat = raw$theta_hat) ~ GEN + Nsamp + FIT,
  data = raw,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
)

unpack <- function(matcol, prefix) {
  data.frame(
    setNames(list(matcol[, "mean"]), paste0(prefix, "_mean")),
    setNames(list(matcol[, "sd"]),   paste0(prefix, "_sd"))
  )
}
summary_df <- cbind(
  summ[, c("GEN","Nsamp","FIT")],
  unpack(summ$conv, "conv_rate"),
  unpack(summ$rmse_p, "rmse_p"),
  unpack(summ$theta_hat, "theta_hat")
)

# parameter bias (per component)
for (k in 1:K) {
  nm <- paste0("p", k, "_hat")
  bias_df <- aggregate(
    raw[[nm]] ~ GEN + Nsamp + FIT,
    data = raw,
    FUN = function(x) mean(x - p_true[k], na.rm = TRUE)
  )
  names(bias_df)[4] <- paste0("bias_p", k)
  summary_df <- merge(summary_df, bias_df, by = c("GEN","Nsamp","FIT"), all.x = TRUE, sort = FALSE)
}

# ----------------------------- write outputs --------------------------------

raw_path <- paste0(out_prefix, "_raw.csv")
sum_path <- paste0(out_prefix, "_summary.csv")
write.csv(raw, raw_path, row.names = FALSE)
write.csv(summary_df, sum_path, row.names = FALSE)
cat("Wrote: ", raw_path, "\n", sep = "")
cat("Wrote: ", sum_path, "\n", sep = "")

# Plot: RMSE by FIT, faceted by GEN and Nsamp
p <- ggplot(raw, aes(x = FIT, y = rmse_p)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_grid(GEN ~ Nsamp, scales = "free_y") +
  labs(title = "Self-test: RMSE of p-hat (MN vs DM fits)",
       subtitle = sprintf("Truth p=(%s), G=%d, R=%d", paste(round(p_true, 3), collapse=", "), G, R),
       x = "Fit likelihood", y = "RMSE(p_hat, p_true)") +
  theme_bw()

plot_path <- paste0(out_prefix, "_rmse_boxplot.png")
ggsave(plot_path, p, width = 9, height = 5.5, dpi = 300)
cat("Wrote: ", plot_path, "\n", sep = "")

# Print a compact summary to console / list file
cat("\n--- Summary (means) ---\n")
print(summary_df)

# Show plot
print(p)

on.exit(sink())
