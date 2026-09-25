## file = validate_hake_selftest.R
## Purpose: Validate DM optimization + diagnose large theta estimates via
##   (1) adding a "pure DM" generator (DMpure) that matches the fitted DML form,
##   (2) bounding theta during DM fitting,
##   (3) recording theta profiles for a subset of replicates,
##   (4) reporting how often theta hits the upper bound (suggesting MN-like data).
##
## Run in RStudio:
##   source("validate_hake_selftest.R", echo = FALSE, print.eval = TRUE)
##
## Input file (default): validate_hake_selftest.inp
## Backward compatible with hake_selftest.R inputs. Optional new keys:
##   theta_lower = 1e-6
##   theta_upper = 1e4
##   profile_reps = 5
##   profile_grid_n = 40
##   include_DMpure = 1
##   include_DMstrong = 1
##
## Outputs (prefix set in .inp):
##  - CSV (raw):            <out_prefix>_raw.csv
##  - CSV (summary):        <out_prefix>_summary.csv
##  - CSV (theta profiles): <out_prefix>_theta_profile.csv
##  - Plot:                 <out_prefix>_rmse_boxplot.png
##  - List:                 <out_prefix>.lst

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
  # MN has no theta; keep theta_hat as NA
  list(p_hat = p_hat, theta_hat = NA_real_, conv = TRUE, nll = NA_real_, hit_upper = NA)
}

fit_dml_bounded <- function(X, N_row, z0, t0, theta_lower, theta_upper, maxit = 500) {
  # MLE of (p, theta) under alpha0_g = theta * N_g, with bounds on theta
  K <- ncol(X)
  t_lo <- log(theta_lower)
  t_hi <- log(theta_upper)

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

  lower <- c(rep(-Inf, K-1), t_lo)
  upper <- c(rep( Inf, K-1), t_hi)

  fit <- try(optim(c(z0, t0), nll, grad,
                   method = "L-BFGS-B",
                   lower = lower, upper = upper,
                   control = list(maxit = maxit, factr = 1e7)),
             silent = TRUE)

  if (inherits(fit, "try-error")) {
    return(list(p_hat = rep(NA_real_, K), theta_hat = NA_real_, conv = FALSE,
                nll = NA_real_, hit_upper = NA))
  }

  conv <- isTRUE(fit$convergence == 0)
  zhat <- fit$par[1:(K-1)]
  p_hat <- softmax_from_z(zhat)
  theta_hat <- exp(fit$par[K])

  # bound-hit flag (tolerance in log space)
  hit_upper <- isTRUE(abs(fit$par[K] - t_hi) <= 1e-6)

  list(p_hat = p_hat, theta_hat = theta_hat, conv = conv, nll = fit$value, hit_upper = hit_upper)
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

  if (gen == "DMpure") {
    # Pure DM generator that matches the fitted mean structure (shared p_true)
    # Uses a0g = (theta_true * od_mult) * N_g, consistent with the earlier generator scaling
    for (g in 1:G) {
      a0g <- (theta_true * od_mult) * N_row[g]
      alpha_vec <- a0g * p_true
      phi <- rgamma(K, shape = alpha_vec, rate = 1); phi <- phi / sum(phi)
      X[g, ] <- as.vector(rmultinom(1, size = N_row[g], prob = phi))
    }
    return(X)
  }

  if (gen == "DMstrong") {
    # Mis-specified "strong DM": adds between-group mean variation via lognormal noise (sigma)
    for (g in 1:G) {
      a0g <- (theta_true * od_mult) * N_row[g]
      z <- log(p_true) + rnorm(K, 0, sigma)
      p_g <- exp(z); p_g <- p_g / sum(p_g)
      alpha_vec <- a0g * p_g
      phi <- rgamma(K, shape = alpha_vec, rate = 1); phi <- phi / sum(phi)
      X[g, ] <- as.vector(rmultinom(1, size = N_row[g], prob = phi))
    }
    return(X)
  }

  stop("Unknown gen. Use GEN in {MN, DMpure, DMstrong}.")
}

theta_profile_conditional <- function(X, N_row, p_fixed, theta_grid) {
  # Conditional profile: nll(theta | p_fixed)
  vapply(theta_grid, function(th) {
    -sum(vapply(1:nrow(X), function(g) dm_loglik_group(X[g,], p_fixed, th * N_row[g]), 0.0))
  }, 0.0)
}

# ----------------------------- read input ----------------------------------

args <- commandArgs(trailingOnly = TRUE)
inp_path <- if (length(args) >= 1) args[1] else "validate_hake_selftest.inp"
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

out_prefix <- get_chr("out_prefix", "validate_hake_selftest")
p_true     <- get_vec_num("p_true", K_expected = K, default = rep(1/K, K))
nsamp_levels <- get_vec_num("nsamp_levels", default = c(25, 125, 250))
nsamp_levels <- as.integer(nsamp_levels)

# New optional inputs
theta_lower <- get_num("theta_lower", 1e-6)
theta_upper <- get_num("theta_upper", 1e4)
profile_reps <- get_int("profile_reps", 5L)
profile_grid_n <- get_int("profile_grid_n", 40L)
include_DMpure <- get_int("include_DMpure", 1L)
include_DMstrong <- get_int("include_DMstrong", 1L)

# validation
stopifnot(is.finite(K), K >= 2, is.finite(G), G >= 1)
stopifnot(all(p_true > 0), abs(sum(p_true) - 1) < 1e-8)
stopifnot(all(nsamp_levels >= 1))
stopifnot(is.finite(theta_true) && theta_true > 0)
stopifnot(is.finite(sigma) && sigma > 0)
stopifnot(is.finite(od_mult) && od_mult > 0)
stopifnot(R >= 1)
stopifnot(is.finite(theta_lower) && theta_lower > 0)
stopifnot(is.finite(theta_upper) && theta_upper > theta_lower)
stopifnot(profile_reps >= 0)
stopifnot(profile_grid_n >= 10)

set.seed(random_seed)

gen_levels <- c("MN")
if (isTRUE(include_DMpure == 1L)) gen_levels <- c(gen_levels, "DMpure")
if (isTRUE(include_DMstrong == 1L)) gen_levels <- c(gen_levels, "DMstrong")

params <- list(
  K = K, G = G, p_true = p_true,
  theta_true = theta_true, sigma = sigma, od_mult = od_mult,
  GEN = gen_levels, Nsamp = nsamp_levels, R = R,
  theta_lower = theta_lower, theta_upper = theta_upper,
  profile_reps = profile_reps, profile_grid_n = profile_grid_n,
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
  GEN = gen_levels,
  Nsamp = nsamp_levels,
  rep = seq_len(R),
  stringsAsFactors = FALSE
)

# preallocate results list (faster than rbind)
res <- vector("list", nrow(scen) * 2L)  # *2 fits per replicate
idx <- 0L

# preallocate profile rows (upper bound)
n_profile_scen <- length(gen_levels) * length(nsamp_levels) * min(profile_reps, R)
theta_grid <- exp(seq(log(theta_lower), log(theta_upper), length.out = profile_grid_n))
prof_list <- vector("list", n_profile_scen)
prof_idx <- 0L

# starting values for p (in z-space), centered on truth
z0 <- log(p_true[1:(K-1)]) - log(p_true[K])

# IMPORTANT: the fitted-model "truth" under our generator scaling is theta_true * od_mult
theta_fit_truth <- theta_true * od_mult
t0 <- log(theta_fit_truth)

# ----------------------------- main loop ------------------------------------

for (ii in seq_len(nrow(scen))) {
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
    hit_upper = mn$hit_upper,
    stringsAsFactors = FALSE
  )
  for (k in 1:K) res[[idx]][[paste0("p", k, "_hat")]] <- mn$p_hat[k]

  # FIT = DM (bounded DML with theta estimated)
  dm <- fit_dml_bounded(X, N_row, z0 = z0, t0 = t0,
                       theta_lower = theta_lower, theta_upper = theta_upper, maxit = 500)
  idx <- idx + 1L
  res[[idx]] <- data.frame(
    GEN = gen, Nsamp = Ns, rep = rr, FIT = "DM",
    conv = dm$conv,
    rmse_p = rmse(dm$p_hat, p_true),
    theta_hat = dm$theta_hat,
    hit_upper = dm$hit_upper,
    stringsAsFactors = FALSE
  )
  for (k in 1:K) res[[idx]][[paste0("p", k, "_hat")]] <- dm$p_hat[k]

  # Theta profile subset: first profile_reps reps per (GEN, Nsamp), only for DM fits that converged
  if (profile_reps > 0 && rr <= profile_reps && isTRUE(dm$conv) && all(is.finite(dm$p_hat))) {
    prof_idx <- prof_idx + 1L
    nll_prof <- theta_profile_conditional(X, N_row, p_fixed = dm$p_hat, theta_grid = theta_grid)
    prof_list[[prof_idx]] <- data.frame(
      GEN = gen, Nsamp = Ns, rep = rr,
      theta = theta_grid,
      nll_cond = nll_prof,
      theta_hat = dm$theta_hat,
      hit_upper = dm$hit_upper,
      stringsAsFactors = FALSE
    )
  }
}

raw <- do.call(rbind, res)

# profiles (drop NULLs)
profiles <- if (prof_idx > 0) do.call(rbind, prof_list[seq_len(prof_idx)]) else data.frame()

# ----------------------------- summaries ------------------------------------

# convergence rate + RMSE summaries + theta summaries + upper-bound hit rate (DM only)
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

# upper-bound hit rate for DM fits only
dm_only <- subset(raw, FIT == "DM")
hit_tbl <- aggregate(
  as.numeric(hit_upper) ~ GEN + Nsamp,
  data = dm_only,
  FUN = function(x) mean(x, na.rm = TRUE)
)
names(hit_tbl)[3] <- "theta_upper_hit_rate"
summary_df <- merge(summary_df, hit_tbl, by = c("GEN","Nsamp"), all.x = TRUE, sort = FALSE)

# add the fitted-model truth for reference
summary_df$theta_fit_truth <- theta_fit_truth

# ----------------------------- write outputs --------------------------------

raw_path <- paste0(out_prefix, "_raw.csv")
sum_path <- paste0(out_prefix, "_summary.csv")
prof_path <- paste0(out_prefix, "_theta_profile.csv")

write.csv(raw, raw_path, row.names = FALSE)
write.csv(summary_df, sum_path, row.names = FALSE)
if (nrow(profiles) > 0) write.csv(profiles, prof_path, row.names = FALSE)

cat("Wrote: ", raw_path, "\n", sep = "")
cat("Wrote: ", sum_path, "\n", sep = "")
if (nrow(profiles) > 0) cat("Wrote: ", prof_path, "\n", sep = "")

# Plot: RMSE by FIT, faceted by GEN and Nsamp
p <- ggplot(raw, aes(x = FIT, y = rmse_p)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_grid(GEN ~ Nsamp, scales = "free_y") +
  labs(title = "Self-test: RMSE of p-hat (MN vs bounded DM fits)",
       subtitle = sprintf("Truth p=(%s), G=%d, R=%d; DM bounds [%.1e, %.1e]; theta_fit_truth=%.3f",
                          paste(round(p_true, 3), collapse=", "), G, R, theta_lower, theta_upper, theta_fit_truth),
       x = "Fit likelihood", y = "RMSE(p_hat, p_true)") +
  theme_bw()

plot_path <- paste0(out_prefix, "_rmse_boxplot.png")
ggsave(plot_path, p, width = 9, height = 5.5, dpi = 300)
cat("Wrote: ", plot_path, "\n", sep = "")

# Print a compact summary to console / list file
cat("\n--- Summary (means) ---\n")
print(summary_df[order(summary_df$GEN, summary_df$Nsamp, summary_df$FIT), ])

cat("\nInterpretation tip:\n")
cat("  - If theta_upper_hit_rate is high for a GEN scenario, the fitted DM likelihood is effectively preferring MN.\n")
cat("  - DMPure should recover theta around theta_fit_truth (~theta_true*od_mult) if the optimizer is behaving.\n")
cat("  - DMstrong is deliberately mis-specified relative to the fitted mean structure; large theta can occur even when overdispersion exists.\n")

# Show plot
print(p)

on.exit(sink())
