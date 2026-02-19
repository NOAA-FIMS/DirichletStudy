# hake_ev1.R parameter mapping for models (i)–(iv)

This note maps **estimated parameters** to the **exact objects** and **optimization blocks** in `hake_ev1.R` using **line numbers from the script**.

## Key objects used by all models

### Softmax parameterization for simplex p

`softmax_from_z()` converts an unconstrained vector `z` (length K-1) to a simplex vector `p` (length K).

```text
0157: softmax_from_z <- function(z) { v <- c(z, 0); ev <- exp(v - max(v)); ev / sum(ev) }
```

### Dirichlet–multinomial log-likelihood for one group

Used by models (ii) and (iii) via `dm_loglik_group(x, p, a0)`.

```text
0162: dm_loglik_group <- function(x, p, a0) {
0163:   N <- sum(x)
0164:   lgamma(a0) - lgamma(N + a0) + sum(lgamma(x + a0 * p) - lgamma(a0 * p))
0165: }
```

### Gradient (score) helpers for DM models

- (i): gradient w.r.t. `t = log(alpha_vec)`
- (ii): gradient w.r.t. `z` (for p) and `t_g = log(a0_g)`
- (iii): gradient w.r.t. `z` (for p) and `t = log(theta)`

```text
0167: # scores (i)
0168: score_param_i <- function(t, X) {
0169:   alpha_vec <- exp(t)
0170:   a0    <- sum(alpha_vec)
0171:   g_alpha <- rep(0.0, length(alpha_vec))
0172:   for (g in 1:G) {
0173:     xg <- X[g, ]
0174:     Ng <- sum(xg)
0175:     common <- digamma(a0) - digamma(Ng + a0)         # scalar
0176:     g_alpha <- g_alpha +
0177:       common + (digamma(xg + alpha_vec) - digamma(alpha_vec))# vector add
0178:   }
0179:   g_t <- alpha_vec * g_alpha
0180:   -g_t   # gradient of nll
0181: }
0182: 
0183: ## scores (ii): p shared; a0_g free (opt over z (length K-1) and t_g = log a0_g)
0184: score_param_ii <- function(X, p, a0_vec) {
0185:   G <- nrow(X); K <- length(p)
0186:   g_p <- numeric(K); g_t <- numeric(G)
0187:   for (g in 1:G) {
0188:     xg <- X[g, ]; a0 <- a0_vec[g]; Ng <- sum(xg)
0189:     g_p <- g_p + a0 * (digamma(xg + a0 * p) - digamma(a0 * p))
0190:     g_t[g] <- a0 * ( digamma(a0) - digamma(Ng + a0)
0191:                      + sum(p * (digamma(xg + a0 * p) - digamma(a0 * p))) )
0192:   }
0193:   J <- diag(p) - tcrossprod(p)                     # projection to simplex tangent
0194:   g_z <- as.vector(crossprod(J[, 1:(K-1), drop=FALSE], g_p))
0195:   list(g_z = g_z, g_t = g_t)
0196: }
0197: 
0198: ## scores (iii): p shared; theta shared (opt over z and t = log theta)
0199: score_param_iii <- function(X, p, theta) {
0200:   G <- nrow(X); K <- length(p)
0201:   g_p <- numeric(K); g_theta <- 0
0202:   for (g in 1:G) {
0203:     xg <- X[g, ]; Ng <- sum(xg); a0 <- theta * Ng
0204:     g_p <- g_p + a0 * (digamma(xg + a0 * p) - digamma(a0 * p))
0205:     g_theta <- g_theta + Ng * ( digamma(a0) - digamma(Ng + a0)
0206:                                 + sum(p * (digamma(xg + a0 * p) - digamma(a0 * p))) )
0207:   }
0208:   J <- diag(p) - tcrossprod(p)
0209:   g_z <- as.vector(crossprod(J[, 1:(K-1), drop=FALSE], g_p))
0210:   g_t <- theta * g_theta
0211:   list(g_z = g_z, g_t = g_t)
0212: }
```

## Simulation: where p_true enters alpha_vec

Within each group `g`, the script builds `alpha_vec` and then draws `phi ~ Dirichlet(alpha_vec)` via normalized gamma draws, then draws counts `X[g,] ~ Mult(N_row[g], phi)`.

```text
0280:     p_g <- rgamma(K, shape = kappa_g * p_true, rate = 1)
0281:     p_g <- p_g / sum(p_g)
0282:     
0283:     # propagate into alpha_vec
0284:     alpha_vec <- a0g * p_g
0285:     
0286:     # existing DM step
0287:     phi <- rgamma(K, shape = alpha_vec, rate = 1)
0288:     phi <- phi / sum(phi)
0289:     
0290:     # multinomial observation
0291:     X[g, ] <- as.vector(rmultinom(1, size = N_row[g], prob = phi))
0292:   }
```

Parameter meaning in that block:

- `a0g <- theta_true * N_row[g]`: total Dirichlet concentration for group g (α0_g)
- `alpha_vec <- a0g * p_true`: Dirichlet concentration vector (α_gk)
- `phi <- rgamma(...); phi <- phi/sum(phi)`: a Dirichlet draw on the simplex
- `rmultinom(..., prob=phi)`: multinomial sample given that draw

## Model (i): DM with unconstrained concentration vector α (shared across groups)

### Estimated parameters

- **α (alpha_vec)**: a length-K positive vector, shared for all groups.

### Parameterization used in optimization

- Optimized parameter is **`t`**, where `t = log(alpha_vec)` (unconstrained in R).

### Optimization block and outputs

```text
0295: ## MLE (i): unconstrained alpha_vec (shared across groups) -------------------
0296: ## Start at alpha_vec = alpha0 * p_init with a small positive floor
0297: alpha0_init <- 1.0
0298: pooled <- colSums(X); p_init <- pmax(pooled, 1) / sum(pmax(pooled, 1))
0299: t0_i <- log(pmax(alpha0_init * p_init, 1e-4))   # t = log(alpha_vec), length K
0300: 
0301: ## Negative joint log-likelihood for (i)
0302: nll_i <- function(t,X) {
0303:   alpha_vec <- exp(t)                 # length K
0304:   a0    <- sum(alpha_vec)
0305:   # Sum over groups
0306:   -sum(vapply(1:G, function(g) {
0307:     xg <- X[g, ]
0308:     Ng <- sum(xg)
0309:     lgamma(a0) - lgamma(Ng + a0) +
0310:       sum(lgamma(xg + alpha_vec) - lgamma(alpha_vec))
0311:   }, 0.0))
0312: }
0313:  
0314: fit_i <- optim(t0_i, nll_i, score_param_i, method = "BFGS",
0315:                control = list(maxit = 500, reltol = 1e-10), 
0316:                X = X)
0317: 
0318: alpha_hat_i <- exp(fit_i$par)
0319: p_hat_i     <- alpha_hat_i / sum(alpha_hat_i)
0320: 
```

Line-by-line mapping:

- **0299**: `t0_i` is the starting value for `t = log(alpha_vec)`.

- **0302–0312**: `nll_i(t, X)` computes the **negative log-likelihood** for the DM with α = exp(t).

- **0314–0316**: `optim(t0_i, nll_i, score_param_i, ...)` runs BFGS; the **gradient** is `score_param_i` (defined at **0168–0181**).

- **0318**: `alpha_hat_i <- exp(fit_i$par)` converts the optimizer result back to α.

- **0319**: `p_hat_i <- alpha_hat_i / sum(alpha_hat_i)` is the implied mean composition estimate.

## Model (ii): DM with shared mean p and free α0_g per group (p shared; a0_g free)

### Estimated parameters

- **p**: a length-K simplex vector shared across groups
- **a0_vec**: length-G vector of group-specific total concentrations (α0_g)

### Parameterization used in optimization

- `z` (length K-1): unconstrained parameters mapped to `p = softmax_from_z(z)`
- `t` (length G): unconstrained parameters mapped to `a0_vec = exp(t)`
- Combined into one vector `par = c(z, t)`

### Optimization block and outputs

```text
0321:   ## MLE (ii): (p, a0_g)
0322:   pooled <- colSums(X); p_init <- pmax(pooled, 1) / sum(pmax(pooled, 1))
0323:   z0 <- log(p_init[1:(K-1)]) - log(p_init[K])
0324:   t0 <- rep(0, G)  # log a0_g = 0 => a0_g = 1
0325:   
0326:   nll_ii <- function(par) {
0327:     z <- par[1:(K-1)]; t <- par[K:(K-1+G)]
0328:     p <- softmax_from_z(z); a0 <- exp(t)
0329:     -sum(vapply(1:G, function(g)
0330:       dm_loglik_group(X[g,], p, a0[g]), 0.0))
0331:   }
0332:   grad_ii <- function(par) {
0333:     z <- par[1:(K-1)]; t <- par[K:(K-1+G)]
0334:     p <- softmax_from_z(z); a0 <- exp(t)
0335:     sc <- score_param_ii(X, p, a0)
0336:     -c(sc$g_z, sc$g_t)
0337:   }
0338:   fit_ii <- try(optim(c(z0, t0), nll_ii, grad_ii, method = "BFGS",
0339:                       control = list(maxit = 500, reltol = 1e-10)), silent = TRUE)
0340:   
0341:   if (inherits(fit_ii, "try-error")) {
0342:     p_hat_ii <- rep(NA_real_, K)
0343:   } else {
0344:     zhat <- fit_ii$par[1:(K-1)]
0345:     p_hat_ii <- softmax_from_z(zhat)
0346:   }
0347:   
```

Line-by-line mapping:

- **0323**: `z0` is the starting value for `z` (log-ratio vs last category).

- **0324**: `t0` starts all log(α0_g) at 0 (so α0_g=1).

- **0326–0331**: `nll_ii(par)` builds `p <- softmax_from_z(z)` and `a0 <- exp(t)` and sums DM loglik over g.

- **0332–0337**: `grad_ii(par)` uses `score_param_ii()` (lines **0184–0196**) to return the gradient for BFGS.

- **0338–0339**: `optim(c(z0,t0), nll_ii, grad_ii, ...)` performs the numerical MLE.

- **0344–0346**: `p_hat_ii` is extracted from the optimized `z` using `softmax_from_z()`.

Note: the script records `p_hat_ii`; α0_g estimates live in `fit_ii$par[K:(K-1+G)]` (log-scale) if you want them.

## Model (iii): DM with shared mean p and shared θ (α0_g = θ N_g)

### Estimated parameters

- **p**: length-K simplex vector shared across groups
- **theta**: shared scalar controlling α0_g through α0_g = θ N_g

### Parameterization used in optimization

- `z` (length K-1) → `p = softmax_from_z(z)`
- `t` (scalar) → `theta = exp(t)`
- Combined as `par = c(z, t)`

### Optimization block and outputs

```text
0348:   ## MLE (iii): (p, theta), α0_g = θ N_g
0349:   nll_iii <- function(par) {
0350:     z <- par[1:(K-1)]; t <- par[K]
0351:     p <- softmax_from_z(z); theta <- exp(t)
0352:     -sum(vapply(1:G, function(g)
0353:       dm_loglik_group(X[g,], p, theta * N_vec[g]), 0.0))
0354:   }
0355:   grad_iii <- function(par) {
0356:     z <- par[1:(K-1)]; t <- par[K]
0357:     p <- softmax_from_z(z); theta <- exp(t)
0358:     sc <- score_param_iii(X, p, theta)
0359:     -c(sc$g_z, sc$g_t)
0360:   }
0361:   fit_iii <- try(optim(c(z0, 0), nll_iii, grad_iii, method = "BFGS",
0362:                        control = list(maxit = 500, reltol = 1e-10)), silent = TRUE)
0363:   
0364:   if (inherits(fit_iii, "try-error")) {
0365:     p_hat_iii <- rep(NA_real_, K)
0366:   } else {
0367:     zhat <- fit_iii$par[1:(K-1)]
0368:     p_hat_iii <- softmax_from_z(zhat)
0369:   }
```

Line-by-line mapping:

- **0349–0354**: `nll_iii(par)` uses `a0_g = theta * N_vec[g]` inside `dm_loglik_group()`.

- **0355–0360**: `grad_iii(par)` uses `score_param_iii()` (lines **0199–0212**) to return gradients in z and log(theta).

- **0361–0362**: `optim(c(z0,0), ...)` fits; starting `t=0` gives theta=1.

- **0367–0368**: `p_hat_iii` is extracted from optimized `z`.

Note: `theta_hat <- exp(fit_iii$par[K])` if you want the fitted theta explicitly.

## Model (iv): Multinomial MLE with shared p (no DM overdispersion)

### Estimated parameters

- **p**: length-K simplex vector shared across groups

### Estimation block (closed-form; no optimization)

```text
0370:   ## MLE (iv): multinomial with shared p across groups
0371:   pooled <- colSums(X)
0372:   p_hat_iv <- pooled / sum(pooled)
0373: 
0374:   rmse_i   = rmse(p_hat_i,   p_true)
0375:   rmse_ii  = rmse(p_hat_ii,  p_true)
0376:   rmse_iii = rmse(p_hat_iii, p_true)
0377:   rmse_iv  = rmse(p_hat_iv,  p_true)
```

Line-by-line mapping:

- **0371–0372**: `pooled <- colSums(X)` and `p_hat_iv <- pooled / sum(pooled)` is the multinomial MLE for shared p.

## Where RMSE is computed/stored for all four models

The per-simplex-point routine returns `p_true` followed by `rmse_i..rmse_iv`.

```text
0374:   rmse_i   = rmse(p_hat_i,   p_true)
0375:   rmse_ii  = rmse(p_hat_ii,  p_true)
0376:   rmse_iii = rmse(p_hat_iii, p_true)
0377:   rmse_iv  = rmse(p_hat_iv,  p_true)
0378:   c(p_true,rmse_i,rmse_ii,rmse_iii,rmse_iv)
0379: }
0380: 
0381: # res_mat <- t(apply(mesh, 1, run_one))
0382: res_mat <- do.call(
0383:   rbind,
0384:   lapply(seq_len(nmesh), function(i) {
0385:     run_one(p_true = mesh[i, ], N_row = N_vec[i, ])
0386:   })
0387: )
0388: # Replace your block with this:
0389: 
0390: out <- as.data.frame(res_mat)
0391: 
0392: # K = number of true proportions assumed to be the first K columns of res_mat
0393: K <- ncol(out) - 4L  # (p1..pK, then rmse_i, rmse_ii, rmse_iii, rmse_iv)
0394: 
0395: # Ensure N_vec is a matrix with one row per simulation and G columns (groups)
0396: N_mat <- if (is.null(dim(N_vec))) {
0397:   matrix(N_vec, nrow = nrow(out), byrow = TRUE)
0398: } else {
0399:   N_vec
0400: }
0401: 
0402: # Coerce sample sizes to integers (safe rounding, then integer cast)
0403: N_int <- matrix(as.integer(round(N_mat)),
0404:                 nrow = nrow(N_mat),
0405:                 ncol = ncol(N_mat))
0406: colnames(N_int) <- paste0("N", seq_len(ncol(N_int)))  # N1..NG
0407: 
0408: # Bind G sample-size columns to out
0409: out <- cbind(out, N_int)
0410: 
0411: # Name the existing columns for proportions and RMSEs
0412: colnames(out)[seq_len(K)] <- paste0("p", seq_len(K))
0413: colnames(out)[K + seq_len(4)] <- c("rmse_i","rmse_ii","rmse_iii","rmse_iv")
```

## Where the repeated-measures tests include all 4 models

This block constructs the long data frame with methods i–iv and runs Friedman + pairwise Wilcoxon and other summaries.

```text
0470: n <- length(out$rmse_i)
0471: stopifnot(length(out$rmse_ii) == n, length(out$rmse_iii) == n, length(out$rmse_iv) == n)
0472: 
0473: df <- data.frame(
0474:   id     = rep(seq_len(n), times = 4),                 # block (simplex point)
0475:   method = factor(rep(c("i","ii","iii","iv"), each = n),
0476:                   levels = c("i","ii","iii","iv")),
0477:   rmse   = c(out$rmse_i, out$rmse_ii, out$rmse_iii, out$rmse_iv)
0478: )
0479: 
0480: ## Overall nonparametric repeated-measures test
0481: ft <- friedman.test(rmse ~ method | id, data = df)
0482: print(ft)
0483: 
0484: ## Kendall’s W effect size (approx. from Friedman chi-square)
0485: k <- nlevels(df$method)
0486: W <- as.numeric(ft$statistic) / (n * (k - 1))
0487: cat(sprintf("Kendall's W ≈ %.3f\n", W))
0488: 
0489: ## Post-hoc paired comparisons with multiplicity control
0490: pw <- pairwise.wilcox.test(df$rmse, df$method,
0491:                            paired = TRUE, p.adjust.method = "holm",
0492:                            exact = FALSE)
0493: print(pw)
0494: 
0495: ## Convert to wide for paired, within-simplex-point comparisons
0496: wide <- reshape(df, idvar = "id", timevar = "method", direction = "wide")
0497: 
0498: ## 1) Paired, robust effect: Hodges–Lehmann (median) difference + 95% CI
0499: methods <- c("i","ii","iii","iv")
0500: pair_list <- combn(methods, 2, simplify = FALSE)
0501: 
0502: wilcox_ci <- lapply(pair_list, function(pr) {
0503:   a <- wide[[paste0("rmse.", pr[1])]]
0504:   b <- wide[[paste0("rmse.", pr[2])]]
0505:   wt <- wilcox.test(b, a, paired = TRUE, conf.int = TRUE, exact = FALSE)
0506:   list(comp = paste0(pr[2], " - ", pr[1]), test = wt)
0507: })
0508: 
0509: for (obj in wilcox_ci) {
0510:   cat("
0511: Paired Wilcoxon (", obj$comp, "):
0512: ", sep = "")
0513:   print(obj$test)
0514: }
0515: 
0516: ## 2) Practical impact: percent change in RMSE (median & IQR)
0517: pct <- transform(wide,
0518:                  pct_ii_i    = 100*(rmse.ii  - rmse.i)/rmse.i,
0519:                  pct_iii_i   = 100*(rmse.iii - rmse.i)/rmse.i,
0520:                  pct_iv_i    = 100*(rmse.iv  - rmse.i)/rmse.i,
0521:                  pct_iii_ii  = 100*(rmse.iii - rmse.ii)/rmse.ii,
0522:                  pct_iv_ii   = 100*(rmse.iv  - rmse.ii)/rmse.ii,
0523:                  pct_iv_iii  = 100*(rmse.iv  - rmse.iii)/rmse.iii
0524: )
0525: 
0526: tmp <- sapply(pct[c("pct_ii_i","pct_iii_i","pct_iv_i","pct_iii_ii","pct_iv_ii","pct_iv_iii")],
0527:               \(x) c(median = median(x), IQR = IQR(x),
0528:                      q25 = quantile(x, .25), q75 = quantile(x, .75)))
0529: print(tmp)
0530: 
0531: ## 3) Standardized paired effect r from Wilcoxon Z
0532: r_from_wilcox <- function(x, y){
0533:   wt <- wilcox.test(x, y, paired = TRUE, exact = FALSE)
0534:   z  <- qnorm(wt$p.value/2, lower.tail = FALSE) * sign(median(y - x, na.rm = TRUE))
0535:   r  <- as.numeric(z) / sqrt(sum(!is.na(x) & !is.na(y)))
0536:   return(r)
0537: }
0538: 
0539: r_pairs <- lapply(pair_list, function(pr){
0540:   a <- wide[[paste0("rmse.", pr[1])]]
0541:   b <- wide[[paste0("rmse.", pr[2])]]
0542:   r <- r_from_wilcox(a, b)
0543:   data.frame(pair = paste0(pr[2], " vs ", pr[1]), r = r)
0544: })
0545: r_tbl <- do.call(rbind, r_pairs)
0546: print(r_tbl)
0547: 
0548: m <- lmer(rmse ~ method + (1 | id), data = df)
0549: print(anova(m))                               # overall method effect)
0550: print(emmeans(m, pairwise ~ method, adjust = "holm"))  # post-hoc
0551: 
0552: on.exit(sink()) 
0553: # on.exit(sink(), add = TRUE) 
```
