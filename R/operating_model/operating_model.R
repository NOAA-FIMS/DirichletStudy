# operating_model.R
# Age-structured operating model with plus-group and constant fully-selected F
# Reads parameters from a CSV file with columns: parameter,value,comment
# Required parameters (see om_parameters.csv): Y, A, Linf, K, t0, lw_a, lw_b,
# mat_slope, mat_a50, M, pfem, R0, h, fish_sel_slope, fish_sel_a50, sur_sel_slope, sur_sel_a50,
# F_full, q_survey (optional in this script)

# --- Utilities ----
logistic <- function(x, a50, slope) 1 / (1 + exp(-slope * (x - a50)))

# Survivors per recruit l_a given age-specific total mortality Z_a
spr_survivors <- function(Z) {
  A <- length(Z)
  l <- numeric(A)
  l[1] <- 1.0
  if (A > 1) {
    for (a in 1:(A - 2)) l[a + 1] <- l[a] * exp(-Z[a])
    l[A] <- l[A-1]*exp(-Z[A-1])/(1-exp(-Z[A]))
  }
  l
}

# Spawning biomass per recruit given Z, maturity m, weight W, with plus-group
spr_given_Z <- function(Z, m, W) {
  A <- length(Z)
  l <- spr_survivors(Z)
  spr_no_plus <- sum(l[1:( A - 1)] * m[1:(A - 1)] * W[1:(A - 1)])
  # plus-group contribution at age A
  plus_term <- (l[A] * m[A] * W[A]) / (1 - exp(-Z[A]))
  spr_no_plus + plus_term
}

# Baranov catch in numbers at age
catch_numbers <- function(N_at_age, F_full, sel, Z_at_age) {
  F_at_age <- F_full * sel
  N_at_age * (F_at_age / Z_at_age) * (1 - exp(-Z_at_age))
}

# Beverton-Holt parameters (alpha, beta) from (R0, h, phi0)
bh_from_steepness <- function(R0, h, phi0) {
  # alpha and beta satisfying R = (alpha * SSB) / (1 + beta * SSB)
  alpha <- (4 * h) / (phi0 * (1 - h))
  beta  <- (5 * h - 1) / (R0*phi0 * (1 - h))
  list(alpha = alpha, beta = beta)
}

# Beverton-Holt expected recruitment given SSB and (alpha, beta)
bh_expectation <- function(SSB, alpha, beta) (alpha * SSB) / (1 + beta * SSB)

# --- Main OM function ----
run_operating_model <- function(param_csv = "om_parameters.csv") {
  # Read and coerce parameters
  df <- utils::read.csv(param_csv, stringsAsFactors = FALSE)
  df$value <- utils::type.convert(df$value, as.is = TRUE)
  p <- stats::setNames(as.list(df$value), df$parameter)

  # Pull required parameters with minimal checks
  Y <- as.integer(p$Y); A <- as.integer(p$A)
  Linf <- p$Linf; K <- p$K; t0 <- p$t0
  lw_a <- p$lw_a; lw_b <- p$lw_b
  mat_slope <- p$mat_slope; mat_a50 <- p$mat_a50
  M <- p$M; pfem <- p$pfem
  R0 <- p$R0; h <- p$h
  fish_sel_slope <- p$fish_sel_slope; fish_sel_a50 <- p$fish_sel_a50
  sur_sel_slope <- p$sur_sel_slope; sur_sel_a50 <- p$sur_sel_a50
  F_full <- p$F_full
  q_survey <- if (!is.null(p$q_survey)) p$q_survey else NA_real_

  ages <- seq_len(A)
  years <- seq_len(Y)

  # --- Life history at age (Equations 1.1–1.3) ---
  L_at_age <- Linf * (1 - exp(-K * (ages - t0)))         # 1.1 von Bertalanffy
  W_at_age <- lw_a * (L_at_age ^ lw_b)                   # 1.2 weight-at-length
  mat_at_age <- logistic(ages, mat_a50, mat_slope)       # 1.3 maturity ogive

  # --- Selectivity at age (logistic) ---
  sel_fish <- logistic(ages, fish_sel_a50, fish_sel_slope)
  sel_survey <- logistic(ages, sur_sel_a50, sur_sel_slope)

  # --- Unfished per-recruit quantities (2.1) ---
  Z0 <- rep(M, A)
  phi0 <- spr_given_Z(Z0, mat_at_age, W_at_age) * pfem    # SPR0 in female biomass
  SSB0 <- R0 * phi0                                       # Unfished spawning biomass

  # --- Beverton-Holt recruitment mapping (2.2) ---
  bh <- bh_from_steepness(R0, h, phi0)
  rec_exp <- function(SSB) bh_expectation(SSB, bh$alpha, bh$beta)

  # --- Equilibrium at constant F (3.1–3.4) ---
  ZF <- M + F_full * sel_fish                              # 3.1 total mortality at age
  phiF <- spr_given_Z(ZF, mat_at_age, W_at_age) * pfem     # 3.3 SBR(F) in female biomass
  # Equilibrium recruitment given F (closed form under BH with SSB = R * phiF)
  Req <- max(0, (bh$alpha * phiF - 1) / (bh$beta * phiF))  # 3.4 equilibrium recruitment

  # Survivors per recruit under ZF (3.2)
  lF <- spr_survivors(ZF)

  # Initial numbers at age with plus-group (3.5)
  N_init <- lF * Req
  N_init[A] <- (lF[A] * Req) / (1 - exp(-ZF[A]))           # plus-group aggregate at age A

  # Initial spawning biomass (3.6)
  SSB_init <- pfem * sum(N_init * mat_at_age * W_at_age)

  # --- Time dynamics (Section 4) ---
  # Store numbers at age at start of each year t (columns 1..Y)
  N <- matrix(0.0, nrow = A, ncol = Y)
  Z <- matrix(0.0, nrow = A, ncol = Y)
  C <- matrix(0.0, nrow = A, ncol = Y)                     # fishery catch at age
  P_fish <- matrix(0.0, nrow = A, ncol = Y)    # Proportion of fishery catch at age
  P_survey <- matrix(0.0, nrow = A, ncol = Y)    # Proportion of survey catch at age
  S <- matrix(0.0, nrow = A, ncol = Y)                     # survey catch at age
  SSB <- numeric(Y)
  Catch_weight <- numeric(Y)              # Catch weight
  Survey_index_numbers <- numeric(Y)      # Survey index in numbers
  Survey_index_weight <- numeric(Y)       # Survey index in weight
  TotN <- numeric(Y)
  TotB <- numeric(Y)
  TotC <- numeric(Y)
  TotS <- numeric(Y)
  R <- numeric(Y + 1)                                      # recruitment from t=1..Y+1

  N[, 1] <- N_init
  R[1] <- Req

  for (t in years) {
    # 4.2 Total mortality (constant over years here but stored for clarity)
    Z[, t] <- ZF

    # 5.2 Fishery catch numbers at age a and time t
    C[, t] <- catch_numbers(N[, t], F_full, sel_fish, Z[, t])
    
    # 5.3 Catch weight at time t
    Catch_weight[t] <- sum(C[, t] * W_at_age)
    
    # 6.3 Survey catch numbers at age
    S[, t] <- N[, t] * sel_survey * q_survey
    
    # 6.4 Survey numbers index
    Survey_index_numbers[t] <- sum(S[, t])
    
    # 6.4 Survey weight index
    Survey_index_weight[t] <- sum(S[, t] * W_at_age)

    # 4.5 Spawning biomass at time t
    SSB[t] <- pfem * sum(N[, t] * mat_at_age * W_at_age)

    # 4.6–4.7 Totals
    TotN[t] <- sum(N[, t])
    TotB[t] <- sum(N[, t] * W_at_age)

    # 8.2 Proportion of fishery catch at age    
    TotC[t] <- sum(C[, t])
    if (TotC[t] > 0) P_fish[, t] <- C[, t] / TotC[t]
    
    # 8.5 Proportion of survey catch at age    
    TotS[t] <- sum(S[, t])
    if (TotS[t] > 0) P_survey[, t] <- S[, t] / TotS[t] 

    # 4.1 Expected recruitment at time t+1 (Beverton-Holt)
    R[t + 1] <- rec_exp(SSB[t])

    # 4.3–4.4 Numbers update to start of year t+1 (if not last year)
    if (t < Y) {
      N_next <- numeric(A)
      # age 1 receives recruits
      N_next[1] <- R[t + 1]
      # ages 2..A-1
      if (A > 2) {
        for (a in 2:(A - 1)) N_next[a] <- N[a - 1, t] * exp(-Z[a - 1, t])
      } else if (A == 2) {
        N_next[2] <- N[1, t] * exp(-Z[1, t])
      }
      # plus-group at age A receives survivors from A-1 and stays as a plus-group
      N_next[A] <- N[A - 1, t] * exp(-Z[A - 1, t]) + N[A, t] * exp(-Z[A, t])

      N[, t + 1] <- N_next
    }
  }

  list(
    input_parameters = p,
    ages = ages, years = years,
    life_history = list(L_at_age = L_at_age, W_at_age = W_at_age, 
                        mat_at_age = mat_at_age),
    selectivity = list(fishery = sel_fish, survey = sel_survey),
    per_recruit = list(phi0 = phi0, SSB0 = SSB0, phiF = phiF, 
                       Req = Req, lF = lF),
    initial_state = list(N_at_age = N_init, SSB = SSB_init),
    dynamics = list(R = R, Z_at_age = Z, N_at_age = N, C_at_age = C,
                    Catch_weight = Catch_weight, S_at_age = S, SSB = SSB, 
                    TotN = TotN, TotB = TotB, Survey_N = Survey_index_numbers,
                    Survey_W = Survey_index_weight, Prop_C_at_age = P_fish, 
                    Prop_S_at_age = P_survey)
  )
}

# --- Pretty printer for quick validation ---
print_om_summary <- function(om) {
  cat("Operating Model Summary\n")
  cat(sprintf("  Years (Y): %d | Ages (A): %d\n", length(om$years), length(om$ages)))
  cat(sprintf("  SSB0 (unfished): %.4f\n", om$per_recruit$SSB0))
  cat(sprintf("  Unfished SSB per recruit (phi0): %.4f\n", om$per_recruit$phi0))
  cat(sprintf("  SBR(F): %.4f\n", om$per_recruit$phiF))
  cat(sprintf("  Equilibrium recruitment Req(F): %.4f\n", om$per_recruit$Req))
  cat(sprintf("  Initial SSB: %.4f\n", om$initial_state$SSB))
  cat("  First-year totals:\n")
  cat(sprintf("    SSB[1]=%.4f, TotN[1]=%.4f, TotB[1]=%.4f\n",
              om$dynamics$SSB[1], om$dynamics$TotN[1], om$dynamics$TotB[1]))
}

# Example (uncomment to run in R):
# om <- run_operating_model("om_parameters.csv")
# print_om_summary(om)
# str(om$dynamics, max.level = 1)
