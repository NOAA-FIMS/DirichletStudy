# file = test_DM_simplex.R
# test_DM_simplex.R applies the functional analysis tool to a grid
# of proportions vectors sampled from a uniform simplex with mesh=h
# A total of nsamples (nsamples = 1, in this hard-wired example) 
# are randomly generated from a Dirichlet
# multinomial (DM) distribution for each vector in the simplex
# using the manual approach sampling approach with
# the proportion vectors setting the Dirichlet probabilities P[i]
# and the total concentration alpha0 setting the alpha[i] values
# as alpha[i] = alpha0 * P[i]
# The user needs to specify:
# (1) the dimension (K) of the proportion vectors
# (2) the mesh of the proportion simplex (h)
# (3) the sample size (nsamples, fixed at 1) of DM proportion vectors to draw
# (4) the total concentration parameter of the DM distribution (alpha0)

set.seed(26267)

library(DirichletStudy)

# Load required package
if (!requireNamespace("ggtern", quietly = TRUE)) {
  install.packages("ggtern") }
library(ggtern)  # for ternary plot

# Load nexcom function
nexcom_path <- "nexcom.r"
if (!file.exists(nexcom_path)) stop("Required file 'nexcom.r' not found.")
source(nexcom_path)

# Prompt for simplex type and validate
choice <- tolower(readline("Type 'open' or 'closed' for simplex type: "))
if (!choice %in% c("open", "closed")) stop("Invalid choice. Must be 'open' or 'closed'.")

if (choice == "open") {
  K <- as.integer(readline("Enter number of categories (K): "))
  h <- as.numeric(readline("Enter simplex mesh (h in (0,1)): "))
  if (is.na(K) || is.na(h) || K < 1 || h <= 0 || h >= 1) stop("Invalid input. K > 0, 0 < h < 1.")
  
  N <- as.integer(1 / h) - K
  if (N <= 0) stop("Invalid mesh. Must have h < 1/K.")
  
} else {
  N <- as.integer(readline("Enter integer N: "))
  K <- as.integer(readline("Enter number of parts K: "))
  if (is.na(N) || is.na(K) || N < 1 || K < 1) stop("N and K must be positive integers.")
}

nsamples <- as.integer(readline(prompt = "Enter the number of Dirichlet multinomial samples to draw for each point: "))

alpha0 <- as.integer(readline(prompt = "Enter the total concentration parameter: "))

# Generate compositions
compositions <- nexcom(N, K)
compositions_df <- as.data.frame(compositions)
# colnames(compositions_df) <- c("R", "J", "A")
cat(sprintf("Generated %d K-compositions of %d\n", nrow(compositions_df), N))

# Normalize compositions to proportions
if (choice == "open") {
  compositions <- compositions + 1
  proportions <- compositions / (N + K)
  filename <- sprintf("simplex-O_N%d_K%d.csv", N + K, K)
} else {
  proportions <- compositions / N
  filename <- sprintf("simplex-C_N%d_K%d.csv", N, K)
}

ncompositions <- choose((N+K-1),(K-1))

# Manually sample from a Dirichlet multinomial
# Generate gamma samples from concentration parameters
gamma_samples <- matrix(0, nrow = ncompositions, ncol = K)
for (i in 1:K) {
  alpha_col <- alpha0 * proportions[, i]
  gamma_samples[, i] <- rgamma(ncompositions, shape = alpha_col, rate = 1)
}

# Normalize gamma samples to Dirichlet probabilities
row_sums <- rowSums(gamma_samples)
P <- gamma_samples / row_sums

n <- ncompositions

# Return an n x K matrix of Dirichlet multinomial counts
dirichlet_multinomial_counts <- t(vapply(
  1:n,
  function(i) as.vector(rmultinom(n = 1, size = if (length(n) == 1) n else n[i],
                                  prob = P[i, ])),
  numeric(K)
))

# Normalize Dirichlet multinomial counts to proportion vectors
dirichlet_multinomial_proportions <- dirichlet_multinomial_counts / if (length(n) == 1) n else n

# Save Dirichlet probabilities, Dirichlet multinomial counts and proportions
DM_samples <- list(
  dirichlet_probabilities = P,   # Dirichlet probabilities
  dirichlet_multinomial_counts = dirichlet_multinomial_counts, # Dirichlet multinomial counts
  dirichlet_multinomial_proportions  = dirichlet_multinomial_proportions   # normalized DM proportions
)

dirichlet_samples <- dirichlet_multinomial_proportions

# Convert to data frame
samples_df <- as.data.frame(dirichlet_samples)
colnames(samples_df) <- c("R", "J", "A")

if (K == 3) {
  # Plot samples on ternary diagram
  ternary_plot <- ggtern(data = samples_df, aes(x = R, y = J, z = A)) +
    geom_point(alpha = 0.5, size = 0.5) +
    theme_bw(base_size=7) +
    labs(title = "Age Composition of Observed Fishery Samples", T = "Juvenile", L = "Recruit", R = "Adult")
  print(ternary_plot)
}

print("creating study")
# Create instance of DirichletStudy
dirichlet_study <- new(DirichletStudyInterface)

print(dirichlet_samples)

dirichlet_linear <- new(DirichletLinearInterface)
dirichlet_linear$setSimplexData(dirichlet_samples)
dirichlet_study$addStudy(dirichlet_linear$getId())

dirichlet_fisch <- new(DirichletFischInterface)
dirichlet_fisch$setSimplexData(dirichlet_samples)
dirichlet_study$addStudy(dirichlet_fisch$getId())

#dirichlet_thorson <- new(DirichletThorsonInterface)
#dirichlet_thorson$setSimplexData(dirichlet_samples)
#dirichlet_study$addStudy(dirichlet_thorson$getId())

print(dirichlet_study$runAnalysis())



