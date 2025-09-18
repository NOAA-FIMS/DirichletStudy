# file = test.R
# test.R applies the functional analysis tool to a set
# of random proportions vectors sampled from a Dirichlet
# multinomial distribution using the manual approach
# with gamma distributions setting the Dirichlet probabilities. 
# The user needs to specify:
# (1) the dimension (K) of the proportion vectors
# (2) the sample size (n) of Dirichlet multinomial proportion vectors to generate
# (3) the concentration parameters of the Dirichlet multinomial distribution

library(DirichletStudy)

# # Load required package
# if (!requireNamespace("ggtern", quietly = TRUE)) {
#   install.packages("ggtern")
# }
# library(ggtern)  # for ternary plot

# Prompt user for number of categories (K), number of samples (n), and alpha parameters
K <- as.integer(readline(prompt = "Enter the number of categories K (must be 3 for plotting): "))
if (K != 3) stop("ggtern visualization only works for K = 3.")

n <- as.integer(readline(prompt = "Enter the number of Dirichlet samples to draw: "))

# Ask for Dirichlet alpha parameters
cat("Enter 3 Dirichlet alpha parameters separated by space (default = 1 1 1): ")
alpha_input <- scan(what = numeric(), nmax = 3, quiet = FALSE)

# Use default alpha = 1 if user enters nothing
if (length(alpha_input) == 0) {
  alpha <- rep(1, 3)
} else {
  alpha <- alpha_input
}

# Manually sample from Dirichlet using gamma distributions
set.seed(123)
gamma_samples <- matrix(0, nrow = n, ncol = K)
for (i in 1:K) {
  gamma_samples[, i] <- rgamma(n, shape = alpha[i], rate = 1)
}
row_sums <- rowSums(gamma_samples)
P <- gamma_samples / row_sums

# Returns an M x K matrix of Dirichlet multinomial counts
dirichlet_multinomial_counts <- t(vapply(
  1:n,
  function(i) as.vector(rmultinom(n = 1, size = if (length(n) == 1) n else n[i],
                                  prob = P[i, ])),
  numeric(K)
))

# Normalize Dirichlet multinomial counts to proportion vectors
dirichlet_multinomial_proportions <- dirichlet_multinomial_counts / if (length(n) == 1) n else n

# quick sanity check:
# stopifnot(all(abs(rowSums(dirichlet_multinomial_proportions) - 1) < 1e-12))

# Save Dirichlet probabilities, Dirichlet multinomial counts and proportions
DM_samples <- list(
  dirichlet_probabilities = P,   # Dirichlet samples you already had
  dirichlet_multinomial_counts = dirichlet_multinomial_counts, # Dirichlet–multinomial counts
  dirichlet_multinomial_proportions  = dirichlet_multinomial_proportions   # normalized DM proportions (what you asked for)
)

# Optionally write to disk
# write.csv(dirichlet_multinomial_proportions,  "dm_props.csv",  row.names = FALSE)
# write.csv(dirichlet_multinomial_counts, "dm_counts.csv", row.names = FALSE)

dirichlet_samples <- dirichlet_multinomial_proportions



# Convert to data frame
samples_df <- as.data.frame(dirichlet_samples)
colnames(samples_df) <- c("A", "B", "C")

print("creating study")
# Create instance of DirichletStudy
dirichlet_study <- new(DirichletStudyInterface)

print(dirichlet_samples)

dirichlet_linear <- new(DirichletLinearInterface)
dirichlet_linear$setSimplexData(dirichlet_samples)
dirichlet_study$addStudy(dirichlet_linear$getId())

dirichlet_saturated <- new(DirichletSaturatedInterface)
dirichlet_saturated$setSimplexData(dirichlet_samples)
dirichlet_study$addStudy(dirichlet_saturated$getId())

print(dirichlet_study$runAnalysis())

dirichlet_linear_results <- dirichlet_saturated$getResults()
dirichlet_linear_results$CentralBoundsDeriviativeCorrelation