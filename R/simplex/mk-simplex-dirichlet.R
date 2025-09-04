# file = mk-simplex-dirichlet.R
# mk-simplex-dirichlet.R creates
# random proportions vectors from a Dirichlet
# multinomial distribution using a Dirichlet prior
# for the multinomial proportions.
# Dirichlet multinomial proportions are plotted
# and output to file=dirichlet_samples.csv

if (!requireNamespace("MCMCpack", quietly = TRUE)) {
  install.packages("MCMCpack")
}
if (!requireNamespace("ggtern", quietly = TRUE)) {
  install.packages("ggtern")
}
library(MCMCpack)  # for rdirichlet()
library(ggtern)    # for ternary plots

# Prompt user for number of categories (K), number of samples (n), and Dirichlet parameters
K <- as.integer(readline(prompt = "Enter the number of categories K (must be 3 for plotting): "))
if (K != 3) stop("ggtern visualization only works for K = 3.")

n <- as.integer(readline(prompt = "Enter the number of Dirichlet samples to draw: "))

# Ask for Dirichlet alpha parameters
cat("Enter 3 Dirichlet alpha parameters separated by space (default = 1 1 1): ")
alpha_input <- scan(what = numeric(), nmax = 3, quiet = TRUE)

# Use default alpha = 1 if user enters nothing
if (length(alpha_input) == 0) {
  alpha <- rep(1, 3)
} else {
  alpha <- alpha_input
}

# Sample from Dirichlet distribution
set.seed(123)
P <- rdirichlet(n, alpha)

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
  dirichlet_probabilities = P,   # Dirichlet samples you already had
  dirichlet_multinomial_counts = dirichlet_multinomial_counts, # Dirichlet–multinomial counts
  dirichlet_multinomial_proportions  = dirichlet_multinomial_proportions   # normalized DM proportions (what you asked for)
)

# Optionally write to disk
# write.csv(dirichlet_multinomial_proportions,  "dm_props.csv",  row.names = FALSE)
# write.csv(dirichlet_multinomial_counts, "dm_counts.csv", row.names = FALSE)

dirichlet_samples <- dirichlet_multinomial_proportions

dirichlet_samples_df <- as.data.frame(dirichlet_samples)
colnames(dirichlet_samples_df) <- c("R", "J", "A")

# Save samples to CSV
write.csv(dirichlet_samples_df, file = "dirichlet_samples.csv", row.names = FALSE)
cat("Saved samples to 'dirichlet_samples.csv'\n")

# Plot samples on ternary diagram
ternary_plot <- ggtern(data = dirichlet_samples_df, aes(x = R, y = J, z = A)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_bw(base_size=7) +
  labs(title = "Age Composition of Observed Fishery Samples", T = "Juvenile", L = "Recruit", R = "Adult")
print(ternary_plot)