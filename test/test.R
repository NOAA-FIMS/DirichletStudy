library(DirichletStudy)

# Load required package
if (!requireNamespace("ggtern", quietly = TRUE)) {
  install.packages("ggtern")
}
library(ggtern)  # for ternary plot

# Prompt user for number of categories (K), number of samples (n), and alpha parameters
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

# Manually sample from Dirichlet using gamma distributions
set.seed(123)
gamma_samples <- matrix(0, nrow = n, ncol = K)
for (i in 1:K) {
  gamma_samples[, i] <- rgamma(n, shape = alpha[i], rate = 1)
}
row_sums <- rowSums(gamma_samples)
dirichlet_samples <- gamma_samples / row_sums

# Convert to data frame
samples_df <- as.data.frame(dirichlet_samples)
colnames(samples_df) <- c("A", "B", "C")

print("creating study")
# Create instance of DirichletStudy
dirichlet_study <- new(DirichletStudyInterface)

# Add studies components to DirichletStudy
dirichlet_default <- new(DirichletDefaultInterface)
dirichlet_default$setSimplexData(dirichlet_samples)
dirichlet_study$addStudy(dirichlet_default$getId())

dirichlet_linear <- new(DirichletLinearInterface)
dirichlet_linear$setSimplexData(dirichlet_samples)
dirichlet_study$addStudy(dirichlet_linear$getId())

dirichlet_fisch <- new(DirichletFischInterface)
dirichlet_fisch$setSimplexData(dirichlet_samples)
dirichlet_study$addStudy(dirichlet_fisch$getId())

print(dirichlet_study$runAnalysis())