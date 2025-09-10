# file = test_uniform_simplex.R
# test_DM_simplex.R applies the functional analysis tool to a grid
# of proportions vectors sampled from a uniform open simplex with mesh=h.
# The user needs to specify:
# (1) the dimension (K) of the proportion vectors
# (2) the mesh of the proportion simplex (h) with h < 1/K

set.seed(26267)

library(DirichletStudy) # for functional analysis tool
library(ggtern)  # for ternary plot

# Source nexcom function to generate uniform simplex with mesh h
nexcom_path <- "nexcom.r"
if (!file.exists(nexcom_path)) stop("Required file 'nexcom.r' not found.")
source(nexcom_path)

# Input the number of categories K
K <- as.integer(readline("Enter number of categories (K): "))

# Input the simplex mesh h, with 0 < h < 1/K
h <- as.numeric(readline("Enter simplex mesh (h in (0,1)): "))
N <- as.integer(1 / h) - K
if (N <= 0) stop("Invalid mesh. Must have h < 1/K")

# Generate compositions
compositions <- nexcom(N, K)
compositions_df <- as.data.frame(compositions)
colnames(compositions_df) <- c("R", "J", "A")
cat(sprintf("Generated %d K-compositions of %d\n", nrow(compositions_df), N))

# Normalize to proportions
compositions <- compositions + 1
proportions <- compositions / (N + K)
filename <- sprintf("simplex_N%d_K%d.csv", N + K, K)

ncompositions <- choose((N+K-1),(K-1))

ncompositions <- length(proportions[[1L]])
K <- length(proportions)  # should be 3

uniform_simplex <- matrix(
  as.double(unlist(proportions, use.names = FALSE)),
  nrow = ncompositions, ncol = K
)  

# Convert to data frame
samples_df <- as.data.frame(uniform_simplex)
colnames(samples_df) <- c("R", "J", "A")

if (K == 3) {
  # Plot samples on ternary diagram
  ternary_plot <- ggtern(data = samples_df, aes(x = R, y = J, z = A)) +
    geom_point(alpha = 0.5, size = 0.5) +
    theme_bw(base_size=7) +
    labs(title = "Age Composition of Fishery Samples", T = "Juvenile", L = "Recruit", R = "Adult")
  print(ternary_plot)
}

print("creating study")
# Create instance of DirichletStudy
dirichlet_study <- new(DirichletStudyInterface)

print(uniform_simplex)

dirichlet_linear <- new(DirichletLinearInterface)
dirichlet_linear$setSimplexData(uniform_simplex)
dirichlet_study$addStudy(dirichlet_linear$getId())

dirichlet_fisch <- new(DirichletFischInterface)
dirichlet_fisch$setSimplexData(uniform_simplex)
dirichlet_study$addStudy(dirichlet_fisch$getId())

#dirichlet_thorson <- new(DirichletThorsonInterface)
#dirichlet_thorson$setSimplexData(uniform_simplex)
#dirichlet_study$addStudy(dirichlet_thorson$getId())

print(dirichlet_study$runAnalysis())



