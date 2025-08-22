# Load required packages
library(gtools)     # for rdirichlet
library(ggtern)     # for ternary plots
library(tibble)
library(dplyr)

# Set seed
set.seed(123)

# Define true proportions
true_pi <- c(Recruits = 0.3, Juveniles = 0.5, Adults = 0.2)
age_classes <- names(true_pi)
N <- 200
alpha <- 20
n_sim <- 100

# Function to simulate one sample
simulate_multinom <- function() {
  p <- rmultinom(1, size = N, prob = true_pi)
  as.numeric(p) / N
}

simulate_dm <- function() {
  theta <- as.numeric(rdirichlet(1, alpha * true_pi))
  p <- rmultinom(1, size = N, prob = theta)
  as.numeric(p) / N
}

# Simulate samples
multinom_samples <- replicate(n_sim, simulate_multinom())
dm_samples <- replicate(n_sim, simulate_dm())

# Combine into a data frame
df_multinom <- as.data.frame(t(multinom_samples))
df_dm <- as.data.frame(t(dm_samples))

colnames(df_multinom) <- colnames(df_dm) <- age_classes

df_multinom$Method <- "Multinomial"
df_dm$Method <- "Dirichlet-Multinomial"

df_all <- bind_rows(df_multinom, df_dm)

# Add true proportions
df_true <- as_tibble_row(true_pi)
df_true$Method <- "True"

# Combine all
df_plot <- bind_rows(df_all, df_true)

# Plot: map shape to Method, draw samples, then overlay a larger star for True
ternplot <- ggtern(df_plot, aes(x = Recruits, y = Juveniles, z = Adults,
                                color = Method, shape = Method)) +
  # samples (exclude the True row)
  geom_point(data = dplyr::filter(df_plot, Method != "True"),
             size = 0.5, alpha = 0.6) +
  # big star for the true composition
  geom_point(data = dplyr::filter(df_plot, Method == "True"),
             size = 5, alpha = 1) +
  scale_color_manual(values = c("True" = "black",
                                "Multinomial" = "blue",
                                "Dirichlet-Multinomial" = "red")) +
  scale_shape_manual(values = c("True" = 8,                 
                                "Multinomial" = 16,         # ● filled circle
                                "Dirichlet-Multinomial" = 16)) +
  labs(title = "Age Composition with Sampling Error",
       subtitle = paste(n_sim, "samples with N =", N, "fish"),
       T = "Recruit", L = "Juv", R = "Adult") +
  theme_bw() +
  theme(
    tern.axis.title.T = element_text(size = rel(0.75)),
    tern.axis.title.L = element_text(size = rel(0.75)),
    tern.axis.title.R = element_text(size = rel(0.75))
  )

print(ternplot)

ggsave(filename = "M vs DM sampling.png", plot = ternplot, width = 7, height = 6.5, dpi = 300)
invisible(ternplot)
