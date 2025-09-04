# file=compare-M-DM-sampling.R
# compare-M-DM-sampling.R produces a ternary plot
# comparing multinomial and dirichlet multinomial
# proportion sampling with observation error for
# for a user-input proportion vector comprised
# of recruits, juvenile, and adult fish

read_positive_number <- function(prompt, default) {
  repeat {
    ans <- readline(paste0(prompt, " [", default, "]: "))
    val <- if (nzchar(ans)) suppressWarnings(as.numeric(ans)) else as.numeric(default)
    if (is.finite(val) && val > 0) return(val)
    cat("Please enter a positive number.\n")
  }
}

parse_true_pi <- function(default_str = "0.3,0.5,0.2") {
  repeat {
    ans <- readline(paste0(
      "Enter true_pi as 3 comma-separated values ",
      "[", default_str, "]: "
    ))
    x <- if (nzchar(ans)) ans else default_str
    x <- gsub("\\s", "", x)
    parts <- strsplit(x, ",", fixed = TRUE)[[1]]
    
    # must be three parts
    if (length(parts) != 3L) {
      cat("Please enter exactly 3 values.\n")
      next
    }
    
    # extract values and (optional) names
    vals <- suppressWarnings(sapply(parts, function(p) {
      if (grepl("=", p, fixed = TRUE)) as.numeric(sub(".*=", "", p)) else as.numeric(p)
    }))
    if (any(is.na(vals))) {
      cat("All entries must be numeric (after '=' if using names).\n")
      next
    }
    if (any(vals < 0)) {
      cat("Values must be nonnegative.\n")
      next
    }
    nms <- sapply(parts, function(p) {
      if (grepl("=", p, fixed = TRUE)) sub("=.*", "", p) else NA_character_
    })
    if (all(is.na(nms))) nms <- c("Recruits","Juveniles","Adults")
    
    # normalize to sum 1 if necessary
    s <- sum(vals)
    if (s <= 0) {
      cat("At least one value must be > 0.\n")
      next
    }
    if (abs(s - 1) > 1e-10) {
      message(sprintf("Note: normalizing true_pi to sum to 1 (current sum = %.6f).", s))
      vals <- vals / s
    }
    names(vals) <- nms
    return(vals)
  }
}

# ---- Gather inputs interactively ----
true_pi <- parse_true_pi("0.3,0.5,0.2")
alpha   <- read_positive_number("Enter alpha (Dirichlet concentration)", 20)
N       <- read_positive_number("Enter N (sample size per draw)", 200)
n_sim   <- read_positive_number("Enter n_sim (number of simulated samples)", 100)

# Load required packages
library(gtools)     # for rdirichlet
library(ggtern)     # for ternary plots
library(tibble)
library(dplyr)

# Set seed
set.seed(3451)

age_classes <- c("Recruits", "Juveniles", "Adults")

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
       T = "Juvenile", L = "Recruit", R = "Adult") +
  theme_bw() +
  theme(
    tern.axis.title.T = element_text(size = rel(0.75)),
    tern.axis.title.L = element_text(size = rel(0.75)),
    tern.axis.title.R = element_text(size = rel(0.75))
  )

print(ternplot)

ggsave(filename = "M vs DM sampling.png", plot = ternplot, width = 7, height = 6.5, dpi = 300)
invisible(ternplot)
