# file = metadata_prep.R

# Prepare experiment-level metadata tables 
# named "all_wide" and "acc_long".
# Table "all_wide" is used in analyses: F and H
# Table "acc_long" is used in analyses:

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)

design <- read_csv("hake_factorial_design_matrix.csv") %>%
  mutate(
    example_id   = as.integer(example_id),
    design_block = factor(design_block),
    dist_code    = factor(dist_code),
    G            = factor(G),
    theta_true   = factor(theta_true),
    theta_CV     = factor(theta_CV),
    sigma        = factor(sigma),
    mean_nsamp   = factor(mean_nsamp),
    ln_sd        = factor(ln_sd),
    nb_size      = factor(nb_size),
    N_range_label = factor(N_range_label)
  )

read_one_experiment <- function(example_id) {
  ex_stem <- paste0("hake_ex", example_id)
  ex_dir  <- paste0("ex", example_id)
  f <- file.path(ex_dir, paste0(ex_stem, ".csv"))

  read_csv(f, show_col_types = FALSE) %>%
    mutate(example_id = example_id)
}

all_wide <- map_dfr(design$example_id, read_one_experiment) %>%
  left_join(design, by = "example_id") %>%
  mutate(
    dataset_id = interaction(example_id, mesh_id, sim_id, drop = TRUE),
    simplex_id = interaction(example_id, mesh_id, drop = TRUE)
  )
  
acc_long <- all_wide %>%
  pivot_longer(
    cols = matches("^(rmse|L1_norm|Linf_norm)_(i|ii|iii|iv|v)$"),
    names_to = c("metric", "method"),
    names_pattern = "(rmse|L1_norm|Linf_norm)_(i|ii|iii|iv|v)",
    values_to = "accuracy"
  ) %>%
  mutate(
    method = factor(method, levels = c("i", "ii", "iii", "iv", "v")),
    metric = factor(metric, levels = c("rmse", "L1_norm", "Linf_norm"))
  )