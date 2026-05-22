# file = code-F.R

acc_profile <- all_wide %>%
  select(
    example_id, mesh_id, sim_id, dataset_id,
    p1, p2, p3, design_block, G, theta_true, theta_CV,
    sigma, mean_nsamp, ln_sd, nb_size, N_range_label,
    matches("^(rmse|L1_norm|Linf_norm)_(i|ii|iii|iv|v)$")
  )

Y <- acc_profile %>%
  select(matches("^(rmse|L1_norm|Linf_norm)_(i|ii|iii|iv|v)$")) %>%
  mutate(across(everything(), log1p))

pca <- prcomp(Y, center = TRUE, scale. = TRUE)

pc_dat <- bind_cols(
  acc_profile %>% select(example_id, mesh_id, sim_id, p1, p2, p3,
                         design_block, G, theta_true, theta_CV,
                         sigma, mean_nsamp, N_range_label),
  as.data.frame(pca$x[, 1:5])
)

fit_pc1 <- lm(PC1 ~ design_block + G + theta_true + theta_CV +
                sigma + mean_nsamp + p1 + p2,
              data = pc_dat)

anova(fit_pc1)
summary(fit_pc1)

