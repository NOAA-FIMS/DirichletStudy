# file = code-D1.R

library(vegan)

acc_mv2 <- acc_long %>%
  group_by(example_id, mesh_id, sim_id, method,
           p1, p2, p3, design_block, G, theta_true,
           theta_CV, sigma, mean_nsamp, ln_sd, nb_size,
           N_range_label) %>%
  summarize(
    rmse = accuracy[metric == "rmse"][1],
    L1   = accuracy[metric == "L1_norm"][1],
    Linf = accuracy[metric == "Linf_norm"][1],
    .groups = "drop"
  )

Y <- acc_mv2 %>% select(rmse, L1, Linf)

fit_perm <- adonis2(
  Y ~ method * design_block +
    method * G +
    method * theta_true +
    method * theta_CV +
    method * sigma +
    method * mean_nsamp +
    method * p1 + method * p2,
  data = acc_mv2,
  permutations = 999,
  method = "euclidean"
)

fit_perm