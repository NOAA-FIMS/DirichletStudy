# file = codeA.R

acc_mv <- acc_long %>%
  group_by(
    example_id, mesh_id, method,
    p1, p2, p3,
    design_block, G, theta_true, theta_CV, sigma,
    mean_nsamp, ln_sd, nb_size, N_range_label
  ) %>%
  summarize(
    rmse = mean(accuracy[metric == "rmse"], na.rm = TRUE),
    L1   = mean(accuracy[metric == "L1_norm"], na.rm = TRUE),
    Linf = mean(accuracy[metric == "Linf_norm"], na.rm = TRUE),
    .groups = "drop"
  )