# file = code-C1.R

make_pairwise_diff <- function(dat, m1, m2) {
  dat %>%
    select(example_id, mesh_id, sim_id, dataset_id,
           p1, p2, p3,
           design_block, G, theta_true, theta_CV, sigma,
           mean_nsamp, ln_sd, nb_size, N_range_label,
           method, metric, accuracy) %>%
    filter(method %in% c(m1, m2)) %>%
    pivot_wider(names_from = method, values_from = accuracy) %>%
    mutate(
      comparison = paste0(m1, "_minus_", m2),
      diff = .data[[m1]] - .data[[m2]]
    ) %>%
    select(-all_of(c(m1, m2)))
}

############################## ii_i ############################## 
diff_ii_i <- make_pairwise_diff(acc_long, "ii", "i")

diff_ii_i_mv <- diff_ii_i %>%
  pivot_wider(names_from = metric, values_from = diff) %>%
  rename(L1 = L1_norm, Linf = Linf_norm)

fit_hotelling <- manova(cbind(rmse, L1, Linf) ~ 1, data = diff_ii_i_mv)
summary(fit_hotelling, test = "Pillai")

############################## iii_i ############################## 
diff_iii_i <- make_pairwise_diff(acc_long, "iii", "i")

diff_iii_i_mv <- diff_iii_i %>%
  pivot_wider(names_from = metric, values_from = diff) %>%
  rename(L1 = L1_norm, Linf = Linf_norm)

fit_hotelling <- manova(cbind(rmse, L1, Linf) ~ 1, data = diff_iii_i_mv)
summary(fit_hotelling, test = "Pillai")

############################## iv_i ############################## 
diff_iv_i <- make_pairwise_diff(acc_long, "iv", "i")

diff_iv_i_mv <- diff_iv_i %>%
  pivot_wider(names_from = metric, values_from = diff) %>%
  rename(L1 = L1_norm, Linf = Linf_norm)

fit_hotelling <- manova(cbind(rmse, L1, Linf) ~ 1, data = diff_iv_i_mv)
summary(fit_hotelling, test = "Pillai")

############################## v_i ############################## 
diff_v_i <- make_pairwise_diff(acc_long, "v", "i")

diff_v_i_mv <- diff_v_i %>%
  pivot_wider(names_from = metric, values_from = diff) %>%
  rename(L1 = L1_norm, Linf = Linf_norm)

fit_hotelling <- manova(cbind(rmse, L1, Linf) ~ 1, data = diff_v_i_mv)
summary(fit_hotelling, test = "Pillai")

############################## iii_ii ############################## 
diff_iii_ii <- make_pairwise_diff(acc_long, "iii", "ii")

diff_iii_ii_mv <- diff_iii_ii %>%
  pivot_wider(names_from = metric, values_from = diff) %>%
  rename(L1 = L1_norm, Linf = Linf_norm)

fit_hotelling <- manova(cbind(rmse, L1, Linf) ~ 1, data = diff_iii_ii_mv)
summary(fit_hotelling, test = "Pillai")

############################## iv_ii ############################## 
diff_iv_ii <- make_pairwise_diff(acc_long, "iv", "ii")

diff_iv_ii_mv <- diff_iv_ii %>%
  pivot_wider(names_from = metric, values_from = diff) %>%
  rename(L1 = L1_norm, Linf = Linf_norm)

fit_hotelling <- manova(cbind(rmse, L1, Linf) ~ 1, data = diff_iv_ii_mv)
summary(fit_hotelling, test = "Pillai")

############################## v_ii ############################## 
diff_v_ii <- make_pairwise_diff(acc_long, "v", "ii")

diff_v_ii_mv <- diff_v_ii %>%
  pivot_wider(names_from = metric, values_from = diff) %>%
  rename(L1 = L1_norm, Linf = Linf_norm)

fit_hotelling <- manova(cbind(rmse, L1, Linf) ~ 1, data = diff_v_ii_mv)
summary(fit_hotelling, test = "Pillai")

############################## iv_iii ############################## 
diff_iv_iii <- make_pairwise_diff(acc_long, "iv", "iii")

diff_iv_iii_mv <- diff_iv_iii %>%
  pivot_wider(names_from = metric, values_from = diff) %>%
  rename(L1 = L1_norm, Linf = Linf_norm)

fit_hotelling <- manova(cbind(rmse, L1, Linf) ~ 1, data = diff_iv_iii_mv)
summary(fit_hotelling, test = "Pillai")

############################## v_iii ############################## 
diff_v_iii <- make_pairwise_diff(acc_long, "v", "iii")

diff_v_iii_mv <- diff_v_iii %>%
  pivot_wider(names_from = metric, values_from = diff) %>%
  rename(L1 = L1_norm, Linf = Linf_norm)

fit_hotelling <- manova(cbind(rmse, L1, Linf) ~ 1, data = diff_v_iii_mv)
summary(fit_hotelling, test = "Pillai")

############################## v_iv ############################## 
diff_v_iv <- make_pairwise_diff(acc_long, "v", "iv")

diff_v_iv_mv <- diff_v_iv %>%
  pivot_wider(names_from = metric, values_from = diff) %>%
  rename(L1 = L1_norm, Linf = Linf_norm)

fit_hotelling <- manova(cbind(rmse, L1, Linf) ~ 1, data = diff_v_iv_mv)
summary(fit_hotelling, test = "Pillai")

##################################################################