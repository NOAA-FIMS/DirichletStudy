# file = code-H.R

err_long <- all_wide %>%
  pivot_longer(
    cols = matches("^p_hat_(i|ii|iii|iv|v)_([123])$"),
    names_to = c("method", "cat"),
    names_pattern = "p_hat_(i|ii|iii|iv|v)_([123])",
    values_to = "p_hat"
  ) %>%
  mutate(
    cat = as.integer(cat),
    p_true = case_when(
      cat == 1 ~ p1,
      cat == 2 ~ p2,
      cat == 3 ~ p3
    ),
    error = p_hat - p_true
  )

err_wide <- err_long %>%
  select(example_id, mesh_id, sim_id, dataset_id,
         p1, p2, p3, design_block, G, theta_true, theta_CV,
         sigma, mean_nsamp, method, cat, error) %>%
  pivot_wider(names_from = cat, values_from = error,
              names_prefix = "err")

fit_err_manova <- manova(
  cbind(err1, err2) ~ method * design_block +
    method * G +
    method * theta_true +
    method * theta_CV +
    method * sigma +
    method * mean_nsamp +
    method * p1 + method * p2,
  data = err_wide
)

summary(fit_err_manova, test = "Pillai")



