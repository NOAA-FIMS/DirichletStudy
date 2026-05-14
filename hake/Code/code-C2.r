# file = code-C2.R

############################## ii_i ##############################
fit_diff_manova <- manova(
  cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
    theta_CV + sigma + mean_nsamp + p1 + p2,
  data = diff_ii_i_mv
)

summary(fit_diff_manova, test = "Pillai")

############################## iii_i ##############################
fit_diff_manova <- manova(
  cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
    theta_CV + sigma + mean_nsamp + p1 + p2,
  data = diff_iii_i_mv
)

summary(fit_diff_manova, test = "Pillai")

############################## iv_i ##############################
fit_diff_manova <- manova(
  cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
    theta_CV + sigma + mean_nsamp + p1 + p2,
  data = diff_iv_i_mv
)

summary(fit_diff_manova, test = "Pillai")

############################## v_i ##############################
fit_diff_manova <- manova(
  cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
    theta_CV + sigma + mean_nsamp + p1 + p2,
  data = diff_v_i_mv
)

summary(fit_diff_manova, test = "Pillai")

############################## iii_ii ##############################
fit_diff_manova <- manova(
  cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
    theta_CV + sigma + mean_nsamp + p1 + p2,
  data = diff_iii_ii_mv
)

summary(fit_diff_manova, test = "Pillai")

############################## iv_ii ##############################
fit_diff_manova <- manova(
  cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
    theta_CV + sigma + mean_nsamp + p1 + p2,
  data = diff_iv_ii_mv
)

summary(fit_diff_manova, test = "Pillai")

############################## v_ii ##############################
fit_diff_manova <- manova(
  cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
    theta_CV + sigma + mean_nsamp + p1 + p2,
  data = diff_v_ii_mv
)

summary(fit_diff_manova, test = "Pillai")

############################## iv_iii ##############################
fit_diff_manova <- manova(
  cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
    theta_CV + sigma + mean_nsamp + p1 + p2,
  data = diff_iv_iii_mv
)

summary(fit_diff_manova, test = "Pillai")

############################## v_iii ##############################
fit_diff_manova <- manova(
  cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
    theta_CV + sigma + mean_nsamp + p1 + p2,
  data = diff_v_iii_mv
)

summary(fit_diff_manova, test = "Pillai")

############################## v_iv ##############################
fit_diff_manova <- manova(
  cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
    theta_CV + sigma + mean_nsamp + p1 + p2,
  data = diff_v_iv_mv
)

summary(fit_diff_manova, test = "Pillai")

##################################################################
