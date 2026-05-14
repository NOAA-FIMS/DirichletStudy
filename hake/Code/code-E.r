# file = code-E.R

library(mgcv)

fit_gam <- gam(
  log1p(accuracy) ~ metric * method +
    s(p1, p2, by = method, k = 30) +
    design_block + G + theta_true + theta_CV + sigma + mean_nsamp,
  data = acc_long,
  method = "REML"
)

summary(fit_gam)
anova(fit_gam)

diff_v_iv <- make_pairwise_diff(acc_long, "v", "iv")

fit_diff_gam <- gam(
  log1p(abs(diff)) ~ metric +
    s(p1, p2, by = metric, k = 30) +
    design_block + G + theta_true + theta_CV + sigma + mean_nsamp,
  data = diff_v_iv,
  method = "REML"
)

summary(fit_diff_gam)

