# file = code-B.R

library(lme4)
library(lmerTest)
library(emmeans)

fit_mv_lmm <- lmer(
  log1p(accuracy) ~ method * metric +
    method * design_block +
    method * G +
    method * theta_true +
    method * theta_CV +
    method * sigma +
    method * mean_nsamp +
    method * p1 + method * p2 +
    (1 | example_id) +
    (1 | simplex_id) +
    (1 | dataset_id),
  data = acc_long
)

anova(fit_mv_lmm)

emmeans(fit_mv_lmm, pairwise ~ method | metric, adjust = "holm")
emmeans(fit_mv_lmm, pairwise ~ method, adjust = "holm")