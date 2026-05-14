# file = code-G.R

rank_dat <- acc_long %>%
  group_by(dataset_id, metric) %>%
  mutate(rank = rank(accuracy, ties.method = "average")) %>%
  ungroup()

rank_summary <- rank_dat %>%
  group_by(metric, method) %>%
  summarize(
    mean_rank = mean(rank, na.rm = TRUE),
    prob_best = mean(rank == 1, na.rm = TRUE),
    .groups = "drop"
  )

rank_summary

fit_rank <- lmer(
  rank ~ method * metric +
    design_block + G + theta_true + theta_CV + sigma + mean_nsamp +
    p1 + p2 +
    (1 | example_id) + (1 | dataset_id),
  data = rank_dat
)

anova(fit_rank)
emmeans(fit_rank, pairwise ~ method | metric, adjust = "holm")



