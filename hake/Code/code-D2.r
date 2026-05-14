# file = code-D2.R

library(vegan)

fit_perm_pair <- adonis2(
  Y ~ method,
  data = acc_mv2,
  permutations = how(blocks = acc_mv2$dataset_id),
  method = "euclidean"
)