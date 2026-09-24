# file = test_runtime_code-A.R

# code-A.R fits a MANOVA with three response variables
# and a relatively large design matrix with interactions 
# between method and several design factors, plus p1 and p2.

# The runtime depends mostly on:
# nrow(acc_mv)
# number of factor levels
# available RAM
# CPU speed
# whether R is 64-bit

library(data.table)

# Load acc_mv
setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN1")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

acc_long <- load_metadata_table("acc_long")
acc_mv   <- make_acc_mv_from_long(acc_long)

# Exact R syntax to estimate runtime, after loading acc_mv, run:
nrow(acc_mv)
object.size(acc_mv)
format(object.size(acc_mv), units = "GB")

# Then check the size of the model matrix:
X <- model.matrix(
  ~ method * design_block +
    method * G +
    method * theta_true +
    method * theta_CV +
    method * sigma +
    method * mean_nsamp +
    method * p1 + method * p2,
  data = acc_mv
)

dim(X)
format(object.size(X), units = "GB")

rm(X)
gc()

# Time a smaller test run. Use a 5% sample first:
set.seed(123)

acc_mv_test <- acc_mv[sample(seq_len(nrow(acc_mv)),
                             size = ceiling(0.05 * nrow(acc_mv))), ]

system.time({
  fit_manova_test <- manova(
    cbind(rmse, L1, Linf) ~ method * design_block +
      method * G +
      method * theta_true +
      method * theta_CV +
      method * sigma +
      method * mean_nsamp +
      method * p1 + method * p2,
    data = acc_mv_test
  )
})

# Full run with timing
gc()

system.time({
  fit_manova <- manova(
    cbind(rmse, L1, Linf) ~ method * design_block +
      method * G +
      method * theta_true +
      method * theta_CV +
      method * sigma +
      method * mean_nsamp +
      method * p1 + method * p2,
    data = acc_mv
  )
})

system.time({
  print(summary(fit_manova, test = "Pillai"))
})

system.time({
  print(summary.aov(fit_manova))
})

# If acc_mv has fewer than about 500,000 rows, 
# code-A.R may finish in a few minutes. 
# If it has millions of rows, it may take tens 
# of minutes and may still hit memory limits. 
# The safest next step is to run the model.matrix() 
# size check and the 5% timing test before running the full MANOVA.

