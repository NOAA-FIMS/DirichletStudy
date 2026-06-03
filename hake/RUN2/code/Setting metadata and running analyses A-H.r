# file = Setting metadata and running analyses A-H.R
####################################################################################
# A key detail: after you have already built metadata_output/
# Do not let metadata.R rerun the build when you only want loader functions.
# Use HAKE_RUN_ON_SOURCE = "false" before sourcing it.

# Below is the exact R syntax to use to access the metadata in R
# after the partitioned metadata files have already been written to metadata_output/

# The key step is to source metadata.R without rerunning the metadata build.
####################################################################################
setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN2")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

# This loads the helper functions, especially:
# load_metadata_table()
# make_acc_mv_from_long()

####################################################################################
# Analysis A: MANOVA on multivariate accuracy
# Analysis A uses acc_mv, not all_wide directly. acc_mv is built from acc_long.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN2")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

acc_long <- load_metadata_table("acc_long")
acc_mv   <- make_acc_mv_from_long(acc_long)

rm(acc_long)
gc()

source("code-A.R")

####################################################################################
# Analysis B: multivariate linear mixed model
# Analysis B uses acc_long.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN2")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

acc_long <- load_metadata_table("acc_long")
# acc_mv   <- make_acc_mv_from_long(acc_long)

# rm(acc_long)
gc()

source("code-B.R")

source("code-B-RE-simplex.R")

####################################################################################
# Analysis C1: pairwise method difference, Hotelling-style MANOVA
# Analysis C1 uses acc_long and creates diff_i_iv and diff_i_iv_mv.
# Analysis C2: MANOVA on pairwise differences
# Analysis C2 depends on diff_i_iv_mv, which is created in C1. Run C1 first, then C2.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN2")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

acc_long <- load_metadata_table("acc_long")

source("code-C1.R")
source("code-C2.R")

####################################################################################
# Analysis D1: PERMANOVA on multivariate accuracy
# Analysis D1 uses acc_long, then builds acc_mv2 and Y.
# Analysis D2: blocked PERMANOVA
# Analysis D2 depends on acc_mv2 and Y, which are created in D1. Run D1 first, then D2.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN2")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

acc_long <- load_metadata_table("acc_long")

source("code-D1.R")
source("code-D2.R")

####################################################################################
# Analysis E: GAM for spatial accuracy patterns
# Analysis E uses acc_long and also needs make_pairwise_diff(), as defined in code-C1.R.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN2")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

acc_long <- load_metadata_table("acc_long")

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

source("code-E.R")

####################################################################################
# Analysis F: PCA of full accuracy profile
# Analysis F uses all_wide, not acc_long.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN2")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

all_wide <- load_metadata_table("all_wide")

source("code-F.R")

####################################################################################
# Analysis G: rank-based method comparison
# Analysis G uses acc_long.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN2")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

acc_long <- load_metadata_table("acc_long")

source("code-G-revised.R")

####################################################################################
# Analysis H: error-vector MANOVA
# Analysis H uses all_wide, not acc_long.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN2")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

all_wide <- load_metadata_table("all_wide")

source("code-H.R")

####################################################################################