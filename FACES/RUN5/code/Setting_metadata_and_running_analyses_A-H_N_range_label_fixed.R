# file = Setting metadata and running analyses A-H.R
####################################################################################
# A key detail: after you have already built metadata_output/
# Do not let metadata.R rerun the build when you only want loader functions.
# Use HAKE_RUN_ON_SOURCE = "false" before sourcing it.

# Below is the exact R syntax to use to access the metadata in R
# after the partitioned metadata files have already been written to metadata_output/

# The key step is to source metadata.R without rerunning the metadata build.
####################################################################################
setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN3")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

# Repair legacy metadata tables that do not contain N_range_label.
# metadata.R functions and several analysis scripts expect this column.
ensure_N_range_label <- function(dat, default_label = "N50_100") {
  if (!("N_range_label" %in% names(dat))) {
    if (all(c("Nmin", "Nmax") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$Nmin, "_", dat$Nmax)
    } else if (all(c("N_min", "N_max") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$N_min, "_", dat$N_max)
    } else if ("N_range" %in% names(dat)) {
      dat$N_range_label <- as.character(dat$N_range)
    } else {
      dat$N_range_label <- default_label
    }
  }
  dat
}


# This loads the helper functions, especially:
# load_metadata_table()
# make_acc_mv_from_long()

####################################################################################
# Analysis A: MANOVA on multivariate accuracy
# Analysis A uses acc_mv, not all_wide directly. acc_mv is built from acc_long.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN3")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

# Repair legacy metadata tables that do not contain N_range_label.
# metadata.R functions and several analysis scripts expect this column.
ensure_N_range_label <- function(dat, default_label = "N50_100") {
  if (!("N_range_label" %in% names(dat))) {
    if (all(c("Nmin", "Nmax") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$Nmin, "_", dat$Nmax)
    } else if (all(c("N_min", "N_max") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$N_min, "_", dat$N_max)
    } else if ("N_range" %in% names(dat)) {
      dat$N_range_label <- as.character(dat$N_range)
    } else {
      dat$N_range_label <- default_label
    }
  }
  dat
}

acc_long <- load_metadata_table("acc_long")
acc_long <- ensure_N_range_label(acc_long)
acc_mv   <- make_acc_mv_from_long(acc_long)

rm(acc_long)
gc()

source("code-A.R")

####################################################################################
# Analysis B: multivariate linear mixed model
# Analysis B uses acc_long.

rm(list = ls())
gc()

# setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN3")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

# Repair legacy metadata tables that do not contain N_range_label.
# metadata.R functions and several analysis scripts expect this column.
ensure_N_range_label <- function(dat, default_label = "N50_100") {
  if (!("N_range_label" %in% names(dat))) {
    if (all(c("Nmin", "Nmax") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$Nmin, "_", dat$Nmax)
    } else if (all(c("N_min", "N_max") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$N_min, "_", dat$N_max)
    } else if ("N_range" %in% names(dat)) {
      dat$N_range_label <- as.character(dat$N_range)
    } else {
      dat$N_range_label <- default_label
    }
  }
  dat
}

acc_long <- load_metadata_table("acc_long")
acc_long <- ensure_N_range_label(acc_long)
# use for code-B acc_mv   <- make_acc_mv_from_long(acc_long)

# use for code-B rm(acc_long)
gc()

source("code-B-run3.R")

source("code-B-RE-simplex-run3.R")

####################################################################################
# Analysis C1: pairwise method difference, Hotelling-style MANOVA
# Analysis C1 uses acc_long and creates diff_i_iv and diff_i_iv_mv.
# Analysis C2: MANOVA on pairwise differences
# Analysis C2 depends on diff_i_iv_mv, which is created in C1. Run C1 first, then C2.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN3")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

# Repair legacy metadata tables that do not contain N_range_label.
# metadata.R functions and several analysis scripts expect this column.
ensure_N_range_label <- function(dat, default_label = "N50_100") {
  if (!("N_range_label" %in% names(dat))) {
    if (all(c("Nmin", "Nmax") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$Nmin, "_", dat$Nmax)
    } else if (all(c("N_min", "N_max") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$N_min, "_", dat$N_max)
    } else if ("N_range" %in% names(dat)) {
      dat$N_range_label <- as.character(dat$N_range)
    } else {
      dat$N_range_label <- default_label
    }
  }
  dat
}


acc_long <- load_metadata_table("acc_long")
acc_long <- ensure_N_range_label(acc_long)

source("code-C1-run3.R")
source("code-C2-run3.R")

####################################################################################
# Analysis D1: PERMANOVA on multivariate accuracy
# Analysis D1 uses acc_long, then builds acc_mv2 and Y.
# Analysis D2: blocked PERMANOVA
# Analysis D2 depends on acc_mv2 and Y, which are created in D1. Run D1 first, then D2.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN3")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

# Repair legacy metadata tables that do not contain N_range_label.
# metadata.R functions and several analysis scripts expect this column.
ensure_N_range_label <- function(dat, default_label = "N50_100") {
  if (!("N_range_label" %in% names(dat))) {
    if (all(c("Nmin", "Nmax") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$Nmin, "_", dat$Nmax)
    } else if (all(c("N_min", "N_max") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$N_min, "_", dat$N_max)
    } else if ("N_range" %in% names(dat)) {
      dat$N_range_label <- as.character(dat$N_range)
    } else {
      dat$N_range_label <- default_label
    }
  }
  dat
}


acc_long <- load_metadata_table("acc_long")
acc_long <- ensure_N_range_label(acc_long)

source("code-D1.R")
source("code-D2.R")

####################################################################################
# Analysis E: GAM for spatial accuracy patterns
# Analysis E uses acc_long and also needs make_pairwise_diff(), as defined in code-C1.R.

rm(list = ls())
gc()

# setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN3")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

# Repair legacy metadata tables that do not contain N_range_label.
# metadata.R functions and several analysis scripts expect this column.
ensure_N_range_label <- function(dat, default_label = "N50_100") {
  if (!("N_range_label" %in% names(dat))) {
    if (all(c("Nmin", "Nmax") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$Nmin, "_", dat$Nmax)
    } else if (all(c("N_min", "N_max") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$N_min, "_", dat$N_max)
    } else if ("N_range" %in% names(dat)) {
      dat$N_range_label <- as.character(dat$N_range)
    } else {
      dat$N_range_label <- default_label
    }
  }
  dat
}

acc_long <- load_metadata_table("acc_long")
acc_long <- ensure_N_range_label(acc_long)

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

source("code-E-run3.R")

####################################################################################
# Analysis F: PCA of full accuracy profile
# Analysis F uses all_wide, not acc_long.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN3")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

# Repair legacy metadata tables that do not contain N_range_label.
# metadata.R functions and several analysis scripts expect this column.
ensure_N_range_label <- function(dat, default_label = "N50_100") {
  if (!("N_range_label" %in% names(dat))) {
    if (all(c("Nmin", "Nmax") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$Nmin, "_", dat$Nmax)
    } else if (all(c("N_min", "N_max") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$N_min, "_", dat$N_max)
    } else if ("N_range" %in% names(dat)) {
      dat$N_range_label <- as.character(dat$N_range)
    } else {
      dat$N_range_label <- default_label
    }
  }
  dat
}


all_wide <- load_metadata_table("all_wide")
all_wide <- ensure_N_range_label(all_wide)

source("code-F.R")

####################################################################################
# Analysis G: rank-based method comparison
# Analysis G uses acc_long.

rm(list = ls())
gc()

# setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN3")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

# Repair legacy metadata tables that do not contain N_range_label.
# metadata.R functions and several analysis scripts expect this column.
ensure_N_range_label <- function(dat, default_label = "N50_100") {
  if (!("N_range_label" %in% names(dat))) {
    if (all(c("Nmin", "Nmax") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$Nmin, "_", dat$Nmax)
    } else if (all(c("N_min", "N_max") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$N_min, "_", dat$N_max)
    } else if ("N_range" %in% names(dat)) {
      dat$N_range_label <- as.character(dat$N_range)
    } else {
      dat$N_range_label <- default_label
    }
  }
  dat
}

acc_long <- load_metadata_table("acc_long")
acc_long <- ensure_N_range_label(acc_long)

source("code-G-run3.R")

####################################################################################
# Analysis H: error-vector MANOVA
# Analysis H uses all_wide, not acc_long.

rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN3")

Sys.setenv(HAKE_RUN_ON_SOURCE = "false")
Sys.setenv(HAKE_METADATA_OUTDIR = file.path(getwd(), "metadata_output"))

source("metadata.R")

# Repair legacy metadata tables that do not contain N_range_label.
# metadata.R functions and several analysis scripts expect this column.
ensure_N_range_label <- function(dat, default_label = "N50_100") {
  if (!("N_range_label" %in% names(dat))) {
    if (all(c("Nmin", "Nmax") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$Nmin, "_", dat$Nmax)
    } else if (all(c("N_min", "N_max") %in% names(dat))) {
      dat$N_range_label <- paste0("N", dat$N_min, "_", dat$N_max)
    } else if ("N_range" %in% names(dat)) {
      dat$N_range_label <- as.character(dat$N_range)
    } else {
      dat$N_range_label <- default_label
    }
  }
  dat
}


all_wide <- load_metadata_table("all_wide")
all_wide <- ensure_N_range_label(all_wide)

source("code-H.R")

####################################################################################