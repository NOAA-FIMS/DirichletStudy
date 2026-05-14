# file = loading_acc_mv_code-A.R

# STEP 1
rm(list = ls())
gc()

setwd("C:/Users/Jon.Brodziak/Desktop/Dirichlet Study/hake/RUN1")

library(data.table)

out_dir <- file.path(getwd(), "metadata_output")

wide_files <- list.files(
  file.path(out_dir, "all_wide_parts"),
  pattern = "^all_wide_part_[0-9]+\\.rds$",
  full.names = TRUE
)

length(wide_files)

# STEP 2
acc_mv_dir <- file.path(out_dir, "acc_mv_parts")

if (dir.exists(acc_mv_dir)) {
  unlink(acc_mv_dir, recursive = TRUE, force = TRUE)
}

dir.create(acc_mv_dir, recursive = TRUE)

method_levels <- c("i", "ii", "iii", "iv", "v")

id_cols <- c(
  "example_id", "mesh_id", "sim_id", "dataset_id", "simplex_id",
  "p1", "p2", "p3",
  "design_block", "G", "theta_true", "theta_CV", "sigma",
  "mean_nsamp", "ln_sd", "nb_size", "N_range_label"
)

manifest <- vector("list", length(wide_files))

for (k in seq_along(wide_files)) {
  
  message("Processing acc_mv partition ", k, " of ", length(wide_files))
  
  wide <- readRDS(wide_files[k])
  setDT(wide)
  
  keep_id_cols <- intersect(id_cols, names(wide))
  
  pieces <- lapply(method_levels, function(m) {
    
    rmse_col <- paste0("rmse_", m)
    L1_col   <- paste0("L1_norm_", m)
    Linf_col <- paste0("Linf_norm_", m)
    
    missing_cols <- setdiff(c(rmse_col, L1_col, Linf_col), names(wide))
    
    if (length(missing_cols) > 0) {
      stop("Missing columns in ", basename(wide_files[k]), ": ",
           paste(missing_cols, collapse = ", "))
    }
    
    data.table(
      wide[, ..keep_id_cols],
      method = factor(m, levels = method_levels),
      rmse   = wide[[rmse_col]],
      L1     = wide[[L1_col]],
      Linf   = wide[[Linf_col]]
    )
  })
  
  acc_mv_part <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  
  out_file <- file.path(acc_mv_dir, sprintf("acc_mv_part_%03d.rds", k))
  
  saveRDS(acc_mv_part, out_file, compress = "gzip")
  
  manifest[[k]] <- data.table(
    part_id = k,
    input_file = wide_files[k],
    output_file = out_file,
    nrow_acc_mv = nrow(acc_mv_part),
    object_gb = as.numeric(object.size(acc_mv_part)) / 1024^3
  )
  
  rm(wide, pieces, acc_mv_part)
  gc()
}

acc_mv_manifest <- rbindlist(manifest)

fwrite(
  acc_mv_manifest,
  file.path(out_dir, "acc_mv_manifest.csv")
)

sum(acc_mv_manifest$nrow_acc_mv)
max(acc_mv_manifest$object_gb)

# STEP 3
rm(list = setdiff(ls(), c("out_dir")))
gc()

library(data.table)

acc_mv_files <- list.files(
  file.path(out_dir, "acc_mv_parts"),
  pattern = "^acc_mv_part_[0-9]+\\.rds$",
  full.names = TRUE
)

sum(file.info(acc_mv_files)$size) / 1024^3

acc_mv <- rbindlist(
  lapply(acc_mv_files, readRDS),
  use.names = TRUE,
  fill = TRUE
)

nrow(acc_mv)
format(object.size(acc_mv), units = "GB")
gc()

