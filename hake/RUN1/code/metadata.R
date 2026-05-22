# file = metadata.R
# Revised metadata preparation for the Dirichlet Study hake simulation outputs.
#
# Goal:
#   Avoid the large allocation created by building all 972 experiments in memory
#   and then expanding all_wide to acc_long in one step.
#
# Default behavior when sourced:
#   1. Reads hake_factorial_design_matrix.csv.
#   2. Processes experiment CSV files in small example_id blocks.
#   3. Writes partitioned all_wide and acc_long tables to metadata_output/.
#   4. Does NOT combine the full tables in memory unless requested.
#
# Before source(), optional controls can be set in R:
#   Sys.setenv(HAKE_METADATA_CHUNK_N = "12")          # lower this if memory is tight
#   Sys.setenv(HAKE_COMBINE_IN_MEMORY = "true")       # only if you have enough RAM
#   Sys.setenv(HAKE_KEEP_ONLY_ANALYSIS_COLUMNS = "true")
#   Sys.setenv(HAKE_WRITE_CSV = "true")
#   Sys.setenv(HAKE_WRITE_PARQUET = "auto")           # auto/true/false; uses arrow if installed
#   Sys.setenv(HAKE_OVERWRITE_METADATA = "true")
#
# Run these options in R to set up metadata.R under approach 1
#   Sys.setenv(HAKE_METADATA_CHUNK_N = "12")          
#   Sys.setenv(HAKE_COMBINE_IN_MEMORY = "false")       
#   Sys.setenv(HAKE_KEEP_ONLY_ANALYSIS_COLUMNS = "true")
#   Sys.setenv(HAKE_WRITE_PARQUET = "false")
#   Sys.setenv(HAKE_WRITE_CSV = "false")  
#
# After source("metadata.R"), useful objects/functions include:
#   metadata_result
#   metadata_manifest
#   load_metadata_table("all_wide")
#   load_metadata_table("acc_long")
#   open_metadata_dataset("all_wide")     # requires arrow and parquet output
#   open_metadata_dataset("acc_long")     # requires arrow and parquet output
#   make_acc_mv_from_long(acc_long)
#
# For analyses A-H:
#   Analyses F and H use all_wide.
#   Analyses B, C, D, E, and G use acc_long.
#   Analysis A needs acc_mv, which can be built with make_acc_mv_from_long(acc_long).

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

# -----------------------------
# User-adjustable configuration
# -----------------------------
truthy <- function(x) tolower(as.character(x)) %in% c("true", "t", "yes", "y", "1")

PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
OUT_DIR <- Sys.getenv("HAKE_METADATA_OUTDIR", unset = file.path(PROJECT_DIR, "metadata_output"))
CHUNK_N <- as.integer(Sys.getenv("HAKE_METADATA_CHUNK_N", unset = "24"))
if (is.na(CHUNK_N) || CHUNK_N < 1L) CHUNK_N <- 24L

COMBINE_IN_MEMORY <- truthy(Sys.getenv("HAKE_COMBINE_IN_MEMORY", unset = "false"))
KEEP_ONLY_ANALYSIS_COLUMNS <- truthy(Sys.getenv("HAKE_KEEP_ONLY_ANALYSIS_COLUMNS", unset = "true"))
WRITE_CSV <- truthy(Sys.getenv("HAKE_WRITE_CSV", unset = "true"))
OVERWRITE_METADATA <- truthy(Sys.getenv("HAKE_OVERWRITE_METADATA", unset = "true"))
WRITE_PARQUET_SETTING <- tolower(Sys.getenv("HAKE_WRITE_PARQUET", unset = "auto"))

ACCURACY_COL_REGEX    <- "^(rmse|L1_norm|Linf_norm)_(i|ii|iii|iv|v)$"
ACCURACY_NAME_PATTERN <- "^(rmse|L1_norm|Linf_norm)_(i|ii|iii|iv|v)$"
PHAT_COL_REGEX        <- "^p_hat_(i|ii|iii|iv|v)_([123])$"
METHOD_LEVELS <- c("i", "ii", "iii", "iv", "v")
METRIC_LEVELS <- c("rmse", "L1_norm", "Linf_norm")

# -----------------------------
# General utilities
# -----------------------------
ensure_dir <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  invisible(x)
}

find_first_existing <- function(candidates, label) {
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Could not find ", label, ". Tried: ", paste(candidates, collapse = "; "))
  }
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

rds_part_files <- function(table_name, out_dir = OUT_DIR) {
  table_name <- match.arg(table_name, c("all_wide", "acc_long"))
  part_dir <- file.path(out_dir, paste0(table_name, "_parts"))
  files <- list.files(part_dir, pattern = paste0("^", table_name, "_part_[0-9]+\\.rds$"),
                      full.names = TRUE)
  sort(files)
}

format_gb <- function(bytes) sprintf("%.2f GB", bytes / 1024^3)

# -----------------------------
# Data readers and transforms
# -----------------------------
read_design <- function(design_file) {
  design <- read_csv(design_file, show_col_types = FALSE)

  if (!"example_id" %in% names(design)) {
    stop("The design matrix must contain an example_id column.")
  }

  factor_cols <- intersect(
    c("design_block", "dist_code", "G", "theta_true", "theta_CV", "sigma",
      "mean_nsamp", "ln_sd", "nb_size", "N_range_label"),
    names(design)
  )

  design %>%
    mutate(example_id = as.integer(example_id)) %>%
    mutate(across(all_of(factor_cols), as.factor))
}

thin_experiment_columns <- function(dat) {
  if (!KEEP_ONLY_ANALYSIS_COLUMNS) return(dat)

  required_base <- c("mesh_id", "sim_id", "p1", "p2", "p3")
  missing_base <- setdiff(required_base, names(dat))
  if (length(missing_base)) {
    stop("Experiment file is missing required columns: ", paste(missing_base, collapse = ", "))
  }

  accuracy_cols <- grep(ACCURACY_COL_REGEX, names(dat), value = TRUE)
  phat_cols <- grep(PHAT_COL_REGEX, names(dat), value = TRUE)
  keep_cols <- unique(c(required_base, accuracy_cols, phat_cols))
  dat %>% select(all_of(keep_cols))
}

read_one_experiment <- function(example_id, base_dir = PROJECT_DIR) {
  ex_stem <- paste0("hake_ex", example_id)
  ex_dir  <- file.path(base_dir, paste0("ex", example_id))
  f <- file.path(ex_dir, paste0(ex_stem, ".csv"))

  if (!file.exists(f)) {
    stop("Missing experiment CSV for example_id=", example_id, ": ", f)
  }

  read_csv(f, show_col_types = FALSE) %>%
    thin_experiment_columns() %>%
    mutate(example_id = as.integer(example_id))
}

make_ids <- function(dat) {
  dat %>%
    mutate(
      dataset_id = paste(example_id, mesh_id, sim_id, sep = "_"),
      simplex_id = paste(example_id, mesh_id, sep = "_")
    )
}

make_acc_long_from_wide <- function(all_wide_chunk) {
  accuracy_cols <- grep(ACCURACY_COL_REGEX, names(all_wide_chunk), value = TRUE)
  
  if (!length(accuracy_cols)) {
    stop("No accuracy columns matched ", ACCURACY_COL_REGEX)
  }
  
  if (length(accuracy_cols) != length(METRIC_LEVELS) * length(METHOD_LEVELS)) {
    warning(
      "Expected ", length(METRIC_LEVELS) * length(METHOD_LEVELS),
      " accuracy columns, but found ", length(accuracy_cols), "."
    )
  }
  
  out <- all_wide_chunk %>%
    pivot_longer(
      cols = all_of(accuracy_cols),
      names_to = c("metric", "method"),
      names_pattern = ACCURACY_NAME_PATTERN,
      values_to = "accuracy"
    ) %>%
    mutate(
      method = factor(method, levels = METHOD_LEVELS),
      metric = factor(metric, levels = METRIC_LEVELS)
    )
  
  method_counts <- table(out$method, useNA = "ifany")
  print(method_counts)
  
  if (any(method_counts[METHOD_LEVELS] == 0)) {
    stop("One or more method levels are missing after pivot_longer(). Check ACCURACY_NAME_PATTERN.")
  }
  
  out
}

make_acc_mv_from_long <- function(acc_long) {
  acc_long %>%
    group_by(example_id, mesh_id, sim_id, method,
             p1, p2, p3, design_block, G, theta_true,
             theta_CV, sigma, mean_nsamp, ln_sd, nb_size,
             N_range_label) %>%
    summarize(
      rmse = accuracy[metric == "rmse"][1],
      L1   = accuracy[metric == "L1_norm"][1],
      Linf = accuracy[metric == "Linf_norm"][1],
      .groups = "drop"
    )
}

# -----------------------------
# Disk writers
# -----------------------------
write_one_chunk <- function(chunk_id, ids, design, out_dir = OUT_DIR) {
  message("Processing chunk ", chunk_id, " with example_id range ",
          min(ids), "-", max(ids), "; n_examples=", length(ids))

  all_wide_chunk <- purrr::map_dfr(ids, read_one_experiment) %>%
    left_join(design, by = "example_id") %>%
    make_ids()

  acc_long_chunk <- make_acc_long_from_wide(all_wide_chunk)

  wide_rds_dir <- ensure_dir(file.path(out_dir, "all_wide_parts"))
  long_rds_dir <- ensure_dir(file.path(out_dir, "acc_long_parts"))
  wide_csv_dir <- ensure_dir(file.path(out_dir, "all_wide_csv_parts"))
  long_csv_dir <- ensure_dir(file.path(out_dir, "acc_long_csv_parts"))

  wide_rds <- file.path(wide_rds_dir, sprintf("all_wide_part_%03d.rds", chunk_id))
  long_rds <- file.path(long_rds_dir, sprintf("acc_long_part_%03d.rds", chunk_id))

  saveRDS(all_wide_chunk, wide_rds, compress = "gzip")
  saveRDS(acc_long_chunk, long_rds, compress = "gzip")

  wide_csv <- NA_character_
  long_csv <- NA_character_
  if (WRITE_CSV) {
    wide_csv <- file.path(wide_csv_dir, sprintf("all_wide_part_%03d.csv.gz", chunk_id))
    long_csv <- file.path(long_csv_dir, sprintf("acc_long_part_%03d.csv.gz", chunk_id))
    write_csv(all_wide_chunk, wide_csv)
    write_csv(acc_long_chunk, long_csv)
  }

  parquet_wide <- NA_character_
  parquet_long <- NA_character_
  use_arrow <- requireNamespace("arrow", quietly = TRUE) &&
    WRITE_PARQUET_SETTING %in% c("auto", "true", "yes", "1")

  if (use_arrow) {
    parquet_wide_dir <- ensure_dir(file.path(out_dir, "parquet", "all_wide"))
    parquet_long_dir <- ensure_dir(file.path(out_dir, "parquet", "acc_long"))
    parquet_wide <- file.path(parquet_wide_dir, sprintf("part_%03d.parquet", chunk_id))
    parquet_long <- file.path(parquet_long_dir, sprintf("part_%03d.parquet", chunk_id))
    arrow::write_parquet(all_wide_chunk, parquet_wide)
    arrow::write_parquet(acc_long_chunk, parquet_long)
  }

  out <- tibble(
    chunk_id = chunk_id,
    n_examples = length(ids),
    example_id_min = min(ids),
    example_id_max = max(ids),
    nrow_all_wide = nrow(all_wide_chunk),
    nrow_acc_long = nrow(acc_long_chunk),
    all_wide_object_gb = as.numeric(object.size(all_wide_chunk)) / 1024^3,
    acc_long_object_gb = as.numeric(object.size(acc_long_chunk)) / 1024^3,
    all_wide_rds = wide_rds,
    acc_long_rds = long_rds,
    all_wide_csv = wide_csv,
    acc_long_csv = long_csv,
    all_wide_parquet = parquet_wide,
    acc_long_parquet = parquet_long
  )

  rm(all_wide_chunk, acc_long_chunk)
  invisible(gc())
  out
}

# -----------------------------
# Loaders and combination helpers
# -----------------------------
load_metadata_table <- function(table_name = c("all_wide", "acc_long"),
                                out_dir = OUT_DIR,
                                part_ids = NULL) {
  table_name <- match.arg(table_name)
  files <- rds_part_files(table_name, out_dir)
  if (!length(files)) stop("No RDS partitions found for ", table_name, " in ", out_dir)

  if (!is.null(part_ids)) {
    stopifnot(all(part_ids >= 1), all(part_ids <= length(files)))
    files <- files[part_ids]
  }

  total_rds_bytes <- sum(file.info(files)$size, na.rm = TRUE)
  message("Reading ", length(files), " ", table_name, " partition(s); compressed input size = ",
          format_gb(total_rds_bytes), ".")

  pieces <- lapply(files, readRDS)
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE))
  } else {
    dplyr::bind_rows(pieces)
  }
}

open_metadata_dataset <- function(table_name = c("all_wide", "acc_long"), out_dir = OUT_DIR) {
  table_name <- match.arg(table_name)
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("The arrow package is not installed. Install it or use load_metadata_table().")
  }
  ds_dir <- file.path(out_dir, "parquet", table_name)
  if (!dir.exists(ds_dir) || !length(list.files(ds_dir, pattern = "\\.parquet$"))) {
    stop("No parquet dataset found at ", ds_dir,
         ". Re-run with Sys.setenv(HAKE_WRITE_PARQUET = 'true') and arrow installed.")
  }
  arrow::open_dataset(ds_dir)
}

combine_partitioned_metadata <- function(table_name = c("all_wide", "acc_long"),
                                         out_dir = OUT_DIR,
                                         combined_file = NULL,
                                         max_compressed_input_gb = 12) {
  table_name <- match.arg(table_name)
  files <- rds_part_files(table_name, out_dir)
  if (!length(files)) stop("No RDS partitions found for ", table_name)

  input_gb <- sum(file.info(files)$size, na.rm = TRUE) / 1024^3
  if (input_gb > max_compressed_input_gb) {
    warning("The compressed RDS partitions sum to ", sprintf("%.2f", input_gb),
            " GB. The uncompressed object can be much larger. ",
            "Consider using open_metadata_dataset() with arrow or load selected partitions.")
  }

  combined <- load_metadata_table(table_name, out_dir)
  if (is.null(combined_file)) {
    combined_file <- file.path(out_dir, paste0(table_name, "_combined.rds"))
  }
  saveRDS(combined, combined_file, compress = "gzip")
  message("Wrote combined ", table_name, " table to: ", combined_file)
  invisible(combined_file)
}

# -----------------------------
# Main runner
# -----------------------------
run_metadata_prep <- function(project_dir = PROJECT_DIR,
                              out_dir = OUT_DIR,
                              chunk_n = CHUNK_N,
                              combine_in_memory = COMBINE_IN_MEMORY,
                              overwrite = OVERWRITE_METADATA) {

  design_file <- find_first_existing(
    c(Sys.getenv("HAKE_DESIGN_FILE", unset = ""),
      file.path(project_dir, "hake_factorial_design_matrix.csv"),
      file.path(project_dir, "hake_factorial_design_matrix(1).csv")),
    label = "design matrix CSV"
  )

  if (overwrite && dir.exists(out_dir)) {
    message("Removing existing metadata output directory: ", out_dir)
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }
  ensure_dir(out_dir)

  design <- read_design(design_file)
  ids <- sort(unique(design$example_id))
  chunks <- split(ids, ceiling(seq_along(ids) / chunk_n))

  message("Design file: ", design_file)
  message("Number of experiments in design: ", length(ids))
  message("Chunk size: ", chunk_n, " experiment(s)")
  message("Output directory: ", out_dir)
  message("Keep only analysis columns: ", KEEP_ONLY_ANALYSIS_COLUMNS)
  message("Write CSV partitions: ", WRITE_CSV)
  message("Write parquet partitions: ", requireNamespace("arrow", quietly = TRUE) &&
            WRITE_PARQUET_SETTING %in% c("auto", "true", "yes", "1"))

  manifest <- purrr::imap_dfr(chunks, ~ write_one_chunk(as.integer(.y), .x, design, out_dir))
  manifest_file <- file.path(out_dir, "metadata_manifest.csv")
  write_csv(manifest, manifest_file)

  message("Finished partitioned metadata preparation.")
  message("Manifest: ", manifest_file)
  message("Total all_wide rows: ", sum(manifest$nrow_all_wide))
  message("Total acc_long rows: ", sum(manifest$nrow_acc_long))
  message("Peak all_wide chunk object size: ", sprintf("%.3f GB", max(manifest$all_wide_object_gb)))
  message("Peak acc_long chunk object size: ", sprintf("%.3f GB", max(manifest$acc_long_object_gb)))

  if (combine_in_memory) {
    message("Combining partitions in memory because HAKE_COMBINE_IN_MEMORY=true.")
    assign("all_wide", load_metadata_table("all_wide", out_dir), envir = .GlobalEnv)
    assign("acc_long", load_metadata_table("acc_long", out_dir), envir = .GlobalEnv)
    assign("acc_mv", make_acc_mv_from_long(get("acc_long", envir = .GlobalEnv)), envir = .GlobalEnv)
  } else {
    message("Full tables were not combined in memory. To load later, run:")
    message("  all_wide <- load_metadata_table('all_wide')")
    message("  acc_long <- load_metadata_table('acc_long')")
    message("  acc_mv <- make_acc_mv_from_long(acc_long)")
  }

  invisible(list(
    project_dir = project_dir,
    out_dir = out_dir,
    design_file = design_file,
    manifest_file = manifest_file,
    manifest = manifest,
    combine_in_memory = combine_in_memory
  ))
}

# Run automatically when sourced, unless disabled.
RUN_ON_SOURCE <- truthy(Sys.getenv("HAKE_RUN_ON_SOURCE", unset = "true"))
if (RUN_ON_SOURCE) {
  metadata_result <- run_metadata_prep()
  metadata_manifest <- metadata_result$manifest
}
