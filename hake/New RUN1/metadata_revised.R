# file = metadata_revised.R
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

# Sys.setenv(
#  HAKE_METADATA_CHUNK_N = "12",
#  HAKE_COMBINE_IN_MEMORY = "false",
#  HAKE_KEEP_ONLY_ANALYSIS_COLUMNS = "true",
#  HAKE_WRITE_PARQUET = "true",
#  HAKE_WRITE_CSV = "false",
#  HAKE_OVERWRITE_METADATA = "true",
#  HAKE_SETUP_DISK_BACKED_OBJECTS = "true"
# )

# source("metadata_revised.R")

#
# After source("metadata.R"), useful objects/functions include:
#   metadata_result
#   metadata_manifest
#   acc_long                         # disk-backed Dataset when parquet is available
#   all_wide                         # disk-backed Arrow Dataset when parquet is available
#   load_metadata_table("acc_long", part_ids = 1:2)  # load only selected RDS parts
#   open_metadata_dataset("acc_long")                 # disk-backed full table
#   collect_metadata_table("acc_long", filter_expr = example_id <= 24)
#   make_acc_mv_partitioned()        # writes acc_mv parts without loading full acc_long
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
SETUP_DISK_BACKED_OBJECTS <- truthy(Sys.getenv("HAKE_SETUP_DISK_BACKED_OBJECTS", unset = "true"))

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
                                part_ids = NULL,
                                columns = NULL,
                                max_parts_in_memory = 4L,
                                force = FALSE) {
  table_name <- match.arg(table_name)
  files <- rds_part_files(table_name, out_dir)
  if (!length(files)) stop("No RDS partitions found for ", table_name, " in ", out_dir)

  if (!is.null(part_ids)) {
    stopifnot(all(part_ids >= 1), all(part_ids <= length(files)))
    files <- files[part_ids]
  }

  if (length(files) > max_parts_in_memory && !force) {
    stop(
      "Refusing to combine ", length(files), " ", table_name,
      " partitions in memory. This previously caused the allocation failure.\n",
      "Use open_metadata_dataset('", table_name, "') for the full disk-backed table,\n",
      "or load selected parts, for example part_ids = 1:4.\n",
      "Set force = TRUE only when sufficient RAM is available."
    )
  }

  total_rds_bytes <- sum(file.info(files)$size, na.rm = TRUE)
  message("Reading ", length(files), " ", table_name,
          " partition(s); compressed input size = ", format_gb(total_rds_bytes), ".")

  if (!requireNamespace("data.table", quietly = TRUE)) {
    pieces <- lapply(files, function(f) {
      x <- readRDS(f)
      if (!is.null(columns)) x <- dplyr::select(x, dplyr::any_of(columns))
      x
    })
    return(dplyr::bind_rows(pieces))
  }

  # Keep the result as a data.table. Avoid as.data.frame(), which makes another
  # full-size copy and raises peak memory substantially.
  pieces <- lapply(files, function(f) {
    x <- data.table::as.data.table(readRDS(f))
    if (!is.null(columns)) {
      keep <- intersect(columns, names(x))
      x <- x[, ..keep]
    }
    x
  })
  data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE)
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


collect_metadata_table <- function(table_name = c("all_wide", "acc_long"),
                                   out_dir = OUT_DIR,
                                   columns = NULL,
                                   filter_expr = NULL) {
  table_name <- match.arg(table_name)
  ds <- open_metadata_dataset(table_name, out_dir)
  qry <- ds
  filter_quo <- rlang::enquo(filter_expr)
  if (!rlang::quo_is_null(filter_quo)) {
    qry <- dplyr::filter(qry, !!filter_quo)
  }
  if (!is.null(columns)) {
    qry <- dplyr::select(qry, dplyr::any_of(columns))
  }
  dplyr::collect(qry)
}

make_acc_mv_partitioned <- function(out_dir = OUT_DIR,
                                    overwrite = TRUE) {
  files <- rds_part_files("acc_long", out_dir)
  if (!length(files)) stop("No acc_long RDS partitions found in ", out_dir)

  target_dir <- file.path(out_dir, "acc_mv_parts")
  if (overwrite && dir.exists(target_dir)) unlink(target_dir, recursive = TRUE, force = TRUE)
  ensure_dir(target_dir)

  out_files <- character(length(files))
  for (i in seq_along(files)) {
    message("Creating acc_mv partition ", i, " of ", length(files), ".")
    x <- readRDS(files[[i]])
    y <- make_acc_mv_from_long(x)
    out_files[[i]] <- file.path(target_dir, sprintf("acc_mv_part_%03d.rds", i))
    saveRDS(y, out_files[[i]], compress = "gzip")
    rm(x, y)
    invisible(gc())
  }
  invisible(out_files)
}

setup_metadata_objects <- function(out_dir = OUT_DIR, envir = .GlobalEnv) {
  parquet_root <- file.path(out_dir, "parquet")
  have_arrow <- requireNamespace("arrow", quietly = TRUE)

  for (nm in c("all_wide", "acc_long")) {
    ds_dir <- file.path(parquet_root, nm)
    if (have_arrow && dir.exists(ds_dir) &&
        length(list.files(ds_dir, pattern = "\\.parquet$"))) {
      assign(nm, arrow::open_dataset(ds_dir), envir = envir)
      message("Created disk-backed object `", nm, "` from ", ds_dir, ".")
    } else {
      assign(paste0(nm, "_parts"), rds_part_files(nm, out_dir), envir = envir)
      message("Parquet not available for `", nm, "`; created `", nm,
              "_parts` as the RDS partition file list.")
    }
  }
  invisible(TRUE)
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

  combined <- load_metadata_table(table_name, out_dir, force = TRUE)
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
    assign("all_wide", load_metadata_table("all_wide", out_dir, force = TRUE), envir = .GlobalEnv)
    assign("acc_long", load_metadata_table("acc_long", out_dir, force = TRUE), envir = .GlobalEnv)
    assign("acc_mv", make_acc_mv_from_long(get("acc_long", envir = .GlobalEnv)), envir = .GlobalEnv)
  } else {
    message("Full tables were not combined in memory.")
    message("Use the disk-backed objects created by setup_metadata_objects(),")
    message("or load a small subset, e.g. load_metadata_table('acc_long', part_ids = 1:4).")
    message("To create acc_mv safely, run make_acc_mv_partitioned().")
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
  if (!COMBINE_IN_MEMORY && SETUP_DISK_BACKED_OBJECTS) {
    setup_metadata_objects(metadata_result$out_dir, envir = .GlobalEnv)
  }
}
