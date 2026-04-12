## file = run_hake_experiments.R

# --- Example Usage ---
# To run the process for hake_ex1.inp through hake_ex3.inp, simply call:
# source("run_hake_experiments.R")
# run_hake_experiments(3)

## Clear environment
rm(list = ls())

run_hake_experiments <- function(n) {

  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n %% 1 != 0) {
    stop("Error: Please provide a positive integer for 'n'.")
  }

  safe_rename <- function(from, to) {
    if (!file.exists(from)) {
      return(invisible(FALSE))
    }

    if (file.exists(to)) {
      file.remove(to)
    }

    ok <- file.rename(from = from, to = to)
    if (!ok) {
      warning(sprintf("Could not rename '%s' to '%s'.", from, to))
    }

    invisible(ok)
  }

  ensure_dir <- function(path) {
    if (!dir.exists(path)) {
      ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
      if (!ok && !dir.exists(path)) {
        stop(sprintf("Error: Could not create directory '%s'.", path))
      }
    }
    invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  safe_move_to_dir <- function(file_name, dir_path) {
    if (!file.exists(file_name)) {
      return(invisible(FALSE))
    }

    ensure_dir(dir_path)
    dest <- file.path(dir_path, basename(file_name))

    if (file.exists(dest)) {
      file.remove(dest)
    }

    ok <- file.rename(from = file_name, to = dest)
    if (!ok) {
      copied <- file.copy(from = file_name, to = dest, overwrite = TRUE)
      if (!copied) {
        warning(sprintf("Could not move '%s' to '%s'.", file_name, dest))
        return(invisible(FALSE))
      }
      file.remove(file_name)
    }

    invisible(TRUE)
  }

  stamp_lst_filename <- function(path) {
    if (!file.exists(path)) {
      return(invisible(FALSE))
    }

    header <- sprintf("file = %s", basename(path))
    existing_lines <- readLines(path, warn = FALSE)

    if (length(existing_lines) > 0 && identical(existing_lines[1], header)) {
      return(invisible(TRUE))
    }

    writeLines(c(header, existing_lines), con = path, useBytes = TRUE)
    invisible(TRUE)
  }

  build_output_map <- function(ex_stem) {
    c(
      "hake.lst" = paste0(ex_stem, ".lst"),
      "hake.csv" = paste0(ex_stem, ".csv"),
      "hake_emmeans_summary.csv" = paste0(ex_stem, "_emmeans_summary.csv"),
      "hake_L1_norm_emmeans_summary.csv" = paste0(ex_stem, "_L1_norm_emmeans_summary.csv"),
      "hake_Linf_norm_emmeans_summary.csv" = paste0(ex_stem, "_Linf_norm_emmeans_summary.csv"),
      "simplex_samples.csv" = paste0(ex_stem, "_simplex_samples.csv"),
      "rmse_i_hake.png" = paste0("rmse_i_", ex_stem, ".png"),
      "rmse_ii_hake.png" = paste0("rmse_ii_", ex_stem, ".png"),
      "rmse_iii_hake.png" = paste0("rmse_iii_", ex_stem, ".png"),
      "rmse_iv_hake.png" = paste0("rmse_iv_", ex_stem, ".png"),
      "rmse_v_hake.png" = paste0("rmse_v_", ex_stem, ".png"),
      "rmse_boxplot_hake.png" = paste0("rmse_boxplot_", ex_stem, ".png"),
      "rmse_emmeans_hake.png" = paste0("rmse_emmeans_", ex_stem, ".png"),
      "L1_norm_i_hake.png" = paste0("L1_norm_i_", ex_stem, ".png"),
      "L1_norm_ii_hake.png" = paste0("L1_norm_ii_", ex_stem, ".png"),
      "L1_norm_iii_hake.png" = paste0("L1_norm_iii_", ex_stem, ".png"),
      "L1_norm_iv_hake.png" = paste0("L1_norm_iv_", ex_stem, ".png"),
      "L1_norm_v_hake.png" = paste0("L1_norm_v_", ex_stem, ".png"),
      "L1_norm_boxplot_hake.png" = paste0("L1_norm_boxplot_", ex_stem, ".png"),
      "L1_norm_emmeans_hake.png" = paste0("L1_norm_emmeans_", ex_stem, ".png"),
      "Linf_norm_i_hake.png" = paste0("Linf_norm_i_", ex_stem, ".png"),
      "Linf_norm_ii_hake.png" = paste0("Linf_norm_ii_", ex_stem, ".png"),
      "Linf_norm_iii_hake.png" = paste0("Linf_norm_iii_", ex_stem, ".png"),
      "Linf_norm_iv_hake.png" = paste0("Linf_norm_iv_", ex_stem, ".png"),
      "Linf_norm_v_hake.png" = paste0("Linf_norm_v_", ex_stem, ".png"),
      "Linf_norm_boxplot_hake.png" = paste0("Linf_norm_boxplot_", ex_stem, ".png"),
      "Linf_norm_emmeans_hake.png" = paste0("Linf_norm_emmeans_", ex_stem, ".png")
    )
  }

  for (i in seq_len(n)) {
    ex_inp_name <- paste0("hake_ex", i, ".inp")
    ex_stem <- tools::file_path_sans_ext(ex_inp_name)
    ex_dir <- paste0("ex", i)

    cat(sprintf("\n=== Processing %s (%d of %d) ===\n", ex_inp_name, i, n))

    if (!file.exists(ex_inp_name)) {
      stop(sprintf("Error: '%s' not found. Ensure the example file is in the working directory.", ex_inp_name))
    }

    ensure_dir(ex_dir)

    backup_exists <- FALSE
    if (file.exists("hake.inp")) {
      ok_backup <- file.copy(from = "hake.inp", to = "hake_backup.inp", overwrite = TRUE)
      if (!ok_backup) {
        stop("Error: Could not create backup file 'hake_backup.inp'.")
      }
      backup_exists <- TRUE
    } else {
      warning("Original 'hake.inp' not found. Skipping backup.")
    }

    ok_swap_in <- file.rename(from = ex_inp_name, to = "hake.inp")
    if (!ok_swap_in) {
      if (backup_exists && file.exists("hake_backup.inp")) {
        file.rename("hake_backup.inp", "hake.inp")
      }
      stop(sprintf("Error: Could not rename '%s' to 'hake.inp'.", ex_inp_name))
    }

    run_error <- NULL
    tryCatch(
      {
        source("hake.R", echo = FALSE, print.eval = TRUE)
      },
      error = function(e) {
        run_error <<- e
      },
      finally = {
        if (file.exists("hake.inp")) {
          safe_rename("hake.inp", ex_inp_name)
        }
        if (backup_exists && file.exists("hake_backup.inp")) {
          safe_rename("hake_backup.inp", "hake.inp")
        }
      }
    )

    if (!is.null(run_error)) {
      stop(sprintf("Error while running hake.R for '%s': %s", ex_inp_name, conditionMessage(run_error)))
    }

    output_map <- build_output_map(ex_stem)
    for (src in names(output_map)) {
      safe_rename(src, output_map[[src]])
    }

    ex_lst <- paste0(ex_stem, ".lst")
    stamp_lst_filename(ex_lst)

    output_files <- unique(unname(output_map))
    for (file_name in output_files) {
      safe_move_to_dir(file_name, ex_dir)
    }

    cat(sprintf("Finished processing %s. Output files moved to %s.\n", ex_inp_name, ex_dir))
  }

  cat(sprintf("\n=== All %d examples processed successfully! ===\n", n))
}
