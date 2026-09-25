# file = run_experiments.R
#
# --- Example Usage ---
# To run faces_ex4.inp through faces_ex10.inp (inclusive), call:
# source("run_experiments.R")
# run_experiments(4, 10)

## Clear environment
rm(list = ls())

run_experiments <- function(start_experiment, end_experiment) {

  is_single_integer <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x %% 1 == 0
  }

  if (!is_single_integer(start_experiment) || start_experiment < 1) {
    stop("Error: 'start_experiment' must be a positive integer.")
  }
  if (!is_single_integer(end_experiment) || end_experiment < 1) {
    stop("Error: 'end_experiment' must be a positive integer.")
  }
  if (start_experiment >= end_experiment) {
    stop("Error: 'start_experiment' must be less than 'end_experiment'.")
  }

  experiment_ids <- seq.int(from = as.integer(start_experiment),
                            to = as.integer(end_experiment))
  n_experiments <- length(experiment_ids)

  if (!file.exists("faces.R")) {
    stop("Error: 'faces.R' not found in the working directory.")
  }

  experiment_input_files <- paste0("faces_ex", experiment_ids, ".inp")
  missing_input_files <- experiment_input_files[!file.exists(experiment_input_files)]
  if (length(missing_input_files) > 0L) {
    stop(sprintf(
      "Error: The following requested input file(s) were not found: %s",
      paste(missing_input_files, collapse = ", ")
    ))
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

  safe_rename <- function(from, to) {
    if (!file.exists(from)) {
      return(invisible(FALSE))
    }
    if (file.exists(to)) {
      file.remove(to)
    }
    ok <- suppressWarnings(file.rename(from = from, to = to))
    if (!ok) {
      warning(sprintf("Could not rename '%s' to '%s'.", from, to))
    }
    invisible(ok)
  }

  safe_copy <- function(from, to) {
    if (!file.exists(from)) {
      return(invisible(FALSE))
    }
    ensure_dir(dirname(to))
    if (file.exists(to)) {
      file.remove(to)
    }
    ok <- suppressWarnings(file.copy(from = from, to = to, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE))
    if (!ok) {
      warning(sprintf("Could not copy '%s' to '%s'.", from, to))
    }
    invisible(ok)
  }

  safe_move <- function(from, to) {
    if (!file.exists(from)) {
      return(invisible(FALSE))
    }
    ensure_dir(dirname(to))
    if (file.exists(to)) {
      file.remove(to)
    }

    ok <- suppressWarnings(file.rename(from = from, to = to))
    if (ok) {
      return(invisible(TRUE))
    }

    copied <- suppressWarnings(file.copy(from = from, to = to, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE))
    if (!copied || !file.exists(to)) {
      warning(sprintf("Could not move '%s' to '%s'.", from, to))
      return(invisible(FALSE))
    }

    removed <- suppressWarnings(file.remove(from))
    if (!removed && file.exists(from)) {
      warning(sprintf("Copied '%s' to '%s' but could not remove the source file.", from, to))
      return(invisible(FALSE))
    }

    invisible(TRUE)
  }

  restore_sinks <- function(output_baseline = 0L, message_baseline = 2L) {
    while (sink.number(type = "output") > output_baseline) {
      sink(type = "output")
    }
    while (sink.number(type = "message") > message_baseline) {
      sink(type = "message")
    }
    invisible(TRUE)
  }

  output_sink_baseline <- sink.number(type = "output")
  message_sink_baseline <- sink.number(type = "message")

  safe_remove_with_retry <- function(path, tries = 20L, wait_sec = 0.25) {
    if (!file.exists(path)) {
      return(invisible(TRUE))
    }

    for (attempt in seq_len(tries)) {
      restore_sinks(output_baseline = output_sink_baseline,
                    message_baseline = message_sink_baseline)
      invisible(try(gc(), silent = TRUE))

      removed <- suppressWarnings(file.remove(path))
      if (isTRUE(removed) || !file.exists(path)) {
        return(invisible(TRUE))
      }

      Sys.sleep(wait_sec)
    }

    warning(sprintf("Could not remove '%s' after %d attempts.", path, tries))
    invisible(FALSE)
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

  build_output_map <- function(ex_stem, ex_dir) {
    c(
      "faces.lst" = file.path(ex_dir, paste0(ex_stem, ".lst")),
      "faces.csv" = file.path(ex_dir, paste0(ex_stem, ".csv")),
      "faces_emmeans_summary.csv" = file.path(ex_dir, paste0(ex_stem, "_emmeans_summary.csv")),
      "faces_L1_norm_emmeans_summary.csv" = file.path(ex_dir, paste0(ex_stem, "_L1_norm_emmeans_summary.csv")),
      "faces_Linf_norm_emmeans_summary.csv" = file.path(ex_dir, paste0(ex_stem, "_Linf_norm_emmeans_summary.csv")),
      "simplex_samples.csv" = file.path(ex_dir, paste0(ex_stem, "_simplex_samples.csv")),
      "rmse_i_faces.png" = file.path(ex_dir, paste0("rmse_i_", ex_stem, ".png")),
      "rmse_ii_faces.png" = file.path(ex_dir, paste0("rmse_ii_", ex_stem, ".png")),
      "rmse_iii_faces.png" = file.path(ex_dir, paste0("rmse_iii_", ex_stem, ".png")),
      "rmse_iv_faces.png" = file.path(ex_dir, paste0("rmse_iv_", ex_stem, ".png")),
      "rmse_v_faces.png" = file.path(ex_dir, paste0("rmse_v_", ex_stem, ".png")),
      "rmse_boxplot_faces.png" = file.path(ex_dir, paste0("rmse_boxplot_", ex_stem, ".png")),
      "rmse_emmeans_faces.png" = file.path(ex_dir, paste0("rmse_emmeans_", ex_stem, ".png")),
      "L1_norm_i_faces.png" = file.path(ex_dir, paste0("L1_norm_i_", ex_stem, ".png")),
      "L1_norm_ii_faces.png" = file.path(ex_dir, paste0("L1_norm_ii_", ex_stem, ".png")),
      "L1_norm_iii_faces.png" = file.path(ex_dir, paste0("L1_norm_iii_", ex_stem, ".png")),
      "L1_norm_iv_faces.png" = file.path(ex_dir, paste0("L1_norm_iv_", ex_stem, ".png")),
      "L1_norm_v_faces.png" = file.path(ex_dir, paste0("L1_norm_v_", ex_stem, ".png")),
      "L1_norm_boxplot_faces.png" = file.path(ex_dir, paste0("L1_norm_boxplot_", ex_stem, ".png")),
      "L1_norm_emmeans_faces.png" = file.path(ex_dir, paste0("L1_norm_emmeans_", ex_stem, ".png")),
      "Linf_norm_i_faces.png" = file.path(ex_dir, paste0("Linf_norm_i_", ex_stem, ".png")),
      "Linf_norm_ii_faces.png" = file.path(ex_dir, paste0("Linf_norm_ii_", ex_stem, ".png")),
      "Linf_norm_iii_faces.png" = file.path(ex_dir, paste0("Linf_norm_iii_", ex_stem, ".png")),
      "Linf_norm_iv_faces.png" = file.path(ex_dir, paste0("Linf_norm_iv_", ex_stem, ".png")),
      "Linf_norm_v_faces.png" = file.path(ex_dir, paste0("Linf_norm_v_", ex_stem, ".png")),
      "Linf_norm_boxplot_faces.png" = file.path(ex_dir, paste0("Linf_norm_boxplot_", ex_stem, ".png")),
      "Linf_norm_emmeans_faces.png" = file.path(ex_dir, paste0("Linf_norm_emmeans_", ex_stem, ".png"))
    )
  }

  for (run_number in seq_along(experiment_ids)) {
    i <- experiment_ids[[run_number]]
    ex_inp_name <- paste0("faces_ex", i, ".inp")
    ex_stem <- tools::file_path_sans_ext(ex_inp_name)
    ex_dir <- paste0("ex", i)

    cat(sprintf("\n=== Processing %s (%d of %d) ===\n",
                ex_inp_name, run_number, n_experiments))

    ensure_dir(ex_dir)

    safe_copy(ex_inp_name, file.path(ex_dir, ex_inp_name))

    backup_exists <- FALSE
    if (file.exists("faces.inp")) {
      ok_backup <- file.copy(from = "faces.inp", to = "faces_backup.inp", overwrite = TRUE)
      if (!ok_backup) {
        stop("Error: Could not create backup file 'faces_backup.inp'.")
      }
      backup_exists <- TRUE
    } else {
      warning("Original 'faces.inp' not found. Skipping backup.")
    }

    ok_swap_in <- file.copy(from = ex_inp_name, to = "faces.inp", overwrite = TRUE)
    if (!ok_swap_in) {
      if (backup_exists && file.exists("faces_backup.inp")) {
        file.rename("faces_backup.inp", "faces.inp")
      }
      stop(sprintf("Error: Could not copy '%s' to 'faces.inp'.", ex_inp_name))
    }

    stale_outputs <- c(
      "faces.lst", "faces.csv", "faces_emmeans_summary.csv", "faces_L1_norm_emmeans_summary.csv",
      "faces_Linf_norm_emmeans_summary.csv", "simplex_samples.csv", "rmse_i_faces.png",
      "rmse_ii_faces.png", "rmse_iii_faces.png", "rmse_iv_faces.png", "rmse_v_faces.png",
      "rmse_boxplot_faces.png", "rmse_emmeans_faces.png", "L1_norm_i_faces.png",
      "L1_norm_ii_faces.png", "L1_norm_iii_faces.png", "L1_norm_iv_faces.png", "L1_norm_v_faces.png",
      "L1_norm_boxplot_faces.png", "L1_norm_emmeans_faces.png", "Linf_norm_i_faces.png",
      "Linf_norm_ii_faces.png", "Linf_norm_iii_faces.png", "Linf_norm_iv_faces.png", "Linf_norm_v_faces.png",
      "Linf_norm_boxplot_faces.png", "Linf_norm_emmeans_faces.png"
    )
    invisible(lapply(stale_outputs[file.exists(stale_outputs)], safe_remove_with_retry))

    run_error <- NULL

    tryCatch(
      {
        source("faces.R", echo = FALSE, print.eval = TRUE, local = new.env(parent = globalenv()))
      },
      error = function(e) {
        run_error <<- e
      },
      finally = {
        restore_sinks(output_baseline = output_sink_baseline,
                      message_baseline = message_sink_baseline)
        invisible(try(future::plan(future::sequential), silent = TRUE))
        invisible(try(gc(), silent = TRUE))
        if (file.exists("faces.inp")) {
          safe_remove_with_retry("faces.inp")
        }
        if (backup_exists && file.exists("faces_backup.inp")) {
          safe_rename("faces_backup.inp", "faces.inp")
        }
      }
    )

    if (!is.null(run_error)) {
      stop(sprintf("Error while running faces.R for '%s': %s", ex_inp_name, conditionMessage(run_error)))
    }

    output_map <- build_output_map(ex_stem, ex_dir)

    restore_sinks(output_baseline = output_sink_baseline,
                  message_baseline = message_sink_baseline)
    invisible(try(future::plan(future::sequential), silent = TRUE))
    invisible(try(gc(), silent = TRUE))

    if (file.exists("faces.lst")) {
      ex_lst_path <- output_map[["faces.lst"]]
      ok_copy_lst <- safe_copy("faces.lst", ex_lst_path)
      if (ok_copy_lst) {
        stamp_lst_filename(ex_lst_path)
      }
      safe_remove_with_retry("faces.lst")
    } else {
      warning(sprintf("Expected log file '%s' was not created for %s.", "faces.lst", ex_inp_name))
    }

    other_outputs <- setdiff(names(output_map), "faces.lst")
    for (src in other_outputs) {
      safe_move(src, output_map[[src]])
    }

    cat(sprintf("Finished processing %s. Output files moved to %s, and a copy of %s was saved there.\n",
                ex_inp_name, ex_dir, ex_inp_name))
  }

  restore_sinks(output_baseline = output_sink_baseline,
                message_baseline = message_sink_baseline)
  invisible(try(future::plan(future::sequential), silent = TRUE))
  invisible(try(gc(), silent = TRUE))

  if (file.exists("faces.lst")) {
    safe_remove_with_retry("faces.lst", tries = 40L, wait_sec = 0.5)
  }

  if (file.exists("faces.lst")) {
    warning("The concatenated 'faces.lst' file still exists in the working directory after final cleanup.")
  }

  cat(sprintf(
    "\n=== All %d experiments (%d through %d) processed successfully! ===\n",
    n_experiments, start_experiment, end_experiment
  ))
}
