# file = run_hake_experiments.R
#
# --- Example Usage ---
# To run the process for hake_ex1.inp through hake_ex3.inp, simply call:
# source("run_hake_experiments.R")
# run_hake_experiments(3)

## Clear environment
rm(list = ls())

run_hake_experiments <- function(n) {

  if (!is.numeric(n) || length(n) != 1 || n < 1 || n %% 1 != 0) {
    stop("Error: Please provide a positive integer for 'n'.")
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
      "hake.lst" = file.path(ex_dir, paste0(ex_stem, ".lst")),
      "hake.csv" = file.path(ex_dir, paste0(ex_stem, ".csv")),
      "hake_emmeans_summary.csv" = file.path(ex_dir, paste0(ex_stem, "_emmeans_summary.csv")),
      "hake_L1_norm_emmeans_summary.csv" = file.path(ex_dir, paste0(ex_stem, "_L1_norm_emmeans_summary.csv")),
      "hake_Linf_norm_emmeans_summary.csv" = file.path(ex_dir, paste0(ex_stem, "_Linf_norm_emmeans_summary.csv")),
      "simplex_samples.csv" = file.path(ex_dir, paste0(ex_stem, "_simplex_samples.csv")),
      "rmse_i_hake.png" = file.path(ex_dir, paste0("rmse_i_", ex_stem, ".png")),
      "rmse_ii_hake.png" = file.path(ex_dir, paste0("rmse_ii_", ex_stem, ".png")),
      "rmse_iii_hake.png" = file.path(ex_dir, paste0("rmse_iii_", ex_stem, ".png")),
      "rmse_iv_hake.png" = file.path(ex_dir, paste0("rmse_iv_", ex_stem, ".png")),
      "rmse_v_hake.png" = file.path(ex_dir, paste0("rmse_v_", ex_stem, ".png")),
      "rmse_boxplot_hake.png" = file.path(ex_dir, paste0("rmse_boxplot_", ex_stem, ".png")),
      "rmse_emmeans_hake.png" = file.path(ex_dir, paste0("rmse_emmeans_", ex_stem, ".png")),
      "L1_norm_i_hake.png" = file.path(ex_dir, paste0("L1_norm_i_", ex_stem, ".png")),
      "L1_norm_ii_hake.png" = file.path(ex_dir, paste0("L1_norm_ii_", ex_stem, ".png")),
      "L1_norm_iii_hake.png" = file.path(ex_dir, paste0("L1_norm_iii_", ex_stem, ".png")),
      "L1_norm_iv_hake.png" = file.path(ex_dir, paste0("L1_norm_iv_", ex_stem, ".png")),
      "L1_norm_v_hake.png" = file.path(ex_dir, paste0("L1_norm_v_", ex_stem, ".png")),
      "L1_norm_boxplot_hake.png" = file.path(ex_dir, paste0("L1_norm_boxplot_", ex_stem, ".png")),
      "L1_norm_emmeans_hake.png" = file.path(ex_dir, paste0("L1_norm_emmeans_", ex_stem, ".png")),
      "Linf_norm_i_hake.png" = file.path(ex_dir, paste0("Linf_norm_i_", ex_stem, ".png")),
      "Linf_norm_ii_hake.png" = file.path(ex_dir, paste0("Linf_norm_ii_", ex_stem, ".png")),
      "Linf_norm_iii_hake.png" = file.path(ex_dir, paste0("Linf_norm_iii_", ex_stem, ".png")),
      "Linf_norm_iv_hake.png" = file.path(ex_dir, paste0("Linf_norm_iv_", ex_stem, ".png")),
      "Linf_norm_v_hake.png" = file.path(ex_dir, paste0("Linf_norm_v_", ex_stem, ".png")),
      "Linf_norm_boxplot_hake.png" = file.path(ex_dir, paste0("Linf_norm_boxplot_", ex_stem, ".png")),
      "Linf_norm_emmeans_hake.png" = file.path(ex_dir, paste0("Linf_norm_emmeans_", ex_stem, ".png"))
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
    if (!file.exists("hake.R")) {
      stop("Error: 'hake.R' not found in the working directory.")
    }

    ensure_dir(ex_dir)

    safe_copy(ex_inp_name, file.path(ex_dir, ex_inp_name))

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

    ok_swap_in <- file.copy(from = ex_inp_name, to = "hake.inp", overwrite = TRUE)
    if (!ok_swap_in) {
      if (backup_exists && file.exists("hake_backup.inp")) {
        file.rename("hake_backup.inp", "hake.inp")
      }
      stop(sprintf("Error: Could not copy '%s' to 'hake.inp'.", ex_inp_name))
    }

    stale_outputs <- c(
      "hake.lst", "hake.csv", "hake_emmeans_summary.csv", "hake_L1_norm_emmeans_summary.csv",
      "hake_Linf_norm_emmeans_summary.csv", "simplex_samples.csv", "rmse_i_hake.png",
      "rmse_ii_hake.png", "rmse_iii_hake.png", "rmse_iv_hake.png", "rmse_v_hake.png",
      "rmse_boxplot_hake.png", "rmse_emmeans_hake.png", "L1_norm_i_hake.png",
      "L1_norm_ii_hake.png", "L1_norm_iii_hake.png", "L1_norm_iv_hake.png", "L1_norm_v_hake.png",
      "L1_norm_boxplot_hake.png", "L1_norm_emmeans_hake.png", "Linf_norm_i_hake.png",
      "Linf_norm_ii_hake.png", "Linf_norm_iii_hake.png", "Linf_norm_iv_hake.png", "Linf_norm_v_hake.png",
      "Linf_norm_boxplot_hake.png", "Linf_norm_emmeans_hake.png"
    )
    invisible(lapply(stale_outputs[file.exists(stale_outputs)], safe_remove_with_retry))

    run_error <- NULL

    tryCatch(
      {
        source("hake.R", echo = FALSE, print.eval = TRUE, local = new.env(parent = globalenv()))
      },
      error = function(e) {
        run_error <<- e
      },
      finally = {
        restore_sinks(output_baseline = output_sink_baseline,
                      message_baseline = message_sink_baseline)
        invisible(try(future::plan(future::sequential), silent = TRUE))
        invisible(try(gc(), silent = TRUE))
        if (file.exists("hake.inp")) {
          safe_remove_with_retry("hake.inp")
        }
        if (backup_exists && file.exists("hake_backup.inp")) {
          safe_rename("hake_backup.inp", "hake.inp")
        }
      }
    )

    if (!is.null(run_error)) {
      stop(sprintf("Error while running hake.R for '%s': %s", ex_inp_name, conditionMessage(run_error)))
    }

    output_map <- build_output_map(ex_stem, ex_dir)

    restore_sinks(output_baseline = output_sink_baseline,
                  message_baseline = message_sink_baseline)
    invisible(try(future::plan(future::sequential), silent = TRUE))
    invisible(try(gc(), silent = TRUE))

    if (file.exists("hake.lst")) {
      ex_lst_path <- output_map[["hake.lst"]]
      ok_copy_lst <- safe_copy("hake.lst", ex_lst_path)
      if (ok_copy_lst) {
        stamp_lst_filename(ex_lst_path)
      }
      safe_remove_with_retry("hake.lst")
    } else {
      warning(sprintf("Expected log file '%s' was not created for %s.", "hake.lst", ex_inp_name))
    }

    other_outputs <- setdiff(names(output_map), "hake.lst")
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

  if (file.exists("hake.lst")) {
    safe_remove_with_retry("hake.lst", tries = 40L, wait_sec = 0.5)
  }

  if (file.exists("hake.lst")) {
    warning("The concatenated 'hake.lst' file still exists in the working directory after final cleanup.")
  }

  cat(sprintf("\n=== All %d examples processed successfully! ===\n", n))
}
