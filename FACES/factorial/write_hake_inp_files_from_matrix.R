# file = write_hake_inp_files_from_matrix.R
# Purpose:
#   Read a factorial design matrix and write hake_ex1.inp through hake_exN.inp
#   files for use with hake.R and run_hake_examples.R.
#
# Usage:
#   source("write_hake_inp_files_from_matrix.R")
#   write_hake_inp_files_from_matrix(
#     matrix_file = "hake_factorial_design_matrix.csv",
#     output_dir = "."
#   )
#
# Notes:
#   - This script expects the design matrix columns created in
#     hake_factorial_design_matrix.csv.
#   - Each output file is written in "key = value" format expected by hake.R.
#   - Character fields p_lbound and p_ubound are preserved as comma-separated
#     vectors, e.g. "0.0, 0.0, 0.0".

write_hake_inp_files_from_matrix <- function(matrix_file = "hake_factorial_design_matrix.csv",
                                             output_dir = ".",
                                             overwrite = TRUE,
                                             verbose = TRUE) {
  if (!file.exists(matrix_file)) {
    stop(sprintf("Input design matrix file not found: '%s'", matrix_file))
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  design <- utils::read.csv(matrix_file, stringsAsFactors = FALSE, check.names = FALSE)

  required_cols <- c(
    "example_id", "inp_file",
    "K", "G", "h", "theta_true", "theta_CV",
    "Nmin", "Nmax", "nsims", "random.seed", "od_mult", "sigma",
    "dist_code", "mean_nsamp", "ln_sd", "nb_size",
    "p_lbound", "p_ubound"
  )

  missing_cols <- setdiff(required_cols, names(design))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Design matrix is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  format_scalar <- function(x) {
    if (length(x) != 1) {
      stop("Expected a scalar value when formatting an input line.")
    }

    if (is.na(x)) {
      return("")
    }

    if (is.character(x)) {
      return(trimws(x))
    }

    if (is.numeric(x)) {
      if (isTRUE(all.equal(x, round(x)))) {
        return(as.character(as.integer(round(x))))
      } else {
        return(formatC(x, digits = 15, format = "fg", flag = "#"))
      }
    }

    as.character(x)
  }

  build_inp_lines <- function(row) {
    keys <- c(
      "K", "G", "h", "theta_true", "theta_CV",
      "Nmin", "Nmax", "nsims", "random.seed", "od_mult", "sigma",
      "dist_code", "mean_nsamp", "ln_sd", "nb_size",
      "p_lbound", "p_ubound"
    )

    values <- vapply(keys, function(k) format_scalar(row[[k]]), character(1))
    paste(keys, "=", values)
  }

  written_files <- character(nrow(design))

  for (i in seq_len(nrow(design))) {
    row <- design[i, , drop = FALSE]
    out_name <- trimws(row$inp_file[[1]])

    if (!nzchar(out_name)) {
      out_name <- sprintf("hake_ex%d.inp", as.integer(row$example_id[[1]]))
    }

    out_path <- file.path(output_dir, out_name)

    if (file.exists(out_path) && !overwrite) {
      stop(sprintf("File already exists and overwrite=FALSE: '%s'", out_path))
    }

    header <- c(
      sprintf("# file = %s", basename(out_path)),
      sprintf("# example_id = %s", format_scalar(row$example_id[[1]])),
      sprintf("# design_block = %s", format_scalar(row$design_block[[1]])),
      ""
    )

    inp_lines <- build_inp_lines(row)
    writeLines(c(header, inp_lines), con = out_path, useBytes = TRUE)

    written_files[i] <- out_path

    if (isTRUE(verbose)) {
      message(sprintf("Wrote %s", out_path))
    }
  }

  summary_df <- data.frame(
    example_id = design$example_id,
    inp_file = basename(written_files),
    output_path = normalizePath(written_files, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )

  summary_file <- file.path(output_dir, "hake_inp_file_manifest.csv")
  utils::write.csv(summary_df, summary_file, row.names = FALSE)

  if (isTRUE(verbose)) {
    message(sprintf("Wrote manifest: %s", summary_file))
    message(sprintf("Finished writing %d .inp files.", nrow(design)))
  }

  invisible(summary_df)
}

# Example:
# write_hake_inp_files_from_matrix(
#   matrix_file = "hake_factorial_design_matrix.csv",
#   output_dir = "."
# )
