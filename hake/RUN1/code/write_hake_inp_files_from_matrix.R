# file = write_hake_inp_files_from_matrix.R
# Purpose:
#   Read a hake.R factorial design matrix and write hake_ex1.inp through hake_exN.inp
#   files for use with hake.R and run_hake_examples.R.
#
# Compatibility:
#   This version is compatible with build_hake_factorial_design_matrix_common_G.R,
#   where G is the common fifth factor for all four sampling distributions:
#     1 = Poisson
#     2 = Lognormal
#     3 = Negative binomial
#     4 = Log-uniform
#
#   The revised common-G design fixes the distribution-specific spread parameters
#   by default at:
#     lognormal ln_sd = 1.0
#     negative-binomial nb_size = 25
#     log-uniform Nmin = 25
#     log-uniform Nmax = 100
#
# Usage:
#   source("write_hake_inp_files_from_matrix.R")
#   write_hake_inp_files_from_matrix(
#     matrix_file = "hake_factorial_design_matrix.csv",
#     output_dir = "."
#   )
#
# Notes:
#   - Each output file is written in "key = value" format expected by hake.R.
#   - Character fields p_lbound and p_ubound are preserved as comma-separated
#     vectors, e.g. "0.0, 0.0, 0.0".
#   - The compatibility checks can be relaxed by setting validate_common_G = FALSE.

write_hake_inp_files_from_matrix <- function(matrix_file = "hake_factorial_design_matrix.csv",
                                             output_dir = ".",
                                             overwrite = TRUE,
                                             verbose = TRUE,
                                             validate_common_G = TRUE,
                                             expected_G_levels = c(2, 4, 8),
                                             expected_lognormal_ln_sd = 1.0,
                                             expected_nb_size = 25,
                                             expected_loguniform_Nmin = 25,
                                             expected_loguniform_Nmax = 100) {
  if (!file.exists(matrix_file)) {
    stop(sprintf("Input design matrix file not found: '%s'", matrix_file))
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  design <- utils::read.csv(matrix_file, stringsAsFactors = FALSE, check.names = FALSE)

  required_cols <- c(
    "example_id", "inp_file", "design_block", "dist_code",
    "factor1", "level1_code", "factor2", "level2_code",
    "factor3", "level3_code", "factor4", "level4_code",
    "factor5", "level5_code",
    "K", "G", "h", "theta_true", "theta_CV",
    "Nmin", "Nmax", "nsims", "random.seed", "od_mult", "sigma",
    "mean_nsamp", "ln_sd", "nb_size", "p_lbound", "p_ubound"
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

  validate_design_matrix <- function(design) {
    if (nrow(design) < 1L) {
      stop("Design matrix has no rows.")
    }

    if (anyDuplicated(design$example_id)) {
      stop("Design matrix contains duplicate example_id values.")
    }

    if (anyDuplicated(trimws(design$inp_file))) {
      stop("Design matrix contains duplicate inp_file values.")
    }

    if (any(!design$dist_code %in% c(1L, 2L, 3L, 4L))) {
      stop("dist_code must contain only values 1, 2, 3, and 4.")
    }

    numeric_cols <- c(
      "K", "G", "h", "theta_true", "theta_CV", "Nmin", "Nmax",
      "nsims", "random.seed", "od_mult", "sigma", "mean_nsamp",
      "ln_sd", "nb_size"
    )

    for (cc in numeric_cols) {
      vals <- suppressWarnings(as.numeric(design[[cc]]))
      if (anyNA(vals)) {
        stop(sprintf("Column '%s' contains non-numeric or missing values.", cc))
      }
      design[[cc]] <- vals
    }

    integer_cols <- c("K", "G", "Nmin", "Nmax", "nsims", "random.seed", "dist_code")
    for (cc in integer_cols) {
      vals <- as.numeric(design[[cc]])
      if (any(abs(vals - round(vals)) > sqrt(.Machine$double.eps))) {
        stop(sprintf("Column '%s' must contain integer-valued entries.", cc))
      }
    }

    if (any(design$K < 2)) stop("K must be at least 2.")
    if (any(design$G < 1)) stop("G must be at least 1.")
    if (any(design$Nmin < 1)) stop("Nmin must be at least 1.")
    if (any(design$Nmin >= design$Nmax)) stop("Each row must have Nmin < Nmax.")
    if (any(design$nsims < 1)) stop("nsims must be at least 1.")
    if (any(design$mean_nsamp <= 0)) stop("mean_nsamp must be positive.")
    if (any(design$sigma <= 0)) stop("sigma must be positive.")
    if (any(design$theta_true <= 0)) stop("theta_true must be positive.")
    if (any(design$theta_CV < 0)) stop("theta_CV must be nonnegative.")
    if (any(design$ln_sd <= 0)) stop("ln_sd must be positive.")
    if (any(design$nb_size <= 0)) stop("nb_size must be positive.")

    if (!isTRUE(validate_common_G)) {
      return(invisible(design))
    }

    expected_blocks <- c("Poisson", "Lognormal", "Negative binomial", "Log-uniform")
    observed_blocks <- sort(unique(design$design_block))
    if (!identical(observed_blocks, sort(expected_blocks))) {
      stop(
        "Common-G validation expected design_block values: ",
        paste(expected_blocks, collapse = ", "),
        "; observed: ", paste(observed_blocks, collapse = ", ")
      )
    }

    expected_n <- length(expected_blocks) * 3L^5L
    if (nrow(design) != expected_n) {
      stop(sprintf(
        "Common-G design should have %d rows = 4 distributions x 3^5, but found %d rows.",
        expected_n, nrow(design)
      ))
    }

    rows_by_block <- table(design$design_block)
    if (!all(rows_by_block[expected_blocks] == 3L^5L)) {
      stop("Each design_block must contain exactly 3^5 = 243 rows.")
    }

    if (!all(trimws(design$factor5) == "G")) {
      stop("For the revised common-G design, factor5 must be 'G' in every row.")
    }

    expected_G_levels_int <- as.integer(round(expected_G_levels))
    if (length(expected_G_levels_int) != 3L) {
      stop("expected_G_levels must contain exactly 3 values.")
    }

    for (blk in expected_blocks) {
      got_G <- sort(unique(as.integer(round(design$G[design$design_block == blk]))))
      if (!identical(got_G, sort(expected_G_levels_int))) {
        stop(sprintf(
          "Block '%s' does not contain the expected G levels: expected %s; observed %s.",
          blk,
          paste(sort(expected_G_levels_int), collapse = ", "),
          paste(got_G, collapse = ", ")
        ))
      }
    }

    if (!all(design$ln_sd[design$dist_code == 2L] == expected_lognormal_ln_sd)) {
      stop(sprintf(
        "For dist_code = 2, ln_sd must be constant at %s. If you changed this baseline in the design builder, pass expected_lognormal_ln_sd accordingly.",
        expected_lognormal_ln_sd
      ))
    }

    if (!all(design$nb_size[design$dist_code == 3L] == expected_nb_size)) {
      stop(sprintf(
        "For dist_code = 3, nb_size must be constant at %s. If you changed this baseline in the design builder, pass expected_nb_size accordingly.",
        expected_nb_size
      ))
    }

    if (!all(design$Nmin[design$dist_code == 4L] == expected_loguniform_Nmin)) {
      stop(sprintf(
        "For dist_code = 4, Nmin must be constant at %s. If you changed this baseline in the design builder, pass expected_loguniform_Nmin accordingly.",
        expected_loguniform_Nmin
      ))
    }

    if (!all(design$Nmax[design$dist_code == 4L] == expected_loguniform_Nmax)) {
      stop(sprintf(
        "For dist_code = 4, Nmax must be constant at %s. If you changed this baseline in the design builder, pass expected_loguniform_Nmax accordingly.",
        expected_loguniform_Nmax
      ))
    }

    invisible(design)
  }

  design <- validate_design_matrix(design)

  format_scalar <- function(x) {
    if (length(x) != 1L) {
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
      sprintf("# dist_code = %s", format_scalar(row$dist_code[[1]])),
      sprintf("# factor5 = %s", format_scalar(row$factor5[[1]])),
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
    design_block = design$design_block,
    dist_code = design$dist_code,
    G = design$G,
    ln_sd = design$ln_sd,
    nb_size = design$nb_size,
    Nmin = design$Nmin,
    Nmax = design$Nmax,
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
#
# If you changed fixed spread baselines in the design builder, use:
# write_hake_inp_files_from_matrix(
#   matrix_file = "hake_factorial_design_matrix.csv",
#   output_dir = ".",
#   expected_lognormal_ln_sd = 1.25,
#   expected_nb_size = 50,
#   expected_loguniform_Nmin = 20,
#   expected_loguniform_Nmax = 120
# )
