# file = write_hake_inp_files_from_matrix_4factorial.R
# Purpose:
#   Read a 4-factor hake.R factorial design matrix and write hake_ex1.inp
#   through hake_exN.inp files for use with the multinomial operating-model
#   version of hake.R and run_hake_examples.R.
#
# Compatibility:
#   This version is compatible with build_hake_multinomial_om_factorial_design_matrix_4factor.R,
#   where the four retained factors are:
#     factor1 = mean_nsamp
#     factor2 = sigma
#     factor3 = theta_true
#     factor4 = G
#
#   The original fourth factor theta_CV is removed from the design matrix and
#   is not written to the .inp files. Because theta_CV is removed, this writer
#   is intended for multinomial_om_flag = 1.
#
#   Optional support for Log-uniform, dist_code = 4, is retained when the
#   matrix includes a Log-uniform block and expected_blocks includes "Log-uniform".
#
# Usage:
# source("write_hake_inp_files_from_matrix_4factorial.R")
# manifest <- write_hake_inp_files_from_matrix(
#  matrix_file = "hake_factorial_design_matrix.csv",
#  output_dir = ".",
#  validate_4factor_design = TRUE,
#  random_seed_generator_seed = 11131
# )
#
# Notes:
#   - The required design-matrix columns are checked exactly; no extra columns
#     such as theta_CV, factor5, level5_code, or N_range_label are allowed.
#   - Each output file is written in "key = value" format expected by hake.R.
#   - Character fields p_lbound and p_ubound are preserved as comma-separated
#     vectors, e.g. "0.36, 0.01, 0.0".
#   - By default, this script replaces the design-matrix random.seed column
#     with a different positive integer random.seed for each experiment before
#     writing the .inp files. Set assign_unique_random_seeds = FALSE to retain
#     the seed values already stored in the design matrix.

write_hake_inp_files_from_matrix <- function(matrix_file = "hake_factorial_design_matrix.csv",
                                             output_dir = ".",
                                             overwrite = TRUE,
                                             verbose = TRUE,
                                             validate_4factor_design = TRUE,
                                             expected_G_levels = c(1, 2, 3),
                                             expected_blocks = c("Poisson", "Lognormal", "Negative binomial"),
                                             expected_lognormal_ln_sd = 1.0,
                                             expected_nb_size = 25,
                                             expected_loguniform_Nmin = 25,
                                             expected_loguniform_Nmax = 100,
                                             multinomial_om_flag = 1L,
                                             assign_unique_random_seeds = TRUE,
                                             random_seed_generator_seed = NULL,
                                             random_seed_max = .Machine$integer.max) {
  if (!file.exists(matrix_file)) {
    stop(sprintf("Input design matrix file not found: '%s'", matrix_file))
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  multinomial_om_flag <- as.integer(multinomial_om_flag)
  if (length(multinomial_om_flag) != 1L || is.na(multinomial_om_flag) || multinomial_om_flag != 1L) {
    stop("multinomial_om_flag must be 1 because theta_CV has been removed from the design and .inp files.")
  }

  design <- utils::read.csv(matrix_file, stringsAsFactors = FALSE, check.names = FALSE)

  required_cols <- c(
    "example_id", "inp_file", "design_block", "dist_code",
    "factor1", "level1_code", "factor2", "level2_code",
    "factor3", "level3_code", "factor4", "level4_code",
    "K", "G", "h", "theta_true",
    "Nmin", "Nmax", "nsims", "random.seed", "od_mult", "sigma",
    "mean_nsamp", "ln_sd", "nb_size", "p_lbound", "p_ubound"
  )

  missing_cols <- setdiff(required_cols, names(design))
  if (length(missing_cols) > 0) {
    stop(sprintf("Design matrix is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  unexpected_cols <- setdiff(names(design), required_cols)
  if (length(unexpected_cols) > 0) {
    stop(sprintf(
      "Design matrix contains unexpected columns for the 4-factor design: %s",
      paste(unexpected_cols, collapse = ", ")
    ))
  }

  forbidden_cols <- c("theta_CV", "factor5", "level5_code", "N_range_label")
  if (any(forbidden_cols %in% names(design))) {
    stop("The 4-factor design matrix must not contain: ",
         paste(intersect(forbidden_cols, names(design)), collapse = ", "))
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
      "K", "G", "h", "theta_true", "Nmin", "Nmax",
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
    if (any(design$random.seed < 1)) stop("random.seed must be a positive integer.")
    if (any(design$mean_nsamp <= 0)) stop("mean_nsamp must be positive.")
    if (any(design$sigma <= 0)) stop("sigma must be positive.")
    if (any(design$theta_true <= 0)) stop("theta_true must be positive.")
    if (any(design$ln_sd <= 0)) stop("ln_sd must be positive.")
    if (any(design$nb_size <= 0)) stop("nb_size must be positive.")

    if (!isTRUE(validate_4factor_design)) {
      return(invisible(design))
    }

    allowed_blocks <- c("Poisson", "Lognormal", "Negative binomial", "Log-uniform")
    expected_blocks <- as.character(expected_blocks)
    if (length(expected_blocks) < 1L) {
      stop("expected_blocks must contain at least one sampling-distribution block.")
    }
    if (any(!expected_blocks %in% allowed_blocks)) {
      stop(
        "expected_blocks contains unsupported block names: ",
        paste(setdiff(expected_blocks, allowed_blocks), collapse = ", "),
        ". Supported names are: ", paste(allowed_blocks, collapse = ", ")
      )
    }

    observed_blocks <- sort(unique(design$design_block))
    if (!identical(observed_blocks, sort(expected_blocks))) {
      stop(
        "4-factor validation expected design_block values: ",
        paste(expected_blocks, collapse = ", "),
        "; observed: ", paste(observed_blocks, collapse = ", ")
      )
    }

    n_factors <- 4L
    expected_n <- length(expected_blocks) * 3L^n_factors
    if (nrow(design) != expected_n) {
      stop(sprintf(
        "4-factor design should have %d rows = %d distributions x 3^4, but found %d rows.",
        expected_n, length(expected_blocks), nrow(design)
      ))
    }

    rows_by_block <- table(design$design_block)
    if (!all(rows_by_block[expected_blocks] == 3L^n_factors)) {
      stop("Each design_block must contain exactly 3^4 = 81 rows.")
    }

    expected_factor_names <- c("mean_nsamp", "sigma", "theta_true", "G")
    got_factor_names <- vapply(seq_len(n_factors), function(j) {
      vals <- unique(trimws(as.character(design[[paste0("factor", j)]])))
      if (length(vals) != 1L) {
        stop(sprintf("factor%d must contain exactly one factor name.", j))
      }
      vals
    }, character(1))

    if (!identical(got_factor_names, expected_factor_names)) {
      stop(sprintf(
        "Expected factor names %s but observed %s.",
        paste(expected_factor_names, collapse = ", "),
        paste(got_factor_names, collapse = ", ")
      ))
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

    if (any(design$dist_code == 4L) &&
        !all(design$Nmin[design$dist_code == 4L] == expected_loguniform_Nmin)) {
      stop(sprintf(
        "For dist_code = 4, Nmin must be constant at %s. If you changed this baseline in the design builder, pass expected_loguniform_Nmin accordingly.",
        expected_loguniform_Nmin
      ))
    }

    if (any(design$dist_code == 4L) &&
        !all(design$Nmax[design$dist_code == 4L] == expected_loguniform_Nmax)) {
      stop(sprintf(
        "For dist_code = 4, Nmax must be constant at %s. If you changed this baseline in the design builder, pass expected_loguniform_Nmax accordingly.",
        expected_loguniform_Nmax
      ))
    }

    invisible(design)
  }

  design <- validate_design_matrix(design)

  generate_unique_positive_random_seeds <- function(n, generator_seed = NULL, max_seed = .Machine$integer.max) {
    if (length(n) != 1L || is.na(n) || n < 1L) {
      stop("n must be a positive scalar when generating random seeds.")
    }

    max_seed_num <- suppressWarnings(as.numeric(max_seed))
    if (length(max_seed_num) != 1L || is.na(max_seed_num) || !is.finite(max_seed_num)) {
      stop("random_seed_max must be a finite positive scalar.")
    }

    max_seed_num <- floor(max_seed_num)
    if (max_seed_num > .Machine$integer.max) {
      max_seed_num <- .Machine$integer.max
    }
    if (max_seed_num < n) {
      stop(sprintf(
        "random_seed_max must be at least the number of experiments. random_seed_max=%s; n=%s.",
        as.character(max_seed_num), as.character(n)
      ))
    }

    if (!is.null(generator_seed)) {
      generator_seed_num <- suppressWarnings(as.numeric(generator_seed))
      if (length(generator_seed_num) != 1L || is.na(generator_seed_num) ||
          !is.finite(generator_seed_num) || generator_seed_num < 1) {
        stop("random_seed_generator_seed must be a positive integer scalar when provided.")
      }
      if (abs(generator_seed_num - round(generator_seed_num)) > sqrt(.Machine$double.eps)) {
        stop("random_seed_generator_seed must be integer-valued when provided.")
      }
      if (generator_seed_num > .Machine$integer.max) {
        stop(sprintf("random_seed_generator_seed must be <= %s.", .Machine$integer.max))
      }

      had_old_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      old_seed <- if (had_old_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
      on.exit({
        if (had_old_seed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }, add = TRUE)

      set.seed(as.integer(round(generator_seed_num)))
    }

    seeds <- sample.int(n = as.integer(max_seed_num), size = as.integer(n), replace = FALSE)
    seeds <- as.integer(seeds)

    if (length(seeds) != n || anyNA(seeds) || any(seeds < 1L) || anyDuplicated(seeds)) {
      stop("Failed to generate a unique positive integer random.seed for each experiment.")
    }

    seeds
  }

  if (isTRUE(assign_unique_random_seeds)) {
    design$random.seed <- generate_unique_positive_random_seeds(
      n = nrow(design),
      generator_seed = random_seed_generator_seed,
      max_seed = random_seed_max
    )

    if (isTRUE(verbose)) {
      message(sprintf(
        "Assigned %d unique positive integer random.seed values before writing .inp files.",
        nrow(design)
      ))
    }
  } else {
    if (anyDuplicated(as.integer(round(design$random.seed)))) {
      warning(
        "assign_unique_random_seeds = FALSE and the design matrix contains duplicate random.seed values. ",
        "The generated .inp files will retain those duplicate seeds."
      )
    }
  }

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
    keys_from_design <- c(
      "K", "G", "h", "theta_true",
      "Nmin", "Nmax", "nsims", "random.seed", "od_mult", "sigma",
      "dist_code", "mean_nsamp", "ln_sd", "nb_size",
      "p_lbound", "p_ubound"
    )

    design_values <- vapply(keys_from_design, function(k) format_scalar(row[[k]]), character(1))
    lines <- paste(keys_from_design, "=", design_values)

    c(
      lines[seq_len(2L)],
      paste("multinomial_om_flag", "=", format_scalar(multinomial_om_flag)),
      lines[-seq_len(2L)]
    )
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
      sprintf("# factor4 = %s", format_scalar(row$factor4[[1]])),
      sprintf("# random.seed = %s", format_scalar(row$random.seed[[1]])),
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
    multinomial_om_flag = multinomial_om_flag,
    random.seed = design$random.seed,
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
# For reproducible unique per-experiment random seeds, use:
# write_hake_inp_files_from_matrix(
#   matrix_file = "hake_factorial_design_matrix.csv",
#   output_dir = ".",
#   random_seed_generator_seed = 8675309
# )
#
# If you need to retain the random.seed values already stored in the design matrix, use:
# write_hake_inp_files_from_matrix(
#   matrix_file = "hake_factorial_design_matrix.csv",
#   output_dir = ".",
#   assign_unique_random_seeds = FALSE
# )
