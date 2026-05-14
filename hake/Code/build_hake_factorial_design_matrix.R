
# file = build_hake_factorial_design_matrix.R
# Build a 3-level factorial design matrix for hake.R across:
#   1 = Poisson
#   2 = Lognormal
#   3 = Negative binomial
#   4 = Log-uniform
#
# Outputs:
#   - hake_factorial_design_matrix.csv
#   - hake_factorial_design_legend.csv
#
# Main interactive entry point:
#   source("build_hake_factorial_design_matrix.R")
#   design <- build_hake_factorial_design_matrix()
#
# Exact recreation of the current 972-row design:
#   design <- build_hake_factorial_design_matrix_default()
#
# Verification:
#   verify_hake_factorial_design_matrix("new_matrix.csv", "hake_factorial_design_matrix.csv")

fmt_num <- function(x) {
  if (length(x) != 1L) stop("fmt_num expects a scalar.")
  if (is.na(x)) return(NA_character_)
  if (isTRUE(all.equal(x, round(x)))) {
    as.character(as.integer(round(x)))
  } else {
    formatC(x, digits = 15, format = "fg", flag = "#")
  }
}

fmt_triplet_label <- function(x) paste(vapply(x, fmt_num, character(1)), collapse = "/")
fmt_bounds <- function(x) paste(vapply(x, fmt_num, character(1)), collapse = ", ")
fmt_pair_triplet_label <- function(a, b) {
  paste(sprintf("(%s,%s)", vapply(a, fmt_num, character(1)), vapply(b, fmt_num, character(1))), collapse = "/")
}

parse_triplet_numeric <- function(txt, factor_name) {
  parts <- unlist(strsplit(trimws(txt), "[,[:space:]]+"))
  parts <- parts[nzchar(parts)]
  if (length(parts) != 3L) stop(sprintf("Factor '%s' must have exactly 3 values.", factor_name))
  vals <- suppressWarnings(as.numeric(parts))
  if (anyNA(vals)) stop(sprintf("Factor '%s' contains a non-numeric entry.", factor_name))
  vals
}

prompt_yes_no <- function(prompt, default = TRUE) {
  suffix <- if (isTRUE(default)) " [Y/n]: " else " [y/N]: "
  ans <- readline(paste0(prompt, suffix))
  ans <- tolower(trimws(ans))
  if (!nzchar(ans)) return(default)
  ans %in% c("y", "yes")
}

prompt_triplet_numeric <- function(factor_name, default_values = NULL, extra_note = NULL) {
  default_txt <- if (!is.null(default_values)) paste(default_values, collapse = ", ") else NULL
  msg <- paste0("Enter 3 values for ", factor_name)
  if (!is.null(extra_note)) msg <- paste0(msg, " (", extra_note, ")")
  if (!is.null(default_txt)) msg <- paste0(msg, " [default: ", default_txt, "]")
  msg <- paste0(msg, ": ")
  ans <- readline(msg)
  if (!nzchar(trimws(ans))) {
    if (is.null(default_values)) stop(sprintf("No values entered for '%s'.", factor_name))
    return(as.numeric(default_values))
  }
  parse_triplet_numeric(ans, factor_name)
}

prompt_triplet_pairs <- function(default_Nmin = c(40, 25, 10), default_Nmax = c(60, 100, 250)) {
  cat("\nLog-uniform requires 3 paired (Nmin, Nmax) levels.\n")
  nmin <- prompt_triplet_numeric("Nmin", default_Nmin)
  nmax <- prompt_triplet_numeric("Nmax", default_Nmax)
  if (any(nmin >= nmax)) stop("Each Nmin value must be strictly less than its paired Nmax value.")
  list(Nmin = nmin, Nmax = nmax)
}

build_hake_factorial_design_matrix_from_spec <- function(
    common_by_block,
    poisson_G,
    lognormal_ln_sd,
    nb_nb_size,
    loguniform_Nmin,
    loguniform_Nmax,
    output_csv = "hake_factorial_design_matrix.csv",
    legend_csv = "hake_factorial_design_legend.csv",
    verify_against = NULL,
    run_verification = TRUE
) {
  if (length(poisson_G) != 3L) stop("poisson_G must have length 3.")
  if (length(lognormal_ln_sd) != 3L) stop("lognormal_ln_sd must have length 3.")
  if (length(nb_nb_size) != 3L) stop("nb_nb_size must have length 3.")
  if (length(loguniform_Nmin) != 3L || length(loguniform_Nmax) != 3L) stop("loguniform_Nmin and loguniform_Nmax must each have length 3.")
  if (any(loguniform_Nmin >= loguniform_Nmax)) stop("Each log-uniform Nmin must be strictly less than its paired Nmax.")

  needed_blocks <- c("Poisson", "Lognormal", "Negative binomial", "Log-uniform")
  if (!identical(sort(names(common_by_block)), sort(needed_blocks))) {
    stop("common_by_block must have exactly these names: Poisson, Lognormal, Negative binomial, Log-uniform")
  }

  for (blk in needed_blocks) {
    req <- c("mean_nsamp", "sigma", "theta_true", "theta_CV")
    if (!identical(sort(names(common_by_block[[blk]])), sort(req))) {
      stop(sprintf("Block '%s' in common_by_block must contain mean_nsamp, sigma, theta_true, and theta_CV.", blk))
    }
    lens <- vapply(common_by_block[[blk]], length, integer(1))
    if (any(lens != 3L)) stop(sprintf("All common-factor triplets for block '%s' must have length 3.", blk))
  }

  base_fixed <- list(
    K = 3L,
    h = 0.05,
    nsims = 10L,
    random.seed = sample.int(.Machine$integer.max, 1L),
    od_mult = 1.0,
    p_lbound = "0.0, 0.0, 0.0",
    p_ubound = "1.0, 1.0, 1.0"
  )

  build_block <- function(design_block, dist_code, factor5_name, factor5_values, common_vals,
                          factor5_is_pair = FALSE, factor5_pair_max = NULL) {
    codes <- c(-1L, 0L, 1L)
    # Match the existing matrix row order exactly:
    # factor5 cycles fastest, then factor4, factor3, factor2, factor1.
    grid <- expand.grid(
      level5_code = codes,
      level4_code = codes,
      level3_code = codes,
      level2_code = codes,
      level1_code = codes,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    grid <- grid[, c("level1_code", "level2_code", "level3_code", "level4_code", "level5_code")]

    idx1 <- match(grid$level1_code, codes)
    idx2 <- match(grid$level2_code, codes)
    idx3 <- match(grid$level3_code, codes)
    idx4 <- match(grid$level4_code, codes)
    idx5 <- match(grid$level5_code, codes)

    out <- data.frame(
      design_block = design_block,
      dist_code = as.integer(dist_code),
      factor1 = "mean_nsamp",
      level1_code = grid$level1_code,
      factor2 = "sigma",
      level2_code = grid$level2_code,
      factor3 = "theta_true",
      level3_code = grid$level3_code,
      factor4 = "theta_CV",
      level4_code = grid$level4_code,
      factor5 = factor5_name,
      level5_code = grid$level5_code,
      K = as.integer(base_fixed$K),
      G = 2L,
      h = base_fixed$h,
      theta_true = common_vals$theta_true[idx3],
      theta_CV = common_vals$theta_CV[idx4],
      Nmin = 50L,
      Nmax = 100L,
      nsims = as.integer(base_fixed$nsims),
      random.seed = as.integer(base_fixed$random.seed),
      od_mult = base_fixed$od_mult,
      sigma = common_vals$sigma[idx2],
      mean_nsamp = common_vals$mean_nsamp[idx1],
      ln_sd = 1.0,
      nb_size = 25,
      N_range_label = NA_character_,
      p_lbound = base_fixed$p_lbound,
      p_ubound = base_fixed$p_ubound,
      stringsAsFactors = FALSE
    )

    if (!factor5_is_pair) {
      vals5 <- factor5_values[idx5]
      if (identical(factor5_name, "G")) out$G <- as.integer(round(vals5))
      if (identical(factor5_name, "ln_sd")) out$ln_sd <- vals5
      if (identical(factor5_name, "nb_size")) out$nb_size <- vals5
    } else {
      nmin_vals <- factor5_values[idx5]
      nmax_vals <- factor5_pair_max[idx5]
      out$Nmin <- as.integer(round(nmin_vals))
      out$Nmax <- as.integer(round(nmax_vals))
      out$N_range_label <- sprintf("%s-%s", vapply(nmin_vals, fmt_num, character(1)), vapply(nmax_vals, fmt_num, character(1)))
    }

    out
  }

  poisson_df <- build_block("Poisson", 1L, "G", poisson_G, common_by_block$Poisson)
  lognormal_df <- build_block("Lognormal", 2L, "ln_sd", lognormal_ln_sd, common_by_block$Lognormal)
  nb_df <- build_block("Negative binomial", 3L, "nb_size", nb_nb_size, common_by_block$`Negative binomial`)
  loguniform_df <- build_block("Log-uniform", 4L, "N_range", loguniform_Nmin, common_by_block$`Log-uniform`,
                               factor5_is_pair = TRUE, factor5_pair_max = loguniform_Nmax)

  design <- rbind(poisson_df, lognormal_df, nb_df, loguniform_df)
  design$example_id <- seq_len(nrow(design))
  design$inp_file <- sprintf("hake_ex%d.inp", design$example_id)

  desired_order <- c(
    "example_id", "inp_file", "design_block", "dist_code",
    "factor1", "level1_code", "factor2", "level2_code", "factor3", "level3_code",
    "factor4", "level4_code", "factor5", "level5_code",
    "K", "G", "h", "theta_true", "theta_CV", "Nmin", "Nmax", "nsims",
    "random.seed", "od_mult", "sigma", "mean_nsamp", "ln_sd", "nb_size",
    "N_range_label", "p_lbound", "p_ubound"
  )
  design <- design[, desired_order]

  legend <- data.frame(
    Block = c("All blocks", "Poisson", "Lognormal", "Negative binomial", "Log-uniform"),
    Type = c("Fixed baseline", "Factors", "Factors", "Factors", "Factors"),
    Definition = c(
      sprintf(
        "K=%s, h=%s, nsims=%s, random.seed=%s, od_mult=%s, p_lbound=(0,0,0), p_ubound=(1,1,1)",
        fmt_num(base_fixed$K), fmt_num(base_fixed$h), fmt_num(base_fixed$nsims),
        fmt_num(base_fixed$random.seed), fmt_num(base_fixed$od_mult)
      ),
      paste0(
        "mean_nsamp: ", fmt_triplet_label(common_by_block$Poisson$mean_nsamp),
        "; sigma: ", fmt_triplet_label(common_by_block$Poisson$sigma),
        "; theta_true: ", fmt_triplet_label(common_by_block$Poisson$theta_true),
        "; theta_CV: ", fmt_triplet_label(common_by_block$Poisson$theta_CV),
        "; G: ", fmt_triplet_label(poisson_G)
      ),
      paste0(
        "mean_nsamp: ", fmt_triplet_label(common_by_block$Lognormal$mean_nsamp),
        "; sigma: ", fmt_triplet_label(common_by_block$Lognormal$sigma),
        "; theta_true: ", fmt_triplet_label(common_by_block$Lognormal$theta_true),
        "; theta_CV: ", fmt_triplet_label(common_by_block$Lognormal$theta_CV),
        "; ln_sd: ", fmt_triplet_label(lognormal_ln_sd)
      ),
      paste0(
        "mean_nsamp: ", fmt_triplet_label(common_by_block$`Negative binomial`$mean_nsamp),
        "; sigma: ", fmt_triplet_label(common_by_block$`Negative binomial`$sigma),
        "; theta_true: ", fmt_triplet_label(common_by_block$`Negative binomial`$theta_true),
        "; theta_CV: ", fmt_triplet_label(common_by_block$`Negative binomial`$theta_CV),
        "; nb_size: ", fmt_triplet_label(nb_nb_size)
      ),
      paste0(
        "mean_nsamp: ", fmt_triplet_label(common_by_block$`Log-uniform`$mean_nsamp),
        "; sigma: ", fmt_triplet_label(common_by_block$`Log-uniform`$sigma),
        "; theta_true: ", fmt_triplet_label(common_by_block$`Log-uniform`$theta_true),
        "; theta_CV: ", fmt_triplet_label(common_by_block$`Log-uniform`$theta_CV),
        "; (Nmin,Nmax): ", fmt_pair_triplet_label(loguniform_Nmin, loguniform_Nmax)
      )
    ),
    stringsAsFactors = FALSE
  )

  utils::write.csv(design, output_csv, row.names = FALSE, na = "")
  utils::write.csv(legend, legend_csv, row.names = FALSE, na = "")

  message(sprintf("Wrote design matrix: %s", normalizePath(output_csv, winslash = "/", mustWork = FALSE)))
  message(sprintf("Wrote legend: %s", normalizePath(legend_csv, winslash = "/", mustWork = FALSE)))
  message(sprintf("Rows in design matrix: %d", nrow(design)))

  verification <- NULL
  if (isTRUE(run_verification) && !is.null(verify_against) && file.exists(verify_against)) {
    verification <- verify_hake_factorial_design_matrix(output_csv, verify_against)
    if (verification$identical) {
      message("Verification result: exact match to reference matrix.")
    } else {
      message("Verification result: NOT an exact match to reference matrix.")
      message(verification$message)
    }
  }
  attr(design, "verification") <- verification
  invisible(design)
}

build_hake_factorial_design_matrix <- function(
    output_csv = "hake_factorial_design_matrix.csv",
    legend_csv = "hake_factorial_design_legend.csv",
    verify_against = output_csv,
    run_verification = TRUE
) {
  common_defaults <- list(
    mean_nsamp = c(25, 50, 100),
    sigma = c(0.25, 0.5, 1.0),
    theta_true = c(0.5, 1.0, 2.0),
    theta_CV = c(0.0, 0.2, 0.6)
  )

  specific_defaults <- list(
    Poisson = list(G = c(2, 4, 8)),
    Lognormal = list(ln_sd = c(0.5, 1.0, 1.5)),
    `Negative binomial` = list(nb_size = c(5, 25, 100)),
    `Log-uniform` = list(Nmin = c(40, 25, 10), Nmax = c(60, 100, 250))
  )

  cat("\nBuild a 3-level factorial design matrix for hake.R\n")
  cat("The output will contain 4 blocks x 3^5 = 972 rows.\n")
  cat("Press Enter to accept the displayed default values.\n")

  use_equal_common <- prompt_yes_no(
    "Use the same triplets for the common factors across all 4 sampling distributions?",
    default = TRUE
  )

  if (use_equal_common) {
    common_triplets <- list(
      mean_nsamp = prompt_triplet_numeric("mean_nsamp", common_defaults$mean_nsamp),
      sigma = prompt_triplet_numeric("sigma", common_defaults$sigma),
      theta_true = prompt_triplet_numeric("theta_true", common_defaults$theta_true),
      theta_CV = prompt_triplet_numeric("theta_CV", common_defaults$theta_CV)
    )
    common_by_block <- list(
      Poisson = common_triplets,
      Lognormal = common_triplets,
      `Negative binomial` = common_triplets,
      `Log-uniform` = common_triplets
    )
  } else {
    blocks <- c("Poisson", "Lognormal", "Negative binomial", "Log-uniform")
    common_by_block <- setNames(vector("list", length(blocks)), blocks)
    for (blk in blocks) {
      cat("\nCommon-factor triplets for ", blk, ":\n", sep = "")
      common_by_block[[blk]] <- list(
        mean_nsamp = prompt_triplet_numeric("mean_nsamp", common_defaults$mean_nsamp),
        sigma = prompt_triplet_numeric("sigma", common_defaults$sigma),
        theta_true = prompt_triplet_numeric("theta_true", common_defaults$theta_true),
        theta_CV = prompt_triplet_numeric("theta_CV", common_defaults$theta_CV)
      )
    }
  }

  cat("\nNow enter the distribution-specific factors.\n")
  poisson_G <- prompt_triplet_numeric("Poisson G", specific_defaults$Poisson$G)
  lognormal_ln_sd <- prompt_triplet_numeric("Lognormal ln_sd", specific_defaults$Lognormal$ln_sd)
  nb_nb_size <- prompt_triplet_numeric("Negative binomial nb_size", specific_defaults$`Negative binomial`$nb_size)
  loguniform_pair <- prompt_triplet_pairs(
    default_Nmin = specific_defaults$`Log-uniform`$Nmin,
    default_Nmax = specific_defaults$`Log-uniform`$Nmax
  )

  build_hake_factorial_design_matrix_from_spec(
    common_by_block = common_by_block,
    poisson_G = poisson_G,
    lognormal_ln_sd = lognormal_ln_sd,
    nb_nb_size = nb_nb_size,
    loguniform_Nmin = loguniform_pair$Nmin,
    loguniform_Nmax = loguniform_pair$Nmax,
    output_csv = output_csv,
    legend_csv = legend_csv,
    verify_against = verify_against,
    run_verification = run_verification
  )
}

build_hake_factorial_design_matrix_default <- function(
    output_csv = "hake_factorial_design_matrix.csv",
    legend_csv = "hake_factorial_design_legend.csv",
    verify_against = output_csv,
    run_verification = TRUE
) {
  common_triplets <- list(
    mean_nsamp = c(25, 50, 100),
    sigma = c(0.25, 0.5, 1.0),
    theta_true = c(0.5, 1.0, 2.0),
    theta_CV = c(0.0, 0.2, 0.6)
  )

  build_hake_factorial_design_matrix_from_spec(
    common_by_block = list(
      Poisson = common_triplets,
      Lognormal = common_triplets,
      `Negative binomial` = common_triplets,
      `Log-uniform` = common_triplets
    ),
    poisson_G = c(2, 4, 8),
    lognormal_ln_sd = c(0.5, 1.0, 1.5),
    nb_nb_size = c(5, 25, 100),
    loguniform_Nmin = c(40, 25, 10),
    loguniform_Nmax = c(60, 100, 250),
    output_csv = output_csv,
    legend_csv = legend_csv,
    verify_against = verify_against,
    run_verification = run_verification
  )
}

verify_hake_factorial_design_matrix <- function(new_file, reference_file = "hake_factorial_design_matrix.csv") {
  if (!file.exists(new_file)) stop(sprintf("New matrix file not found: '%s'", new_file))
  if (!file.exists(reference_file)) stop(sprintf("Reference matrix file not found: '%s'", reference_file))

  new_df <- utils::read.csv(new_file, stringsAsFactors = FALSE, check.names = FALSE)
  ref_df <- utils::read.csv(reference_file, stringsAsFactors = FALSE, check.names = FALSE)

  if (!identical(names(new_df), names(ref_df))) {
    return(list(identical = FALSE, message = "Column names differ."))
  }
  if (nrow(new_df) != nrow(ref_df)) {
    return(list(identical = FALSE, message = sprintf("Row counts differ: new=%d ref=%d.", nrow(new_df), nrow(ref_df))))
  }
  if (ncol(new_df) != ncol(ref_df)) {
    return(list(identical = FALSE, message = sprintf("Column counts differ: new=%d ref=%d.", ncol(new_df), ncol(ref_df))))
  }

  trim_df <- function(df) {
    out <- df
    for (j in seq_along(out)) out[[j]] <- trimws(as.character(out[[j]]))
    out
  }
  new_chr <- trim_df(new_df)
  ref_chr <- trim_df(ref_df)

  if (identical(new_chr, ref_chr)) {
    return(list(identical = TRUE, message = "The recreated matrix exactly matches the reference file."))
  }

  diff_locs <- which(new_chr != ref_chr, arr.ind = TRUE)
  first_diff <- diff_locs[1, , drop = FALSE]
  i <- first_diff[1, "row"]
  j <- first_diff[1, "col"]
  list(
    identical = FALSE,
    message = sprintf(
      "First difference at row %d, column '%s': new='%s', ref='%s'.",
      i, names(new_chr)[j], new_chr[i, j], ref_chr[i, j]
    )
  )
}
