# file = build_hake_factorial_design_matrix_G123.R
#
# Build a 5-factor, 3-level factorial design matrix for hake.R across:
#   1 = Poisson
#   2 = Lognormal
#   3 = Negative binomial
#
# Optional support for Log-uniform, dist_code = 4, is retained by setting
# include_loguniform = TRUE in the interactive or default builder.
#
# Scientific design revision:
#   The five common factorial factors are:
#     factor1 = mean_nsamp
#     factor2 = sigma
#     factor3 = theta_true
#     factor4 = theta_CV
#     factor5 = G
#
#   The distribution-specific spread parameters are fixed baseline values:
#     Lognormal:          ln_sd = 1.0
#     Negative binomial:  nb_size = 25
#     Log-uniform:        Nmin = 25, Nmax = 100
#
#   These fixed spread parameters can be changed interactively in
#   build_hake_factorial_design_matrix().
#
# Outputs:
#   - hake_factorial_design_matrix.csv
#   - hake_factorial_design_legend.csv
#
# Default row count:
#   - 3 sampling distributions x 3^5 = 729 rows.
#   - 4 sampling distributions x 3^5 = 972 rows when include_loguniform = TRUE.
#
# Main interactive entry point:
#   source("build_hake_factorial_design_matrix.R")
#   design <- build_hake_factorial_design_matrix()
#
# Non-interactive default design:
# source("build_hake_factorial_design_matrix_G123.R")
# design <- build_hake_factorial_design_matrix_default(
#  output_csv = "hake_factorial_design_matrix.csv",
#  legend_csv = "hake_factorial_design_legend.csv"
# )
#
# Verification:
#   verify_hake_factorial_design_matrix("new_matrix.csv", "reference_matrix.csv")

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
fmt_scalar_label <- function(x) fmt_num(x)
fmt_pair_label <- function(a, b) sprintf("(%s,%s)", fmt_num(a), fmt_num(b))

parse_triplet_numeric <- function(txt, factor_name) {
  parts <- unlist(strsplit(trimws(txt), "[,[:space:]]+"))
  parts <- parts[nzchar(parts)]
  if (length(parts) != 3L) stop(sprintf("Factor '%s' must have exactly 3 values.", factor_name))
  vals <- suppressWarnings(as.numeric(parts))
  if (anyNA(vals)) stop(sprintf("Factor '%s' contains a non-numeric entry.", factor_name))
  vals
}

parse_scalar_numeric <- function(txt, factor_name) {
  parts <- unlist(strsplit(trimws(txt), "[,[:space:]]+"))
  parts <- parts[nzchar(parts)]
  if (length(parts) != 1L) stop(sprintf("Parameter '%s' must have exactly 1 value.", factor_name))
  val <- suppressWarnings(as.numeric(parts))
  if (anyNA(val)) stop(sprintf("Parameter '%s' contains a non-numeric entry.", factor_name))
  val
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

prompt_scalar_numeric <- function(parameter_name, default_value = NULL, extra_note = NULL) {
  default_txt <- if (!is.null(default_value)) fmt_num(default_value) else NULL
  msg <- paste0("Enter 1 value for ", parameter_name)
  if (!is.null(extra_note)) msg <- paste0(msg, " (", extra_note, ")")
  if (!is.null(default_txt)) msg <- paste0(msg, " [default: ", default_txt, "]")
  msg <- paste0(msg, ": ")
  ans <- readline(msg)
  if (!nzchar(trimws(ans))) {
    if (is.null(default_value)) stop(sprintf("No value entered for '%s'.", parameter_name))
    return(as.numeric(default_value))
  }
  parse_scalar_numeric(ans, parameter_name)
}

coerce_scalar <- function(x, parameter_name) {
  if (length(x) != 1L) {
    stop(sprintf("Parameter '%s' must have length 1 in the common-G design.", parameter_name))
  }
  x <- suppressWarnings(as.numeric(x))
  if (anyNA(x)) stop(sprintf("Parameter '%s' must be numeric.", parameter_name))
  x
}

coerce_triplet <- function(x, factor_name, integer_required = FALSE) {
  if (length(x) != 3L) stop(sprintf("Factor '%s' must have length 3.", factor_name))
  x <- suppressWarnings(as.numeric(x))
  if (anyNA(x)) stop(sprintf("Factor '%s' must be numeric.", factor_name))
  if (integer_required && any(abs(x - round(x)) > sqrt(.Machine$double.eps))) {
    stop(sprintf("Factor '%s' must contain integer-valued levels.", factor_name))
  }
  x
}

build_hake_factorial_design_matrix_from_spec <- function(
    common_by_block,
    G_levels,
    lognormal_ln_sd = 1.0,
    nb_nb_size = 25,
    loguniform_Nmin = 25,
    loguniform_Nmax = 100,
    output_csv = "hake_factorial_design_matrix.csv",
    legend_csv = "hake_factorial_design_legend.csv",
    verify_against = NULL,
    run_verification = TRUE
) {
  G_levels <- coerce_triplet(G_levels, "G_levels", integer_required = TRUE)
  lognormal_ln_sd <- coerce_scalar(lognormal_ln_sd, "lognormal_ln_sd")
  nb_nb_size <- coerce_scalar(nb_nb_size, "nb_nb_size")
  loguniform_Nmin <- coerce_scalar(loguniform_Nmin, "loguniform_Nmin")
  loguniform_Nmax <- coerce_scalar(loguniform_Nmax, "loguniform_Nmax")

  if (lognormal_ln_sd <= 0) stop("lognormal_ln_sd must be positive.")
  if (nb_nb_size <= 0) stop("nb_nb_size must be positive.")
  if (loguniform_Nmin < 1) stop("loguniform_Nmin must be at least 1.")
  if (abs(loguniform_Nmin - round(loguniform_Nmin)) > sqrt(.Machine$double.eps)) {
    stop("loguniform_Nmin must be integer-valued.")
  }
  if (abs(loguniform_Nmax - round(loguniform_Nmax)) > sqrt(.Machine$double.eps)) {
    stop("loguniform_Nmax must be integer-valued.")
  }
  if (loguniform_Nmin >= loguniform_Nmax) {
    stop("loguniform_Nmin must be strictly less than loguniform_Nmax.")
  }

  allowed_blocks <- c("Poisson", "Lognormal", "Negative binomial", "Log-uniform")
  block_dist_codes <- c(
    Poisson = 1L,
    Lognormal = 2L,
    `Negative binomial` = 3L,
    `Log-uniform` = 4L
  )
  supplied_blocks <- names(common_by_block)

  if (is.null(supplied_blocks) || any(!nzchar(supplied_blocks))) {
    stop("common_by_block must be a named list of sampling-distribution blocks.")
  }
  if (anyDuplicated(supplied_blocks)) {
    stop("common_by_block contains duplicate sampling-distribution block names.")
  }
  if (!all(supplied_blocks %in% allowed_blocks)) {
    stop(
      "common_by_block contains unsupported block names: ",
      paste(setdiff(supplied_blocks, allowed_blocks), collapse = ", "),
      ". Supported names are: ", paste(allowed_blocks, collapse = ", ")
    )
  }

  needed_blocks <- allowed_blocks[allowed_blocks %in% supplied_blocks]
  if (length(needed_blocks) < 1L) {
    stop("common_by_block must contain at least one supported sampling-distribution block.")
  }

  for (blk in needed_blocks) {
    req <- c("mean_nsamp", "sigma", "theta_true", "theta_CV")
    if (!identical(sort(names(common_by_block[[blk]])), sort(req))) {
      stop(sprintf("Block '%s' in common_by_block must contain mean_nsamp, sigma, theta_true, and theta_CV.", blk))
    }
    lens <- vapply(common_by_block[[blk]], length, integer(1))
    if (any(lens != 3L)) {
      stop(sprintf("All common-factor triplets for block '%s' must have length 3.", blk))
    }
    common_by_block[[blk]]$mean_nsamp <- coerce_triplet(common_by_block[[blk]]$mean_nsamp, paste(blk, "mean_nsamp"))
    common_by_block[[blk]]$sigma <- coerce_triplet(common_by_block[[blk]]$sigma, paste(blk, "sigma"))
    common_by_block[[blk]]$theta_true <- coerce_triplet(common_by_block[[blk]]$theta_true, paste(blk, "theta_true"))
    common_by_block[[blk]]$theta_CV <- coerce_triplet(common_by_block[[blk]]$theta_CV, paste(blk, "theta_CV"))
  }

  base_fixed <- list(
    K = 3L,
    h = 0.03,
    nsims = 5L,
    random.seed = sample.int(.Machine$integer.max, 1L),
    od_mult = 1.0,
    p_lbound = "0.36, 0.01, 0.0",
    p_ubound = "0.97, 0.59, 0.21"
  )

  build_block <- function(design_block, dist_code, common_vals) {
    codes <- c(-1L, 0L, 1L)

    # Row order:
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
      factor5 = "G",
      level5_code = grid$level5_code,
      K = as.integer(base_fixed$K),
      G = as.integer(round(G_levels[idx5])),
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

    if (dist_code == 2L) {
      out$ln_sd <- lognormal_ln_sd
    }
    if (dist_code == 3L) {
      out$nb_size <- nb_nb_size
    }
    if (dist_code == 4L) {
      out$Nmin <- as.integer(round(loguniform_Nmin))
      out$Nmax <- as.integer(round(loguniform_Nmax))
      out$N_range_label <- sprintf("%s-%s", fmt_num(loguniform_Nmin), fmt_num(loguniform_Nmax))
    }

    out
  }

  block_dfs <- lapply(needed_blocks, function(blk) {
    build_block(
      design_block = blk,
      dist_code = as.integer(block_dist_codes[[blk]]),
      common_vals = common_by_block[[blk]]
    )
  })

  design <- do.call(rbind, block_dfs)
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

  expected_rows <- length(needed_blocks) * 3L^5L
  if (nrow(design) != expected_rows) {
    stop(sprintf("Unexpected row count: expected %d rows, found %d rows.", expected_rows, nrow(design)))
  }

  rows_by_block <- table(design$design_block)
  if (!all(rows_by_block[needed_blocks] == 3L^5L)) {
    stop("Each design block must contain exactly 3^5 = 243 rows.")
  }

  if (!all(design$factor5 == "G")) stop("factor5 must be G for all rows.")
  if (!identical(sort(unique(design$G)), sort(as.integer(round(G_levels))))) {
    stop("The G column does not contain exactly the requested G levels.")
  }
  for (blk in needed_blocks) {
    g_by_block <- sort(unique(design$G[design$design_block == blk]))
    if (!identical(g_by_block, sort(as.integer(round(G_levels))))) {
      stop(sprintf("Block '%s' does not contain all requested G levels.", blk))
    }
  }

  if (!all(design$ln_sd[design$dist_code == 2L] == lognormal_ln_sd)) {
    stop("ln_sd is not constant at lognormal_ln_sd for dist_code = 2.")
  }
  if (!all(design$nb_size[design$dist_code == 3L] == nb_nb_size)) {
    stop("nb_size is not constant at nb_nb_size for dist_code = 3.")
  }
  if (!all(design$Nmin[design$dist_code == 4L] == as.integer(round(loguniform_Nmin)))) {
    stop("Nmin is not constant at loguniform_Nmin for dist_code = 4.")
  }
  if (!all(design$Nmax[design$dist_code == 4L] == as.integer(round(loguniform_Nmax)))) {
    stop("Nmax is not constant at loguniform_Nmax for dist_code = 4.")
  }

  fixed_spread_summary <- sprintf(
    "Fixed distribution-specific spread parameters: Lognormal ln_sd=%s; Negative binomial nb_size=%s; Log-uniform Nmin=%s, Nmax=%s",
    fmt_num(lognormal_ln_sd), fmt_num(nb_nb_size), fmt_num(loguniform_Nmin), fmt_num(loguniform_Nmax)
  )

  block_definition <- function(blk) {
    fixed_suffix <- ""
    if (identical(blk, "Lognormal")) {
      fixed_suffix <- paste0("; fixed ln_sd: ", fmt_scalar_label(lognormal_ln_sd))
    } else if (identical(blk, "Negative binomial")) {
      fixed_suffix <- paste0("; fixed nb_size: ", fmt_scalar_label(nb_nb_size))
    } else if (identical(blk, "Log-uniform")) {
      fixed_suffix <- paste0("; fixed (Nmin,Nmax): ", fmt_pair_label(loguniform_Nmin, loguniform_Nmax))
    }

    paste0(
      "mean_nsamp: ", fmt_triplet_label(common_by_block[[blk]]$mean_nsamp),
      "; sigma: ", fmt_triplet_label(common_by_block[[blk]]$sigma),
      "; theta_true: ", fmt_triplet_label(common_by_block[[blk]]$theta_true),
      "; theta_CV: ", fmt_triplet_label(common_by_block[[blk]]$theta_CV),
      "; G: ", fmt_triplet_label(G_levels),
      fixed_suffix
    )
  }

  legend <- data.frame(
    Block = c("All selected blocks", needed_blocks),
    Type = c("Fixed baseline and common fifth factor", rep("Factors", length(needed_blocks))),
    Definition = c(
      paste0(
        sprintf(
          "K=%s, h=%s, nsims=%s, random.seed=%s, od_mult=%s, p_lbound=(0,0,0), p_ubound=(1,1,1)",
          fmt_num(base_fixed$K), fmt_num(base_fixed$h), fmt_num(base_fixed$nsims),
          fmt_num(base_fixed$random.seed), fmt_num(base_fixed$od_mult)
        ),
        "; selected blocks: ", paste(needed_blocks, collapse = ", "),
        "; common G factor: ", fmt_triplet_label(G_levels),
        "; ", fixed_spread_summary
      ),
      vapply(needed_blocks, block_definition, character(1))
    ),
    stringsAsFactors = FALSE
  )

  utils::write.csv(design, output_csv, row.names = FALSE, na = "")
  utils::write.csv(legend, legend_csv, row.names = FALSE, na = "")

  message(sprintf("Wrote design matrix: %s", normalizePath(output_csv, winslash = "/", mustWork = FALSE)))
  message(sprintf("Wrote legend: %s", normalizePath(legend_csv, winslash = "/", mustWork = FALSE)))
  message(sprintf("Rows in design matrix: %d", nrow(design)))
  message("Rows per design block:")
  print(rows_by_block)
  message(fixed_spread_summary)

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
    run_verification = FALSE,
    include_loguniform = FALSE
) {
  common_defaults <- list(
    mean_nsamp = c(100, 200, 400),
    sigma = c(0.25, 0.5, 1.0),
    theta_true = c(0.5, 1.0, 2.0),
    theta_CV = c(0.2, 0.4, 0.8),
    G = c(1, 2, 3)
  )

  fixed_spread_defaults <- list(
    lognormal_ln_sd = 1.0,
    nb_nb_size = 25,
    loguniform_Nmin = 25,
    loguniform_Nmax = 100
  )

  selected_blocks <- c("Poisson", "Lognormal", "Negative binomial")
  if (isTRUE(include_loguniform)) {
    selected_blocks <- c(selected_blocks, "Log-uniform")
  }

  cat("\nBuild a 5-factor, 3-level factorial design matrix for hake.R\n")
  cat(sprintf(
    "The output will contain %d sampling distributions x 3^5 = %d rows.\n",
    length(selected_blocks), length(selected_blocks) * 3L^5L
  ))
  cat("The common fifth factor is G for all selected sampling distributions.\n")
  cat("Distribution-specific spread parameters are fixed baselines but can be changed below.\n")
  cat("Press Enter to accept the displayed default values.\n")

  use_equal_common <- prompt_yes_no(
    "Use the same triplets for mean_nsamp, sigma, theta_true, theta_CV, and G across all selected sampling distributions?",
    default = TRUE
  )

  if (use_equal_common) {
    common_triplets <- list(
      mean_nsamp = prompt_triplet_numeric("mean_nsamp", common_defaults$mean_nsamp),
      sigma = prompt_triplet_numeric("sigma", common_defaults$sigma),
      theta_true = prompt_triplet_numeric("theta_true", common_defaults$theta_true),
      theta_CV = prompt_triplet_numeric("theta_CV", common_defaults$theta_CV)
    )
    G_levels <- prompt_triplet_numeric("G", common_defaults$G, extra_note = "common fifth factor; integer-valued levels")
    common_by_block <- stats::setNames(
      rep(list(common_triplets), length(selected_blocks)),
      selected_blocks
    )
  } else {
    blocks <- selected_blocks
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
    G_levels <- prompt_triplet_numeric("G", common_defaults$G, extra_note = "common fifth factor; integer-valued levels")
  }

  cat("\nNow enter fixed distribution-specific spread parameters.\n")
  lognormal_ln_sd <- prompt_scalar_numeric(
    "Lognormal ln_sd",
    fixed_spread_defaults$lognormal_ln_sd,
    extra_note = "fixed for dist_code = 2"
  )
  nb_nb_size <- prompt_scalar_numeric(
    "Negative binomial nb_size",
    fixed_spread_defaults$nb_nb_size,
    extra_note = "fixed for dist_code = 3"
  )
  if ("Log-uniform" %in% selected_blocks) {
    loguniform_Nmin <- prompt_scalar_numeric(
      "Log-uniform Nmin",
      fixed_spread_defaults$loguniform_Nmin,
      extra_note = "fixed for dist_code = 4"
    )
    loguniform_Nmax <- prompt_scalar_numeric(
      "Log-uniform Nmax",
      fixed_spread_defaults$loguniform_Nmax,
      extra_note = "fixed for dist_code = 4"
    )
  } else {
    loguniform_Nmin <- fixed_spread_defaults$loguniform_Nmin
    loguniform_Nmax <- fixed_spread_defaults$loguniform_Nmax
  }

  build_hake_factorial_design_matrix_from_spec(
    common_by_block = common_by_block,
    G_levels = G_levels,
    lognormal_ln_sd = lognormal_ln_sd,
    nb_nb_size = nb_nb_size,
    loguniform_Nmin = loguniform_Nmin,
    loguniform_Nmax = loguniform_Nmax,
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
    run_verification = FALSE,
    include_loguniform = FALSE
) {
  common_triplets <- list(
    mean_nsamp = c(25, 50, 100),
    sigma = c(0.25, 0.5, 1.0),
    theta_true = c(0.5, 1.0, 2.0),
    theta_CV = c(0.0, 0.2, 0.6)
  )

  selected_blocks <- c("Poisson", "Lognormal", "Negative binomial")
  if (isTRUE(include_loguniform)) {
    selected_blocks <- c(selected_blocks, "Log-uniform")
  }

  common_by_block <- stats::setNames(
    rep(list(common_triplets), length(selected_blocks)),
    selected_blocks
  )

  build_hake_factorial_design_matrix_from_spec(
    common_by_block = common_by_block,
    G_levels = c(1, 2, 3),
    lognormal_ln_sd = 1.0,
    nb_nb_size = 25,
    loguniform_Nmin = 25,
    loguniform_Nmax = 100,
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
    return(list(identical = TRUE, message = "The matrices exactly match."))
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
