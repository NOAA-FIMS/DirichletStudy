# code-D-K9.r
# 
# 
# rm(list = ls())
# gc()
# source("metadata.R")
#
# Optional controls:
#   Sys.setenv(CODE_D_SAMPLE_FRACTION = "0.02")
#   Sys.setenv(CODE_D_SEED = "65436")
#
# source("code-D-K9.r")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

truthy_D <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "yes", "y", "1")
}

find_acc_long_parts_D <- function() {
  out_dir <- Sys.getenv(
    "METADATA_OUTDIR",
    unset = if (exists("OUT_DIR", inherits = TRUE)) {
      get("OUT_DIR", inherits = TRUE)
    } else {
      file.path(getwd(), "metadata_output")
    }
  )

  part_dir <- file.path(out_dir, "acc_long_parts")
  files <- sort(list.files(
    part_dir,
    pattern = "^acc_long_part_[0-9]+\\.rds$",
    full.names = TRUE
  ))

  if (!length(files)) {
    stop(
      "No acc_long RDS partitions were found in: ", part_dir, "\n",
      "Run metadata_revised.R first, or set METADATA_OUTDIR to the ",
      "correct metadata_output directory."
    )
  }

  files
}

standardize_metric_D <- function(x) {
  z <- tolower(as.character(x))
  fifelse(
    z %chin% c("rmse", "rmse_norm"), "rmse",
    fifelse(
      z %chin% c("l1", "l1_norm", "l1 norm"), "L1",
      fifelse(
        z %chin% c("linf", "l_inf", "linf_norm", "l_inf_norm", "linf norm"),
        "Linf",
        as.character(x)
      )
    )
  )
}

subsample_simplex_units_D <- function(dat,
                                      sample_fraction,
                                      seed) {
  sample_unit_cols <- c("example_id", "mesh_id", "sim_id")
  missing_cols <- setdiff(sample_unit_cols, names(dat))
  if (length(missing_cols)) {
    stop(
      "Simplex-unit subsampling requires column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dt <- as.data.table(dat)
  sample_units <- unique(dt[, ..sample_unit_cols])
  n_units_available <- nrow(sample_units)
  n_units_sampled <- if (n_units_available > 0L) {
    max(1L, as.integer(ceiling(sample_fraction * n_units_available)))
  } else {
    0L
  }

  if (n_units_sampled > 0L && n_units_sampled < n_units_available) {
    set.seed(seed)
    sample_units <- sample_units[
      sample.int(n_units_available, n_units_sampled, replace = FALSE)
    ]
  }

  list(
    data = dt[sample_units, on = sample_unit_cols, nomatch = 0L][],
    n_units_available = n_units_available,
    n_units_sampled = n_units_sampled
  )
}

make_pairwise_diff_partition_D <- function(dat, m_target, m_ref) {
  id_vars <- c(
    "example_id", "mesh_id", "sim_id", "design_block", "G",
    "theta_true", "theta_CV", "sigma", "mean_nsamp",
    "p1", "p2"
  )

  required_cols <- c(id_vars, "method", "metric", "accuracy")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols)) {
    stop(
      "An acc_long partition is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dt <- as.data.table(dat)[, ..required_cols]
  dt <- dt[as.character(method) %chin% c(m_target, m_ref)]

  if (!nrow(dt)) {
    return(data.table())
  }

  dt[, method := as.character(method)]
  dt[, metric := standardize_metric_D(metric)]
  dt[, accuracy := as.numeric(accuracy)]

  dt <- dt[
    is.finite(accuracy) &
      metric %chin% c("rmse", "L1", "Linf")
  ]

  if (!nrow(dt)) {
    return(data.table())
  }

  # Collapse duplicate method-metric cells before widening.
  group_cols <- c(id_vars, "metric", "method")
  dt <- dt[
    ,
    .(accuracy = mean(accuracy, na.rm = TRUE)),
    by = group_cols
  ]

  # Widen methods inside the in-memory partition.
  method_formula <- as.formula(
    paste(paste(c(id_vars, "metric"), collapse = " + "), "~ method")
  )

  wide_method <- dcast(
    dt,
    formula = method_formula,
    value.var = "accuracy",
    fill = NA_real_
  )

  rm(dt)
  invisible(gc())

  if (!all(c(m_target, m_ref) %in% names(wide_method))) {
    return(data.table())
  }

  wide_method <- wide_method[
    !is.na(get(m_target)) & !is.na(get(m_ref))
  ]

  if (!nrow(wide_method)) {
    return(data.table())
  }

  # Preserve the original D contrast direction:
  # target method minus reference method.
  wide_method[, diff := get(m_target) - get(m_ref)]

  # Widen metrics so each row has rmse, L1, and Linf differences.
  metric_formula <- as.formula(
    paste(paste(id_vars, collapse = " + "), "~ metric")
  )

  out <- dcast(
    wide_method[, c(id_vars, "metric", "diff"), with = FALSE],
    formula = metric_formula,
    value.var = "diff",
    fill = NA_real_
  )

  metric_cols <- c("rmse", "L1", "Linf")
  if (!all(metric_cols %in% names(out))) {
    return(data.table())
  }

  out <- out[complete.cases(out[, ..metric_cols])]
  out
}

sample_pairwise_partitions_D <- function(files,
                                          m_target,
                                          m_ref,
                                          sample_fraction,
                                          seed) {
  pieces <- vector("list", length(files))
  sampling_counts <- data.table(
    partition = basename(files),
    units_available = integer(length(files)),
    units_sampled = integer(length(files))
  )

  for (i in seq_along(files)) {
    message(
      "  ", m_target, "_vs_", m_ref,
      ": partition ", i, " of ", length(files)
    )

    part <- readRDS(files[[i]])
    sampled <- subsample_simplex_units_D(
      part,
      sample_fraction = sample_fraction,
      seed = seed + i
    )
    sampling_counts$units_available[i] <- sampled$n_units_available
    sampling_counts$units_sampled[i] <- sampled$n_units_sampled

    diff_part <- make_pairwise_diff_partition_D(
      sampled$data,
      m_target = m_target,
      m_ref = m_ref
    )

    pieces[[i]] <- diff_part
    rm(part, sampled, diff_part)
    invisible(gc())
  }

  out <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  rm(pieces)
  invisible(gc())

  out <- as.data.frame(out)
  attr(out, "sampling_counts") <- sampling_counts
  out
}

# Select a maximal, deterministic subset of response variables whose fitted
# residuals are linearly independent.  summary.manova() requires the residual
# response matrix to have full column rank.  With composition error metrics,
# exact identities (for example, between RMSE and L1 when G = 2) can otherwise
# make a three-response MANOVA undefined.
independent_responses_D <- function(fit,
                                    response_names,
                                    tol = 1e-7) {
  residual_matrix <- as.matrix(stats::residuals(fit))

  if (ncol(residual_matrix) != length(response_names)) {
    stop("The fitted residual matrix does not match the response names.")
  }

  residual_scale <- sqrt(colSums(residual_matrix^2))
  usable <- is.finite(residual_scale) & residual_scale > 0

  if (!any(usable)) {
    stop("All MANOVA responses have zero or non-finite residual variance.")
  }

  scaled_residuals <- residual_matrix[, usable, drop = FALSE]
  scaled_residuals <- sweep(
    scaled_residuals,
    MARGIN = 2,
    STATS = residual_scale[usable],
    FUN = "/"
  )
  usable_names <- response_names[usable]

  # Add responses in the stated order and retain a response only when it
  # increases the numerical rank.  This makes the choice reproducible.
  keep <- integer(0)
  current_rank <- 0L

  for (j in seq_along(usable_names)) {
    candidate <- scaled_residuals[, c(keep, j), drop = FALSE]
    candidate_rank <- qr(candidate, tol = tol, LAPACK = FALSE)$rank

    if (candidate_rank > current_rank) {
      keep <- c(keep, j)
      current_rank <- candidate_rank
    }
  }

  kept_names <- usable_names[keep]

  list(
    rank = current_rank,
    kept = kept_names,
    dropped = setdiff(response_names, kept_names),
    tolerance = tol
  )
}

make_formula_D <- function(response_names, predictor_names) {
  stats::as.formula(paste0(
    "cbind(", paste(response_names, collapse = ", "), ") ~ ",
    paste(predictor_names, collapse = " + ")
  ))
}

run_code_D <- function(data_name = "acc_long",
                        out_file = "results-code-D-K9.txt") {

  sample_fraction <- suppressWarnings(as.numeric(
    Sys.getenv("CODE_D_SAMPLE_FRACTION", unset = "0.02")
  ))
  if (!is.finite(sample_fraction) ||
      sample_fraction <= 0 || sample_fraction > 1) {
    stop(
      "CODE_D_SAMPLE_FRACTION must be greater than 0 and no greater than 1."
    )
  }

  sample_seed <- suppressWarnings(as.integer(
    Sys.getenv("CODE_D_SEED", unset = "65436")
  ))
  if (is.na(sample_seed)) {
    sample_seed <- 65436L
  }

  con <- file(out_file, open = "wt")
  sink_start <- sink.number(type = "output")
  sink(con, type = "output")

  on.exit({
    while (sink.number(type = "output") > sink_start) {
      sink(type = "output")
    }
    if (isOpen(con)) close(con)
  }, add = TRUE)

  tryCatch({
    cat("Code D: Factorial MANOVA on Pairwise Differences\n")
    cat("Output file:", out_file, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    cat(
      "Simplex sampling fraction:", sample_fraction,
      "(", 100 * sample_fraction, "%)\n",
      sep = ""
    )
    cat("Sampling unit: example_id x mesh_id x sim_id\n")
    cat("Sampling seed:", sample_seed, "\n\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' not found.")
    }

    acc_long_in <- get(data_name, envir = parent.frame(), inherits = TRUE)

    is_arrow_input <- inherits(
      acc_long_in,
      c(
        "Dataset", "FileSystemDataset", "UnionDataset",
        "arrow_dplyr_query", "RecordBatchReader"
      )
    )

    if (is_arrow_input) {
      files <- find_acc_long_parts_D()
      cat("Input type: disk-backed Arrow acc_long dataset\n")
      cat("Processing mode: partition-wise RDS processing\n")
      cat("Number of acc_long partitions:", length(files), "\n\n")
    } else if (is.data.frame(acc_long_in) || data.table::is.data.table(acc_long_in)) {
      files <- NULL
      cat("Input type: in-memory acc_long table\n")
      cat("Processing mode: direct in-memory processing\n\n")
    } else {
      stop(
        "Unsupported acc_long object class: ",
        paste(class(acc_long_in), collapse = ", ")
      )
    }

    methods <- c("i", "ii", "iii", "iv", "v")
    pairs <- combn(methods, 2, simplify = FALSE)

    results <- vector("list", length(pairs))
    names(results) <- vapply(
      pairs,
      function(x) paste0(x[2], "_vs_", x[1]),
      character(1)
    )

    for (k in seq_along(pairs)) {
      p <- pairs[[k]]
      m_target <- p[2]
      m_ref <- p[1]
      comp_name <- paste0(m_target, "_vs_", m_ref)

      cat("==================================================================\n")
      cat("Factorial Analysis of Difference:", comp_name, "\n")
      cat("==================================================================\n")

      if (is_arrow_input) {
        diff_mv <- sample_pairwise_partitions_D(
          files,
          m_target = m_target,
          m_ref = m_ref,
          sample_fraction = sample_fraction,
          seed = sample_seed
        )
      } else {
        sampled <- subsample_simplex_units_D(
          acc_long_in,
          sample_fraction = sample_fraction,
          seed = sample_seed
        )
        sampling_counts <- data.table(
          partition = "in_memory",
          units_available = sampled$n_units_available,
          units_sampled = sampled$n_units_sampled
        )
        diff_mv <- as.data.frame(
          make_pairwise_diff_partition_D(
            sampled$data,
            m_target = m_target,
            m_ref = m_ref
          )
        )
        attr(diff_mv, "sampling_counts") <- sampling_counts
        rm(sampled, sampling_counts)
      }

      sampling_counts <- attr(diff_mv, "sampling_counts")
      total_units_available <- sum(sampling_counts$units_available)
      total_units_sampled <- sum(sampling_counts$units_sampled)

      cat("Simplex units available:", total_units_available, "\n")
      cat("Simplex units sampled:", total_units_sampled, "\n")
      cat(
        "Realized simplex sampling percentage:",
        100 * total_units_sampled / total_units_available,
        "\n"
      )

      required_model_cols <- c(
        "rmse", "L1", "Linf",
        "design_block", "G", "theta_true",
        "theta_CV", "sigma", "mean_nsamp", "p1", "p2"
      )

      missing_model_cols <- setdiff(required_model_cols, names(diff_mv))
      if (length(missing_model_cols)) {
        stop(
          "Pairwise table for ", comp_name,
          " is missing column(s): ",
          paste(missing_model_cols, collapse = ", ")
        )
      }

      diff_mv <- diff_mv[
        complete.cases(diff_mv[, required_model_cols, drop = FALSE]),
        ,
        drop = FALSE
      ]

      if (nrow(diff_mv) < 20L) {
        stop(
          "Too few complete rows for ", comp_name,
          ": ", nrow(diff_mv)
        )
      }

      factor_cols <- c(
        "design_block", "G", "theta_true",
        "theta_CV", "sigma", "mean_nsamp"
      )
      for (fc in factor_cols) {
        diff_mv[[fc]] <- droplevels(as.factor(diff_mv[[fc]]))
      }

      response_names <- c("rmse", "L1", "Linf")
      predictor_names <- c(
        "design_block", "G", "theta_true", "theta_CV",
        "sigma", "mean_nsamp", "p1", "p2"
      )
      formula_diff_full <- make_formula_D(
        response_names,
        predictor_names
      )

      cat("Rows analyzed:", nrow(diff_mv), "\n")
      cat("Full multivariate formula:\n")
      print(formula_diff_full)
      cat("\n")

      # Fit all three responses first.  This model remains valid for the
      # separate univariate ANOVAs even if their joint residual matrix is
      # rank deficient.
      fit_diff_full <- stats::manova(
        formula_diff_full,
        data = diff_mv
      )

      response_check <- independent_responses_D(
        fit_diff_full,
        response_names
      )

      if (response_check$rank < 2L) {
        stop(
          "The residual response matrix for ", comp_name,
          " has rank ", response_check$rank,
          "; at least two independent responses are required for MANOVA."
        )
      }

      cat(
        "Residual response rank:", response_check$rank,
        "of", length(response_names), "\n"
      )
      cat(
        "Responses used for the joint MANOVA:",
        paste(response_check$kept, collapse = ", "), "\n"
      )
      if (length(response_check$dropped)) {
        cat(
          "Responses excluded from the joint MANOVA because of exact ",
          "residual linear dependence: ",
          paste(response_check$dropped, collapse = ", "), "\n",
          sep = ""
        )
      }
      cat("\n")

      if (length(response_check$kept) == length(response_names)) {
        formula_diff_manova <- formula_diff_full
        fit_diff_manova <- fit_diff_full
      } else {
        formula_diff_manova <- make_formula_D(
          response_check$kept,
          predictor_names
        )
        fit_diff_manova <- stats::manova(
          formula_diff_manova,
          data = diff_mv
        )
      }

      cat("Joint MANOVA formula:\n")
      print(formula_diff_manova)
      cat("\n")

      cat("MANOVA Table (Pillai):\n")
      manova_summary <- summary(fit_diff_manova, test = "Pillai")
      print(manova_summary)
      cat("\n")

      cat("Univariate Breakdown for each Error Metric:\n")
      aov_summary <- summary.aov(fit_diff_full)
      print(aov_summary)
      cat("\n")

      results[[comp_name]] <- list(
        comparison = comp_name,
        rows_analyzed = nrow(diff_mv),
        sampling_fraction = sample_fraction,
        sampling_counts = sampling_counts,
        formula = formula_diff_manova,
        full_formula = formula_diff_full,
        manova_formula = formula_diff_manova,
        residual_response_rank = response_check$rank,
        manova_responses = response_check$kept,
        dropped_responses = response_check$dropped,
        full_model = fit_diff_full,
        model = fit_diff_manova,
        manova = manova_summary,
        univariate = aov_summary
      )

      rm(diff_mv, fit_diff_full, fit_diff_manova)
      invisible(gc())
    }

    saveRDS(
      results,
      file = "code-D-K9-model-results.rds",
      compress = "gzip"
    )

    summary_rows <- rbindlist(
      lapply(names(results), function(nm) {
        tab <- as.data.frame(results[[nm]]$manova$stats)
        tab$term <- rownames(tab)
        tab$comparison <- nm
        tab$rows_analyzed <- results[[nm]]$rows_analyzed
        tab$residual_response_rank <-
          results[[nm]]$residual_response_rank
        tab$manova_responses <- paste(
          results[[nm]]$manova_responses,
          collapse = ";"
        )
        tab$dropped_responses <- paste(
          results[[nm]]$dropped_responses,
          collapse = ";"
        )
        rownames(tab) <- NULL
        as.data.table(tab)
      }),
      use.names = TRUE,
      fill = TRUE
    )

    setcolorder(
      summary_rows,
      c("comparison", "term", "rows_analyzed",
        "residual_response_rank", "manova_responses",
        "dropped_responses",
        setdiff(names(summary_rows),
                c("comparison", "term", "rows_analyzed",
                  "residual_response_rank", "manova_responses",
                  "dropped_responses")))
    )

    fwrite(summary_rows, "results-code-D-K9-summary.csv")

    cat("Analysis completed successfully.\n")
    cat("Saved detailed results: code-D-K9-model-results.rds\n")
    cat("Saved MANOVA summary: results-code-D-K9-summary.csv\n")

    invisible(results)

  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n", sep = "")
    stop(e)
  })
}

code_D_result <- run_code_D()
