# code-C2-revised_arrow_safe.R
# Memory-safe factorial MANOVA on pairwise method differences.
#
# Compatible with:
#   1. acc_long as an Arrow Dataset created by:
#        acc_long <- open_metadata_dataset("acc_long")
#   2. acc_long as an ordinary data.frame/data.table.
#
# For Arrow input, this script reads acc_long RDS partitions one at a time,
# converts each partition to pairwise multivariate differences in ordinary R,
# and retains a reproducible bounded sample for each MANOVA.
#
# Optional controls:
#   Sys.setenv(HAKE_CODE_C2_MAX_ROWS = "500000")
#   Sys.setenv(HAKE_CODE_C2_SEED = "20260731")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

truthy_C2 <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "yes", "y", "1")
}

find_acc_long_parts_C2 <- function() {
  out_dir <- Sys.getenv(
    "HAKE_METADATA_OUTDIR",
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
      "Run metadata_revised.R first, or set HAKE_METADATA_OUTDIR to the ",
      "correct metadata_output directory."
    )
  }

  files
}

standardize_metric_C2 <- function(x) {
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

make_pairwise_diff_partition_C2 <- function(dat, m_target, m_ref) {
  id_vars <- c(
    "example_id", "design_block", "G",
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
  dt[, metric := standardize_metric_C2(metric)]
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

  # Preserve the original C2 contrast direction:
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

sample_pairwise_partitions_C2 <- function(files,
                                          m_target,
                                          m_ref,
                                          max_rows,
                                          seed) {
  quota <- ceiling(max_rows / length(files) * 1.10)
  pieces <- vector("list", length(files))

  for (i in seq_along(files)) {
    message(
      "  ", m_target, "_vs_", m_ref,
      ": partition ", i, " of ", length(files)
    )

    part <- readRDS(files[[i]])
    diff_part <- make_pairwise_diff_partition_C2(
      part,
      m_target = m_target,
      m_ref = m_ref
    )

    if (nrow(diff_part) > quota) {
      set.seed(seed + i)
      diff_part <- diff_part[
        sample.int(nrow(diff_part), quota, replace = FALSE)
      ]
    }

    pieces[[i]] <- diff_part
    rm(part, diff_part)
    invisible(gc())
  }

  out <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  rm(pieces)
  invisible(gc())

  if (nrow(out) > max_rows) {
    set.seed(seed)
    out <- out[sample.int(nrow(out), max_rows, replace = FALSE)]
  }

  as.data.frame(out)
}

run_code_C2 <- function(data_name = "acc_long",
                        out_file = "results-code-C2.txt") {

  max_rows <- suppressWarnings(as.integer(
    Sys.getenv("HAKE_CODE_C2_MAX_ROWS", unset = "500000")
  ))
  if (is.na(max_rows) || max_rows < 1000L) {
    max_rows <- 500000L
  }

  sample_seed <- suppressWarnings(as.integer(
    Sys.getenv("HAKE_CODE_C2_SEED", unset = "20260731")
  ))
  if (is.na(sample_seed)) {
    sample_seed <- 20260731L
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
    cat("Code C2: Factorial MANOVA on Pairwise Differences\n")
    cat("Output file:", out_file, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    cat("Maximum rows per pairwise MANOVA:", max_rows, "\n")
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
      files <- find_acc_long_parts_C2()
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
        diff_mv <- sample_pairwise_partitions_C2(
          files,
          m_target = m_target,
          m_ref = m_ref,
          max_rows = max_rows,
          seed = sample_seed + 1000L * k
        )
      } else {
        diff_mv <- as.data.frame(
          make_pairwise_diff_partition_C2(
            acc_long_in,
            m_target = m_target,
            m_ref = m_ref
          )
        )

        if (nrow(diff_mv) > max_rows) {
          set.seed(sample_seed + 1000L * k)
          diff_mv <- diff_mv[
            sample.int(nrow(diff_mv), max_rows, replace = FALSE),
            ,
            drop = FALSE
          ]
        }
      }

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

      formula_diff <- cbind(rmse, L1, Linf) ~
        design_block + G + theta_true +
        theta_CV + sigma + mean_nsamp + p1 + p2

      cat("Rows analyzed:", nrow(diff_mv), "\n")
      cat("MANOVA formula:\n")
      print(formula_diff)
      cat("\n")

      fit_diff_manova <- stats::manova(
        formula_diff,
        data = diff_mv
      )

      cat("MANOVA Table (Pillai):\n")
      manova_summary <- summary(fit_diff_manova, test = "Pillai")
      print(manova_summary)
      cat("\n")

      cat("Univariate Breakdown for each Error Metric:\n")
      aov_summary <- summary.aov(fit_diff_manova)
      print(aov_summary)
      cat("\n")

      results[[comp_name]] <- list(
        comparison = comp_name,
        rows_analyzed = nrow(diff_mv),
        formula = formula_diff,
        model = fit_diff_manova,
        manova = manova_summary,
        univariate = aov_summary
      )

      rm(diff_mv, fit_diff_manova)
      invisible(gc())
    }

    saveRDS(
      results,
      file = "results-code-C2.rds",
      compress = "gzip"
    )

    summary_rows <- rbindlist(
      lapply(names(results), function(nm) {
        tab <- as.data.frame(results[[nm]]$manova$stats)
        tab$term <- rownames(tab)
        tab$comparison <- nm
        tab$rows_analyzed <- results[[nm]]$rows_analyzed
        rownames(tab) <- NULL
        as.data.table(tab)
      }),
      use.names = TRUE,
      fill = TRUE
    )

    setcolorder(
      summary_rows,
      c("comparison", "term", "rows_analyzed",
        setdiff(names(summary_rows),
                c("comparison", "term", "rows_analyzed")))
    )

    fwrite(summary_rows, "results-code-C2-summary.csv")

    cat("Analysis completed successfully.\n")
    cat("Saved detailed results: results-code-C2.rds\n")
    cat("Saved MANOVA summary: results-code-C2-summary.csv\n")

    invisible(results)

  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n", sep = "")
    stop(e)
  })
}

code_C2_result <- run_code_C2()
