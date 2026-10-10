# code-C.R
#
# source("code-C.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

find_acc_long_parts_C <- function() {
  out_dir <- Sys.getenv(
    "METADATA_OUTDIR",
    unset = if (exists("OUT_DIR", inherits = TRUE)) {
      get("OUT_DIR", inherits = TRUE)
    } else {
      file.path(getwd(), "metadata_output")
    }
  )
  part_dir <- file.path(out_dir, "acc_long_parts")
  files <- sort(list.files(part_dir,
    pattern = "^acc_long_part_[0-9]+\\.rds$", full.names = TRUE))
  if (!length(files)) {
    stop("No acc_long RDS partitions found in: ", part_dir,
         "\nRun metadata_revised.R first, or set METADATA_OUTDIR correctly.")
  }
  files
}

standardize_metric_C <- function(x) {
  z <- tolower(as.character(x))
  fifelse(z %chin% c("rmse", "rmse_norm"), "rmse",
    fifelse(z %chin% c("l1", "l1_norm", "l1 norm"), "L1",
      fifelse(z %chin% c("linf", "l_inf", "linf_norm", "l_inf_norm", "linf norm"),
        "Linf", as.character(x))))
}

make_pairwise_diff_partition_C <- function(dat, m_target, m_ref, eps = 1e-6) {
  id_vars <- c("example_id", "design_block", "G", "theta_true", "theta_CV",
               "sigma", "mean_nsamp", "p1", "p2")
  required_cols <- c(id_vars, "method", "metric", "accuracy")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols)) {
    stop("An acc_long partition is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  dt <- as.data.table(dat)[, ..required_cols]
  dt <- dt[as.character(method) %chin% c(m_target, m_ref)]
  if (!nrow(dt)) return(data.table())

  dt[, method := as.character(method)]
  dt[, metric := standardize_metric_C(metric)]
  dt[, accuracy := as.numeric(accuracy)]
  dt[, accuracy_log := log(accuracy + eps)]
  dt <- dt[is.finite(accuracy_log) & metric %chin% c("rmse", "L1", "Linf")]
  if (!nrow(dt)) return(data.table())

  group_cols <- c(id_vars, "metric", "method")
  dt <- dt[, .(accuracy_log = mean(accuracy_log, na.rm = TRUE)), by = group_cols]

  method_formula <- as.formula(
    paste(paste(c(id_vars, "metric"), collapse = " + "), "~ method"))
  wide_method <- dcast(dt, method_formula, value.var = "accuracy_log", fill = NA_real_)
  rm(dt); invisible(gc())

  if (!all(c(m_target, m_ref) %in% names(wide_method))) return(data.table())
  wide_method <- wide_method[!is.na(get(m_target)) & !is.na(get(m_ref))]
  if (!nrow(wide_method)) return(data.table())
  wide_method[, diff := get(m_target) - get(m_ref)]

  metric_formula <- as.formula(paste(paste(id_vars, collapse = " + "), "~ metric"))
  out <- dcast(wide_method[, c(id_vars, "metric", "diff"), with = FALSE],
               metric_formula, value.var = "diff", fill = NA_real_)

  metric_cols <- c("rmse", "L1", "Linf")
  if (!all(metric_cols %in% names(out))) return(data.table())
  out[complete.cases(out[, ..metric_cols])]
}

new_accumulator_C <- function() {
  list(n = 0L, sum = numeric(3), crossprod = matrix(0, 3, 3))
}

update_accumulator_C <- function(acc, y) {
  if (!nrow(y)) return(acc)
  y <- as.matrix(y)
  storage.mode(y) <- "double"
  y <- y[complete.cases(y), , drop = FALSE]
  if (!nrow(y)) return(acc)
  acc$n <- acc$n + nrow(y)
  acc$sum <- acc$sum + colSums(y)
  acc$crossprod <- acc$crossprod + crossprod(y)
  acc
}

finalize_hotelling_C <- function(acc, response_cols = c("rmse", "L1", "Linf")) {
  n <- acc$n; p <- length(response_cols)
  if (n <= p) stop("Hotelling's T^2 requires more complete rows than responses. n=", n)
  ybar <- acc$sum / n
  s <- (acc$crossprod - n * tcrossprod(ybar)) / (n - 1)
  solved <- tryCatch(solve(s, ybar), error = function(e) qr.solve(s, ybar))
  t2 <- as.numeric(n * crossprod(ybar, solved))
  num_df <- p; den_df <- n - p
  approx_f <- ((n - p) / (p * (n - 1))) * t2
  p_value <- pf(approx_f, num_df, den_df, lower.tail = FALSE)
  pillai <- t2 / (t2 + n - 1)

  manova_table <- data.frame(Df = 1L, Pillai = pillai, `approx F` = approx_f,
    `num Df` = num_df, `den Df` = den_df, `Pr(>F)` = p_value,
    check.names = FALSE)
  rownames(manova_table) <- "Intercept"

  mean_table <- data.frame(metric = response_cols,
    mean_log_difference = as.numeric(ybar),
    error_ratio = exp(as.numeric(ybar)),
    percent_error_change = 100 * (exp(as.numeric(ybar)) - 1),
    check.names = FALSE)

  list(n = n, p = p, T2 = t2, covariance = s,
       manova_table = manova_table, mean_table = mean_table)
}

process_pair_from_partitions_C <- function(files, m_target, m_ref, eps) {
  acc <- new_accumulator_C()
  for (i in seq_along(files)) {
    message("  ", m_target, "_vs_", m_ref, ": partition ", i, " of ", length(files))
    part <- readRDS(files[[i]])
    diff_part <- make_pairwise_diff_partition_C(part, m_target, m_ref, eps)
    if (nrow(diff_part)) {
      acc <- update_accumulator_C(acc,
        diff_part[, c("rmse", "L1", "Linf"), with = FALSE])
    }
    rm(part, diff_part); invisible(gc())
  }
  finalize_hotelling_C(acc)
}

process_pair_from_memory_C <- function(dat, m_target, m_ref, eps) {
  diff_mv <- make_pairwise_diff_partition_C(dat, m_target, m_ref, eps)
  acc <- update_accumulator_C(new_accumulator_C(),
    diff_mv[, c("rmse", "L1", "Linf"), with = FALSE])
  finalize_hotelling_C(acc)
}

run_code_C <- function(data_name = "acc_long", out_file = "results-code-C.txt", eps = 1e-6) {
  con <- file(out_file, open = "wt")
  sink_start <- sink.number(type = "output")
  sink(con, type = "output")
  on.exit({
    while (sink.number(type = "output") > sink_start) sink(type = "output")
    if (isOpen(con)) close(con)
  }, add = TRUE)

  tryCatch({
    cat("Code C: Pairwise Hotelling's T^2 / MANOVA results\n")
    cat("Output file:", out_file, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    cat("Transformation: log(accuracy + eps), eps =", eps, "\n")
    cat("Contrast direction: target method minus reference method\n")
    cat("Negative log-differences imply lower error for the target method.\n\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' not found.")
    }
    acc_long_in <- get(data_name, envir = parent.frame(), inherits = TRUE)

    is_arrow_input <- inherits(acc_long_in,
      c("Dataset", "FileSystemDataset", "UnionDataset", "arrow_dplyr_query", "RecordBatchReader"))

    if (is_arrow_input) {
      files <- find_acc_long_parts_C()
      cat("Input type: disk-backed Arrow acc_long dataset\n")
      cat("Processing mode: partition-wise RDS processing\n")
      cat("Number of acc_long partitions:", length(files), "\n\n")
    } else if (is.data.frame(acc_long_in) || is.data.table(acc_long_in)) {
      files <- NULL
      cat("Input type: in-memory acc_long table\n")
      cat("Processing mode: direct in-memory processing\n\n")
    } else {
      stop("Unsupported acc_long object class: ", paste(class(acc_long_in), collapse = ", "))
    }

    methods <- c("i", "ii", "iii", "iv", "v")
    pairs <- combn(methods, 2, simplify = FALSE)
    all_results <- vector("list", length(pairs))
    names(all_results) <- vapply(pairs, function(x) paste0(x[2], "_vs_", x[1]), character(1))

    for (k in seq_along(pairs)) {
      m_ref <- pairs[[k]][1]; m_target <- pairs[[k]][2]
      comp_name <- paste0(m_target, "_vs_", m_ref)
      cat("==================================================================\n")
      cat("Comparison:", comp_name, "\n")
      cat("==================================================================\n")
      cat("Equivalent MANOVA formula: cbind(rmse, L1, Linf) ~ 1\n")
      cat("Null hypothesis: mean log-difference vector = (0, 0, 0)\n\n")

      hotelling <- if (is_arrow_input) {
        process_pair_from_partitions_C(files, m_target, m_ref, eps)
      } else {
        process_pair_from_memory_C(acc_long_in, m_target, m_ref, eps)
      }
      all_results[[comp_name]] <- hotelling

      cat("Pairwise one-sample Hotelling's T^2 / Pillai test\n")
      cat("Complete paired rows:", hotelling$n, "\n")
      cat("Response dimension:", hotelling$p, "\n")
      cat("Hotelling T^2:", format(hotelling$T2, scientific = TRUE, digits = 8), "\n\n")
      print(hotelling$manova_table, digits = 8)
      cat("\nMean log-differences and relative error changes\n")
      print(hotelling$mean_table, digits = 8, row.names = FALSE)
      cat("\n")
    }

    saveRDS(all_results, "code-C-model-results.rds", compress = "gzip")

    summary_rows <- rbindlist(lapply(names(all_results), function(nm) {
      x <- all_results[[nm]]
      ans <- as.data.table(x$mean_table)
      ans[, comparison := nm]
      ans[, complete_paired_rows := x$n]
      ans[, Hotelling_T2 := x$T2]
      ans[, Pillai := x$manova_table$Pillai]
      ans[, approx_F := x$manova_table[["approx F"]]]
      ans[, p_value := x$manova_table[["Pr(>F)"]]]
      ans
    }), use.names = TRUE, fill = TRUE)

    setcolorder(summary_rows, c("comparison", "metric", "complete_paired_rows",
      "mean_log_difference", "error_ratio", "percent_error_change",
      "Hotelling_T2", "Pillai", "approx_F", "p_value"))
    fwrite(summary_rows, "results-code-C-summary.csv")

    cat("Analysis completed successfully.\n")
    cat("Saved detailed results: results-code-C.rds\n")
    cat("Saved summary table: results-code-C-summary.csv\n")
    invisible(all_results)
  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n", sep = "")
    stop(e)
  })
}

code_C_result <- run_code_C()
