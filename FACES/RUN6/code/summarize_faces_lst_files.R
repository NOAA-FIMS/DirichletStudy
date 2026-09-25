# file = summarize_faces_lst_files.R
# Summarize faces_ex#.lst files written by run_experiments.R
# NOTE: Change run# to match RUN# folder name
# source("summarize_faces_lst_files.R")

## ---- user settings ----------------------------------------------------------
root_dir <- "."
n_experiments <- 729L

distribution_labels <- c(
  "1" = "Poisson",
  "2" = "Lognormal",
  "3" = "Negative binomial"
)

method_levels <- c("(i)", "(ii)", "(iii)", "(iv)", "(v)")

numeric_pattern <- paste0(
  "[-+]?(?:(?:[0-9]+(?:\\.[0-9]*)?)|(?:\\.[0-9]+))",
  "(?:[eE][-+]?[0-9]+)?"
)

extract_first <- function(x, pattern, group = 1L) {
  m <- regexec(pattern, x, perl = TRUE, ignore.case = TRUE)
  r <- regmatches(x, m)[[1]]
  if (length(r) >= group + 1L) r[group + 1L] else NA_character_
}

empty_emmeans_table <- function() {
  data.frame(method = method_levels, emmean = NA_real_, stringsAsFactors = FALSE)
}

parse_emmeans_metric <- function(lines, metric_index) {
  # faces.R prints three emmeans objects in this order: RMSE, L1, and Linf.
  markers <- grep("^\\s*\\$emmeans\\s*$", lines, perl = TRUE)
  if (length(markers) < metric_index) return(empty_emmeans_table())

  start <- markers[metric_index] + 1L
  if (start > length(lines)) return(empty_emmeans_table())

  later_markers <- grep(
    "^\\s*\\$(?:contrasts|emmeans)\\s*$",
    lines[start:length(lines)],
    perl = TRUE
  )
  end <- if (length(later_markers) == 0L) {
    min(length(lines), start + 30L)
  } else {
    start + later_markers[1L] - 2L
  }
  if (end < start) return(empty_emmeans_table())

  row_pattern <- paste0(
    "^\\s*\\(?([ivx]+)\\)?\\s+(", numeric_pattern, "|NA)\\b"
  )
  values <- setNames(rep(NA_real_, length(method_levels)), method_levels)

  for (line in lines[start:end]) {
    match <- regexec(row_pattern, line, perl = TRUE, ignore.case = TRUE)
    fields <- regmatches(line, match)[[1]]
    if (length(fields) < 3L) next

    method <- switch(
      tolower(fields[2L]),
      "i" = "(i)",
      "ii" = "(ii)",
      "iii" = "(iii)",
      "iv" = "(iv)",
      "v" = "(v)",
      NA_character_
    )
    if (!is.na(method)) values[method] <- suppressWarnings(as.numeric(fields[3L]))
  }

  data.frame(
    method = method_levels,
    emmean = unname(values[method_levels]),
    stringsAsFactors = FALSE
  )
}

parse_p_value <- function(line) {
  match <- regexec(
    "p(?:-value)?\\s*([<=>]+)\\s*([^,[:space:]]+)",
    line,
    perl = TRUE,
    ignore.case = TRUE
  )
  fields <- regmatches(line, match)[[1]]
  if (length(fields) < 3L) return(NA_character_)
  paste0(fields[2L], fields[3L])
}

parse_test_stats <- function(lines, metric_index) {
  out <- list(
    friedman_stat = NA_real_,
    friedman_df = NA_real_,
    friedman_p = NA_character_,
    kendall_w = NA_real_,
    anova_f = NA_real_,
    anova_p = NA_character_
  )

  friedman_lines <- grep("Friedman chi-squared", lines, value = TRUE, ignore.case = TRUE)
  if (length(friedman_lines) >= metric_index) {
    line <- friedman_lines[metric_index]
    out$friedman_stat <- suppressWarnings(as.numeric(extract_first(
      line, paste0("chi-squared\\s*=\\s*(", numeric_pattern, ")")
    )))
    out$friedman_df <- suppressWarnings(as.numeric(extract_first(
      line, paste0("df\\s*=\\s*(", numeric_pattern, ")")
    )))
    out$friedman_p <- parse_p_value(line)
  }

  kendall_lines <- grep("Kendall'?s W", lines, value = TRUE, ignore.case = TRUE)
  if (length(kendall_lines) >= metric_index) {
    out$kendall_w <- suppressWarnings(as.numeric(extract_first(
      kendall_lines[metric_index],
      paste0("W.*?=\\s*(", numeric_pattern, ")")
    )))
  }

  anova_markers <- grep(
    "Type III Analysis of Variance Table",
    lines,
    ignore.case = TRUE
  )
  if (length(anova_markers) >= metric_index) {
    start <- anova_markers[metric_index] + 1L
    next_marker <- anova_markers[anova_markers > anova_markers[metric_index]]
    end <- if (length(next_marker) == 0L) {
      min(length(lines), start + 20L)
    } else {
      next_marker[1L] - 1L
    }

    if (start <= end) {
      method_rows <- grep("^\\s*method\\s+", lines[start:end], value = TRUE, perl = TRUE)
      if (length(method_rows) > 0L) {
        fields <- strsplit(trimws(method_rows[1L]), "\\s+")[[1L]]
        # lmerTest::anova supplies: term, Sum Sq, Mean Sq, NumDF, DenDF,
        # F value, and Pr(>F).
        if (length(fields) >= 7L) {
          out$anova_f <- suppressWarnings(as.numeric(fields[6L]))
          out$anova_p <- fields[7L]
        }
      }
    }
  }

  out
}

detect_distribution <- function(lines) {
  # Preferred source: the diagnostics line written by faces.R, for example:
  # Sample Size Distribution = 2 (Lognormal)
  diagnostic_lines <- grep(
    "Sample Size Distribution\\s*=",
    lines,
    value = TRUE,
    ignore.case = TRUE
  )

  distribution_code <- NA_integer_
  if (length(diagnostic_lines) > 0L) {
    code_text <- extract_first(
      diagnostic_lines[1L],
      "Sample Size Distribution\\s*=\\s*([123])\\b"
    )
    distribution_code <- suppressWarnings(as.integer(code_text))
  }

  # Secondary source: the printed params list from faces.R.
  if (is.na(distribution_code)) {
    marker <- grep("^\\s*\\$dist_code\\s*$", lines, perl = TRUE)
    if (length(marker) > 0L && marker[1L] < length(lines)) {
      candidate <- lines[(marker[1L] + 1L):min(length(lines), marker[1L] + 3L)]
      candidate <- paste(candidate, collapse = " ")
      code_text <- extract_first(candidate, "\\[1\\]\\s*([123])\\b")
      distribution_code <- suppressWarnings(as.integer(code_text))
    }
  }

  distribution <- if (is.na(distribution_code)) {
    NA_character_
  } else {
    unname(distribution_labels[as.character(distribution_code)])
  }

  list(code = distribution_code, label = distribution)
}

extract_experiment_number <- function(path) {
  file_number <- suppressWarnings(as.integer(sub(
    "^faces_ex([0-9]+)\\.lst$",
    "\\1",
    basename(path),
    ignore.case = TRUE
  )))
  folder_number <- suppressWarnings(as.integer(sub(
    "^ex([0-9]+)$",
    "\\1",
    basename(dirname(path)),
    ignore.case = TRUE
  )))

  if (is.na(file_number) || is.na(folder_number) || file_number != folder_number) {
    stop("Folder and file names do not form the required ex#/faces_ex#.lst pair.")
  }
  file_number
}

best_method <- function(table) {
  valid <- which(is.finite(table$emmean))
  if (length(valid) == 0L) return(NA_character_)
  table$method[valid[which.min(table$emmean[valid])]]
}

get_metric_values <- function(table, prefix) {
  values <- setNames(table$emmean, table$method)
  out <- unname(values[method_levels])
  names(out) <- paste0(prefix, "_", c("i", "ii", "iii", "iv", "v"))
  out
}

parse_one_lst <- function(path) {
  lines <- tryCatch(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    error = function(e) readLines(path, warn = FALSE)
  )

  experiment <- extract_experiment_number(path)
  distribution <- detect_distribution(lines)

  rmse_table <- parse_emmeans_metric(lines, 1L)
  l1_table <- parse_emmeans_metric(lines, 2L)
  linf_table <- parse_emmeans_metric(lines, 3L)

  rmse_test <- parse_test_stats(lines, 1L)
  l1_test <- parse_test_stats(lines, 2L)
  linf_test <- parse_test_stats(lines, 3L)

  data.frame(
    experiment = experiment,
    folder = basename(dirname(path)),
    file = basename(path),
    distribution_code = distribution$code,
    distribution = distribution$label,
    rmse_best = best_method(rmse_table),
    l1_best = best_method(l1_table),
    linf_best = best_method(linf_table),

    t(get_metric_values(rmse_table, "rmse")),
    t(get_metric_values(l1_table, "l1")),
    t(get_metric_values(linf_table, "linf")),

    rmse_friedman_chisq = rmse_test$friedman_stat,
    rmse_friedman_df = rmse_test$friedman_df,
    rmse_friedman_p = rmse_test$friedman_p,
    rmse_kendall_w = rmse_test$kendall_w,
    rmse_anova_f = rmse_test$anova_f,
    rmse_anova_p = rmse_test$anova_p,

    l1_friedman_chisq = l1_test$friedman_stat,
    l1_friedman_df = l1_test$friedman_df,
    l1_friedman_p = l1_test$friedman_p,
    l1_kendall_w = l1_test$kendall_w,
    l1_anova_f = l1_test$anova_f,
    l1_anova_p = l1_test$anova_p,

    linf_friedman_chisq = linf_test$friedman_stat,
    linf_friedman_df = linf_test$friedman_df,
    linf_friedman_p = linf_test$friedman_p,
    linf_kendall_w = linf_test$kendall_w,
    linf_anova_f = linf_test$anova_f,
    linf_anova_p = linf_test$anova_p,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

make_frequency_table <- function(data, metric_column) {
  counts <- table(
    distribution = factor(
      data$distribution,
      levels = unname(distribution_labels)
    ),
    method = factor(data[[metric_column]], levels = method_levels),
    useNA = "no"
  )

  output <- as.data.frame.matrix(counts, stringsAsFactors = FALSE)
  output <- data.frame(
    distribution = rownames(output),
    output,
    row.names = NULL,
    check.names = FALSE
  )
  names(output)[-1L] <- method_levels
  output
}

calculate_distance_correlations <- function(data) {
  method_suffixes <- c("i", "ii", "iii", "iv", "v")
  method_names <- setNames(method_levels, method_suffixes)
  metric_names <- c(rmse = "RMSE", l1 = "L1", linf = "Linf")
  metric_pairs <- combn(names(metric_names), 2L, simplify = FALSE)

  results <- list()
  result_index <- 0L

  for (method_suffix in method_suffixes) {
    for (metric_pair in metric_pairs) {
      x_column <- paste0(metric_pair[1L], "_", method_suffix)
      y_column <- paste0(metric_pair[2L], "_", method_suffix)

      x <- data[[x_column]]
      y <- data[[y_column]]
      keep <- is.finite(x) & is.finite(y)
      x <- x[keep]
      y <- y[keep]
      n_complete <- length(x)

      spearman_rho <- NA_real_
      spearman_s <- NA_real_
      p_value <- NA_real_

      # At least three complete experiments and variation in both metrics are
      # required for a meaningful correlation test.
      if (n_complete >= 3L && length(unique(x)) > 1L && length(unique(y)) > 1L) {
        correlation_test <- suppressWarnings(cor.test(
          x,
          y,
          method = "spearman",
          exact = FALSE
        ))
        spearman_rho <- unname(correlation_test$estimate)
        spearman_s <- unname(correlation_test$statistic)
        p_value <- correlation_test$p.value
      }

      result_index <- result_index + 1L
      results[[result_index]] <- data.frame(
        method = unname(method_names[method_suffix]),
        metric_1 = unname(metric_names[metric_pair[1L]]),
        metric_2 = unname(metric_names[metric_pair[2L]]),
        n_complete = n_complete,
        spearman_rho = spearman_rho,
        spearman_S = spearman_s,
        p_value = p_value,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  output <- do.call(rbind, results)
  row.names(output) <- NULL

  # Retain the unadjusted significance test and also report a Holm adjustment
  # across all 15 within-method correlations.
  output$p_value_holm <- p.adjust(output$p_value, method = "holm")
  output$significant_at_0.05 <- !is.na(output$p_value) & output$p_value < 0.05
  output$significant_holm_at_0.05 <-
    !is.na(output$p_value_holm) & output$p_value_holm < 0.05

  output
}

calculate_distance_correlations_by_experiment <- function(data) {
  method_suffixes <- c("i", "ii", "iii", "iv", "v")
  metric_names <- c(rmse = "RMSE", l1 = "L1", linf = "Linf")
  metric_pairs <- combn(names(metric_names), 2L, simplify = FALSE)

  mean_across_methods <- function(metric_prefix) {
    metric_columns <- paste0(metric_prefix, "_", method_suffixes)
    metric_matrix <- as.matrix(data[, metric_columns, drop = FALSE])

    apply(metric_matrix, 1L, function(values) {
      if (all(is.finite(values))) mean(values) else NA_real_
    })
  }

  # Each experiment is one independent unit. The method-specific values are
  # averaged within experiment before correlations are calculated.
  experiment_metrics <- data.frame(
    experiment = data$experiment,
    rmse = mean_across_methods("rmse"),
    l1 = mean_across_methods("l1"),
    linf = mean_across_methods("linf"),
    stringsAsFactors = FALSE
  )

  results <- list()

  for (pair_index in seq_along(metric_pairs)) {
    metric_pair <- metric_pairs[[pair_index]]
    x <- experiment_metrics[[metric_pair[1L]]]
    y <- experiment_metrics[[metric_pair[2L]]]
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    n_complete <- length(x)

    spearman_rho <- NA_real_
    spearman_s <- NA_real_
    p_value <- NA_real_

    if (n_complete >= 3L && length(unique(x)) > 1L && length(unique(y)) > 1L) {
      correlation_test <- suppressWarnings(cor.test(
        x,
        y,
        method = "spearman",
        exact = FALSE
      ))
      spearman_rho <- unname(correlation_test$estimate)
      spearman_s <- unname(correlation_test$statistic)
      p_value <- correlation_test$p.value
    }

    results[[pair_index]] <- data.frame(
      method = "All methods (experiment mean)",
      metric_1 = unname(metric_names[metric_pair[1L]]),
      metric_2 = unname(metric_names[metric_pair[2L]]),
      n_complete = n_complete,
      spearman_rho = spearman_rho,
      spearman_S = spearman_s,
      p_value = p_value,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  output <- do.call(rbind, results)
  row.names(output) <- NULL

  # These Holm adjustments apply to the three method-averaged correlations.
  output$p_value_holm <- p.adjust(output$p_value, method = "holm")
  output$significant_at_0.05 <- !is.na(output$p_value) & output$p_value < 0.05
  output$significant_holm_at_0.05 <-
    !is.na(output$p_value_holm) & output$p_value_holm < 0.05

  output
}

## ---- locate and parse experiment logs --------------------------------------
if (!is.numeric(n_experiments) || length(n_experiments) != 1L ||
    is.na(n_experiments) || n_experiments < 1 || n_experiments %% 1 != 0) {
  stop("n_experiments must be one positive integer.")
}

experiment_ids <- seq_len(as.integer(n_experiments))
expected_files <- file.path(
  root_dir,
  paste0("ex", experiment_ids),
  paste0("faces_ex", experiment_ids, ".lst")
)

file_exists <- file.exists(expected_files)
lst_files <- expected_files[file_exists]
missing_files <- expected_files[!file_exists]

cat(sprintf(
  "Found %d of %d expected files named ex#/faces_ex#.lst.\n",
  length(lst_files),
  length(expected_files)
))

if (length(missing_files) > 0L) {
  warning(sprintf(
    "%d expected .lst files are missing. First missing path(s): %s",
    length(missing_files),
    paste(head(missing_files, 10L), collapse = ", ")
  ))
}
if (length(lst_files) == 0L) {
  stop("No expected faces_ex#.lst files were found.")
}

parsed_results <- lapply(lst_files, function(path) {
  cat("Parsing:", path, "\n")
  tryCatch(
    parse_one_lst(path),
    error = function(e) {
      message("Failed on ", path, ": ", conditionMessage(e))
      NULL
    }
  )
})
parsed_results <- Filter(Negate(is.null), parsed_results)

if (length(parsed_results) == 0L) {
  stop("None of the located .lst files could be parsed.")
}

summary_df <- do.call(rbind, parsed_results)
summary_df <- summary_df[order(summary_df$experiment), , drop = FALSE]
row.names(summary_df) <- NULL

unknown_distribution <- is.na(summary_df$distribution_code)
if (any(unknown_distribution)) {
  warning(sprintf(
    "The sampling distribution could not be read for %d experiment(s).",
    sum(unknown_distribution)
  ))
}

distance_correlations <- calculate_distance_correlations(summary_df)
distance_correlations_by_experiment <-
  calculate_distance_correlations_by_experiment(summary_df)
distance_correlations <- rbind(
  distance_correlations,
  distance_correlations_by_experiment
)

write.csv(summary_df, "run6_summary.csv", row.names = FALSE, na = "")
write.csv(
  make_frequency_table(summary_df, "rmse_best"),
  "run6_rmse_best_frequency.csv",
  row.names = FALSE
)
write.csv(
  make_frequency_table(summary_df, "l1_best"),
  "run6_l1_best_frequency.csv",
  row.names = FALSE
)
write.csv(
  make_frequency_table(summary_df, "linf_best"),
  "run6_linf_best_frequency.csv",
  row.names = FALSE
)
write.csv(
  distance_correlations,
  "run6_distance_correlations.csv",
  row.names = FALSE,
  na = ""
)

cat("Wrote:\n")
cat("  run6_summary.csv\n")
cat("  run6_rmse_best_frequency.csv\n")
cat("  run6_l1_best_frequency.csv\n")
cat("  run6_linf_best_frequency.csv\n")
cat("  run6_distance_correlations.csv\n")
