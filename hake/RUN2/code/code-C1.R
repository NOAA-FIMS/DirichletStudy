# code-C1.R
# Pairwise one-sample Hotelling's T^2 / MANOVA tests for method differences.
#
# Purpose:
#   For each pair of methods, compute the vector of log accuracy differences
#   for rmse, L1, and Linf, then test whether the mean difference vector is
#   zero. The test is the one-sample Hotelling's T^2 test, which is equivalent
#   to a MANOVA test of the intercept vector against zero.
#
# Interpretation:
#   diff = log(accuracy_target + eps) - log(accuracy_reference + eps)
#   Negative mean differences indicate lower error for the target method.
#   Positive mean differences indicate higher error for the target method.

library(dplyr)
library(tidyr)

run_code_C1 <- function(data_name = "acc_long",
                        out_file = "results-code-C1.txt",
                        eps = 1e-6) {

  con <- file(out_file, open = "wt")
  sink_start <- sink.number(type = "output")
  sink(con, type = "output")

  on.exit({
    while (sink.number(type = "output") > sink_start) {
      sink(type = "output")
    }
    close(con)
  }, add = TRUE)

  tryCatch({
    cat("Code C1: Pairwise Hotelling's T^2 / MANOVA results\n")
    cat("Output file:", out_file, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    cat("Transformation: log(accuracy + eps), eps =", eps, "\n")
    cat("Contrast direction: target method minus reference method\n")
    cat("Negative log-differences imply lower error for the target method.\n\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' not found.")
    }

    acc_long <- get(data_name, envir = parent.frame(), inherits = TRUE)

    required_cols <- c(
      "example_id", "design_block", "G",
      "theta_true", "theta_CV", "sigma", "mean_nsamp",
      "p1", "p2", "method", "metric", "accuracy"
    )

    missing_cols <- setdiff(required_cols, names(acc_long))
    if (length(missing_cols) > 0L) {
      stop(
        "The input data object is missing required column(s): ",
        paste(missing_cols, collapse = ", ")
      )
    }

    make_pairwise_diff <- function(dat, m_target, m_ref, eps = 1e-6) {

      id_vars <- c(
        "example_id", "design_block", "G",
        "theta_true", "theta_CV", "sigma", "mean_nsamp",
        "p1", "p2"
      )

      out <- dat %>%
        select(
          all_of(id_vars),
          method, metric, accuracy
        ) %>%
        filter(method %in% c(m_target, m_ref)) %>%
        mutate(
          accuracy = as.numeric(accuracy),

          # Standardize metric names so cbind(rmse, L1, Linf) is valid.
          metric = case_when(
            tolower(as.character(metric)) %in% c("rmse", "rmse_norm") ~ "rmse",
            tolower(as.character(metric)) %in% c("l1", "l1_norm", "l1 norm") ~ "L1",
            tolower(as.character(metric)) %in% c(
              "linf", "l_inf", "linf_norm", "l_inf_norm", "linf norm"
            ) ~ "Linf",
            TRUE ~ as.character(metric)
          ),

          accuracy_log = log(accuracy + eps)
        ) %>%
        filter(
          is.finite(accuracy_log),
          metric %in% c("rmse", "L1", "Linf")
        ) %>%

        # Collapse duplicate cells before widening.
        group_by(across(all_of(c(id_vars, "metric", "method")))) %>%
        summarise(
          accuracy_log = mean(accuracy_log, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        pivot_wider(
          id_cols = c(all_of(id_vars), metric),
          names_from = method,
          values_from = accuracy_log
        ) %>%
        filter(
          !is.na(.data[[m_target]]),
          !is.na(.data[[m_ref]])
        ) %>%
        mutate(
          diff = .data[[m_target]] - .data[[m_ref]]
        ) %>%
        select(
          all_of(id_vars),
          metric,
          diff
        ) %>%
        pivot_wider(
          id_cols = all_of(id_vars),
          names_from = metric,
          values_from = diff
        )

      metric_cols <- c("rmse", "L1", "Linf")
      missing_metrics <- setdiff(metric_cols, names(out))
      if (length(missing_metrics) > 0L) {
        stop(
          "Pairwise difference table for ", m_target, "_vs_", m_ref,
          " is missing metric column(s): ",
          paste(missing_metrics, collapse = ", ")
        )
      }

      out %>%
        filter(if_all(all_of(metric_cols), ~ !is.na(.x)))
    }

    safe_solve <- function(a, b) {
      out <- tryCatch(
        solve(a, b),
        error = function(e) NULL
      )

      if (!is.null(out)) {
        return(out)
      }

      # Fallback for near-singular covariance matrices.
      qr.solve(a, b)
    }

    one_sample_hotelling <- function(diff_mv,
                                     response_cols = c("rmse", "L1", "Linf")) {

      y <- as.matrix(diff_mv[, response_cols, drop = FALSE])
      storage.mode(y) <- "double"
      y <- y[stats::complete.cases(y), , drop = FALSE]

      n <- nrow(y)
      p <- ncol(y)

      if (n <= p) {
        stop(
          "Hotelling's T^2 requires more complete rows than response variables. ",
          "Complete rows = ", n, "; response variables = ", p, "."
        )
      }

      ybar <- colMeans(y)
      s <- stats::cov(y)

      solved <- safe_solve(s, ybar)
      t2 <- as.numeric(n * crossprod(ybar, solved))

      # One-sample Hotelling's T^2 F approximation.
      num_df <- p
      den_df <- n - p
      approx_f <- ((n - p) / (p * (n - 1))) * t2
      p_value <- stats::pf(approx_f, df1 = num_df, df2 = den_df, lower.tail = FALSE)

      # Pillai trace for the one-sample intercept test.
      # H = n * ybar %*% t(ybar), E = (n - 1) * cov(y).
      # For rank-one H, Pillai = T2 / (T2 + n - 1).
      pillai <- t2 / (t2 + n - 1)

      manova_table <- data.frame(
        Df = 1L,
        Pillai = pillai,
        `approx F` = approx_f,
        `num Df` = num_df,
        `den Df` = den_df,
        `Pr(>F)` = p_value,
        check.names = FALSE
      )

      rownames(manova_table) <- "Intercept"

      mean_table <- data.frame(
        metric = response_cols,
        mean_log_difference = as.numeric(ybar),
        error_ratio = exp(as.numeric(ybar)),
        percent_error_change = 100 * (exp(as.numeric(ybar)) - 1),
        check.names = FALSE
      )

      list(
        n = n,
        p = p,
        T2 = t2,
        manova_table = manova_table,
        mean_table = mean_table
      )
    }

    methods <- c("i", "ii", "iii", "iv", "v")
    pairs <- combn(methods, 2, simplify = FALSE)

    for (pair in pairs) {
      m_ref <- pair[1]
      m_target <- pair[2]
      comp_name <- paste0(m_target, "_vs_", m_ref)

      cat("==================================================================\n")
      cat("Comparison:", comp_name, "\n")
      cat("==================================================================\n")

      diff_mv <- make_pairwise_diff(acc_long, m_target, m_ref, eps = eps)

      # Formula specification for the equivalent intercept-only MANOVA:
      #   cbind(rmse, L1, Linf) ~ 1
      # The default summary(manova(... ~ 1)) prints only residuals, so the
      # intercept test H0: mean difference vector = 0 is computed explicitly.
      cat("Equivalent MANOVA formula: cbind(rmse, L1, Linf) ~ 1\n")
      cat("Null hypothesis: mean log-difference vector = (0, 0, 0)\n\n")

      hotelling <- one_sample_hotelling(
        diff_mv = diff_mv,
        response_cols = c("rmse", "L1", "Linf")
      )

      cat("Pairwise one-sample Hotelling's T^2 / Pillai test\n")
      cat("Complete paired rows:", hotelling$n, "\n")
      cat("Response dimension:", hotelling$p, "\n")
      cat("Hotelling T^2:", format(hotelling$T2, scientific = TRUE, digits = 8), "\n\n")

      print(hotelling$manova_table, digits = 8)
      cat("\n")

      cat("Mean log-differences and relative error changes\n")
      print(hotelling$mean_table, digits = 8, row.names = FALSE)
      cat("\n")
    }

    cat("Analysis completed successfully.\n")
  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n")
    stop(e)
  })
}

run_code_C1()
