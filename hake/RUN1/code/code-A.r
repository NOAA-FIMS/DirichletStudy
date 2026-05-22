# code-A.R
# Run with source("code-A.R")
#
# Assumption:
#   acc_mv already exists in the R session before this script is sourced.
#
# Purpose:
#   Conduct factorial MANOVA/ANOVA analyses of transformed accuracy metrics
#   and add formal post-hoc analyses using estimated marginal means (EMMs).
#
# Optional controls before source():
#   Sys.setenv(HAKE_RUN_CODE_A_POSTHOC = "true")
#   Sys.setenv(HAKE_POSTHOC_ADJUST = "holm")
#   Sys.setenv(HAKE_POSTHOC_BY_FACTORS = "design_block,G,sigma,mean_nsamp")

run_code_A <- function(data_name = "acc_mv", out_file = "results-code-A.txt") {

  # -----------------------------
  # Small utilities
  # -----------------------------
  truthy <- function(x) {
    tolower(as.character(x)) %in% c("true", "t", "yes", "y", "1")
  }

  add_emm_backtransform <- function(x) {
    x <- as.data.frame(x)
    estimate_col <- if ("emmean" %in% names(x)) {
      "emmean"
    } else if ("estimate" %in% names(x)) {
      "estimate"
    } else {
      NA_character_
    }

    if (!is.na(estimate_col)) {
      x$accuracy_scale_estimate <- expm1(x[[estimate_col]])
    }
    if ("lower.CL" %in% names(x)) {
      x$accuracy_scale_lower.CL <- expm1(x$lower.CL)
    }
    if ("upper.CL" %in% names(x)) {
      x$accuracy_scale_upper.CL <- expm1(x$upper.CL)
    }
    x
  }

  add_pairwise_ratio <- function(x) {
    x <- as.data.frame(x)
    if ("estimate" %in% names(x)) {
      # For a contrast on log1p(error), exp(estimate) is the ratio of
      # geometric means on the (1 + error) scale for the two methods.
      x$ratio_1p_error <- exp(x$estimate)
      x$percent_change_1p_error <- 100 * (x$ratio_1p_error - 1)
    }
    x
  }

  print_table <- function(x) {
    print(as.data.frame(x), row.names = FALSE)
  }

  run_emmeans_posthoc <- function(fit, metric_label, by_factors, adjust_method = "holm") {
    cat("\n")
    cat("Post-hoc estimated marginal means for ", metric_label, "\n", sep = "")
    cat(strrep("-", 72), "\n", sep = "")
    cat("Scale: model estimates are on the log1p(error) scale.\n")
    cat("Back-transformed estimates use expm1(estimate) to return to the original error scale.\n")
    cat("Pairwise contrasts are Holm-adjusted. Negative contrasts indicate lower error for the first method.\n")
    cat("For pairwise contrasts, ratio_1p_error = exp(estimate), the ratio on the (1 + error) scale.\n\n")

    emm_method <- emmeans::emmeans(
      fit,
      specs = ~ method,
      weights = "equal",
      cov.reduce = mean
    )

    cat("Marginal EMMs by method\n")
    cat("~~~~~~~~~~~~~~~~~~~~~~~\n")
    emm_method_sum <- summary(emm_method, infer = c(TRUE, TRUE))
    print_table(add_emm_backtransform(emm_method_sum))
    cat("\n")

    cat("Pairwise method comparisons\n")
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    method_pairs <- pairs(emm_method, adjust = adjust_method)
    method_pairs_sum <- summary(method_pairs, infer = c(TRUE, TRUE))
    print_table(add_pairwise_ratio(method_pairs_sum))
    cat("\n")

    if (length(by_factors) > 0) {
      cat("Selected method-by-factor post-hoc comparisons\n")
      cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
      cat("These comparisons diagnose whether method rankings change across important experimental conditions.\n\n")
    }

    for (by_factor in by_factors) {
      if (!by_factor %in% names(stats::model.frame(fit))) {
        cat("Skipping by-factor '", by_factor, "' because it is not in the fitted model frame.\n\n", sep = "")
        next
      }

      specs_by <- stats::as.formula(paste("~ method |", by_factor))

      emm_by <- emmeans::emmeans(
        fit,
        specs = specs_by,
        weights = "equal",
        cov.reduce = mean
      )

      cat("EMMs by method within ", by_factor, "\n", sep = "")
      cat(strrep("~", nchar(paste("EMMs by method within", by_factor))), "\n", sep = "")
      emm_by_sum <- summary(emm_by, infer = c(TRUE, TRUE))
      print_table(add_emm_backtransform(emm_by_sum))
      cat("\n")

      cat("Pairwise method comparisons within ", by_factor, "\n", sep = "")
      cat(strrep("~", nchar(paste("Pairwise method comparisons within", by_factor))), "\n", sep = "")
      by_pairs <- pairs(emm_by, adjust = adjust_method)
      by_pairs_sum <- summary(by_pairs, infer = c(TRUE, TRUE))
      print_table(add_pairwise_ratio(by_pairs_sum))
      cat("\n")
    }

    invisible(list(
      emm_method = emm_method,
      method_pairs = method_pairs
    ))
  }

  # -----------------------------
  # Output sink
  # -----------------------------
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

    cat("Code A MANOVA, ANOVA, and EMM post-hoc results\n")
    cat("Output file:", out_file, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

    # Check that acc_mv exists.
    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' was not found. Load or create acc_mv before sourcing this script.")
    }

    acc_mv_in <- get(data_name, envir = parent.frame(), inherits = TRUE)

    if (!is.data.frame(acc_mv_in)) {
      stop("Object '", data_name, "' must be a data.frame, data.table, or tibble-like object.")
    }

    required_cols <- c(
      "rmse", "L1", "Linf",
      "method", "design_block", "G",
      "theta_true", "theta_CV", "sigma", "mean_nsamp",
      "p1", "p2"
    )

    missing_cols <- setdiff(required_cols, names(acc_mv_in))
    if (length(missing_cols) > 0) {
      stop("The following required column(s) are missing from '", data_name, "': ",
           paste(missing_cols, collapse = ", "))
    }

    # Keep only needed columns and drop rows with missing values in the model variables.
    acc_mv_model <- acc_mv_in[, required_cols, drop = FALSE]
    n_before <- nrow(acc_mv_model)
    acc_mv_model <- acc_mv_model[stats::complete.cases(acc_mv_model), , drop = FALSE]
    n_after <- nrow(acc_mv_model)

    if (n_after == 0) {
      stop("No complete rows are available after removing missing values in the model variables.")
    }

    # Treat factorial-design variables as categorical variables.
    factor_cols <- c(
      "method", "design_block", "G", "theta_true",
      "theta_CV", "sigma", "mean_nsamp"
    )
    for (fc in factor_cols) {
      acc_mv_model[[fc]] <- as.factor(acc_mv_model[[fc]])
    }

    # Transform nonnegative scalar accuracy metrics.
    # This improves normality and stabilizes variance for MANOVA/ANOVA.
    accuracy_cols <- c("rmse", "L1", "Linf")

    if (any(acc_mv_model[accuracy_cols] < 0, na.rm = TRUE)) {
      stop("Accuracy metrics rmse, L1, and Linf must be nonnegative before log1p transformation.")
    }

    acc_mv_model$log_rmse <- log1p(acc_mv_model$rmse)
    acc_mv_model$log_L1   <- log1p(acc_mv_model$L1)
    acc_mv_model$log_Linf <- log1p(acc_mv_model$Linf)

    cat("Input checks\n")
    cat("------------\n")
    cat("Rows before complete-case filtering:", n_before, "\n")
    cat("Rows after complete-case filtering: ", n_after, "\n")
    cat("Number of methods:              ", nlevels(acc_mv_model$method), "\n")
    cat("Number of design blocks:        ", nlevels(acc_mv_model$design_block), "\n")
    cat("Number of G levels:             ", nlevels(acc_mv_model$G), "\n")
    cat("Number of sigma levels:         ", nlevels(acc_mv_model$sigma), "\n")
    cat("Number of mean_nsamp levels:    ", nlevels(acc_mv_model$mean_nsamp), "\n\n")

    # MANOVA uses rmse and L1 only because L1 and Linf may be collinear for
    # three-proportion accuracy errors. Linf is analyzed separately below.
    formula_manova <- cbind(log_rmse, log_L1) ~ method * design_block +
      method * G +
      method * theta_true +
      method * theta_CV +
      method * sigma +
      method * mean_nsamp +
      method * p1 + method * p2

    formula_rmse <- log_rmse ~ method * design_block +
      method * G +
      method * theta_true +
      method * theta_CV +
      method * sigma +
      method * mean_nsamp +
      method * p1 + method * p2

    formula_l1 <- log_L1 ~ method * design_block +
      method * G +
      method * theta_true +
      method * theta_CV +
      method * sigma +
      method * mean_nsamp +
      method * p1 + method * p2

    formula_linf <- log_Linf ~ method * design_block +
      method * G +
      method * theta_true +
      method * theta_CV +
      method * sigma +
      method * mean_nsamp +
      method * p1 + method * p2

    cat("MANOVA model formula\n")
    cat("--------------------\n")
    print(formula_manova)
    cat("\n")

    fit_manova <- stats::manova(formula_manova, data = acc_mv_model)

    cat("MANOVA results: Pillai test\n")
    cat("---------------------------\n")
    print(summary(fit_manova, test = "Pillai"))
    cat("\n")

    cat("Univariate ANOVA results for MANOVA responses\n")
    cat("---------------------------------------------\n")
    print(summary.aov(fit_manova))
    cat("\n")

    cat("Fitting separate univariate linear models for EMM post-hoc analyses.\n\n")

    fit_rmse <- stats::lm(formula_rmse, data = acc_mv_model)
    fit_l1   <- stats::lm(formula_l1,   data = acc_mv_model)

    cat("Linf model formula\n")
    cat("------------------\n")
    print(formula_linf)
    cat("\n")

    fit_linf <- stats::lm(formula_linf, data = acc_mv_model)

    cat("Linf linear model summary\n")
    cat("-------------------------\n")
    print(summary(fit_linf))
    cat("\n")

    cat("Linf ANOVA table\n")
    cat("----------------\n")
    print(stats::anova(fit_linf))
    cat("\n")

    # -----------------------------
    # Formal post-hoc analysis using estimated marginal means
    # -----------------------------
    run_posthoc <- truthy(Sys.getenv("HAKE_RUN_CODE_A_POSTHOC", unset = "true"))
    adjust_method <- Sys.getenv("HAKE_POSTHOC_ADJUST", unset = "holm")

    by_factor_string <- Sys.getenv(
      "HAKE_POSTHOC_BY_FACTORS",
      unset = "design_block,G,sigma,mean_nsamp"
    )
    by_factors <- trimws(unlist(strsplit(by_factor_string, ",", fixed = TRUE)))
    by_factors <- by_factors[nzchar(by_factors)]

    if (run_posthoc) {
      if (!requireNamespace("emmeans", quietly = TRUE)) {
        stop(
          "The emmeans package is required for the Code A post-hoc analysis. ",
          "Install it with install.packages('emmeans'), then re-run this script."
        )
      }

      cat("Formal post-hoc analysis using estimated marginal means\n")
      cat("======================================================\n")
      cat("Post-hoc adjustment method:", adjust_method, "\n")
      cat("By-factor diagnostics:", paste(by_factors, collapse = ", "), "\n")
      cat("Factor averaging: equal weights over factorial-design factor levels.\n")
      cat("Numeric covariates p1 and p2 are held at their means through cov.reduce = mean.\n\n")

      posthoc_rmse <- run_emmeans_posthoc(
        fit = fit_rmse,
        metric_label = "log_rmse",
        by_factors = by_factors,
        adjust_method = adjust_method
      )

      posthoc_l1 <- run_emmeans_posthoc(
        fit = fit_l1,
        metric_label = "log_L1",
        by_factors = by_factors,
        adjust_method = adjust_method
      )

      posthoc_linf <- run_emmeans_posthoc(
        fit = fit_linf,
        metric_label = "log_Linf",
        by_factors = by_factors,
        adjust_method = adjust_method
      )
    } else {
      cat("Formal post-hoc analysis skipped because HAKE_RUN_CODE_A_POSTHOC is false.\n\n")
      posthoc_rmse <- NULL
      posthoc_l1 <- NULL
      posthoc_linf <- NULL
    }

    cat("Analysis completed successfully.\n")

    invisible(list(
      manova = fit_manova,
      rmse = fit_rmse,
      l1 = fit_l1,
      linf = fit_linf,
      posthoc_rmse = posthoc_rmse,
      posthoc_l1 = posthoc_l1,
      posthoc_linf = posthoc_linf,
      output_file = out_file
    ))

  }, error = function(e) {
    cat("ERROR\n")
    cat("-----\n")
    cat(conditionMessage(e), "\n")
    stop(e)
  })
}

# Run the analysis when this script is sourced.
run_code_A()
