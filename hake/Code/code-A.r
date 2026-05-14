# code-A.R
# Run with source("code-A.R")
#
# Assumption:
#   acc_mv already exists in the R session before this script is sourced.
#
# Ensure method is treated as a factor

acc_mv$method <- factor(acc_mv$method)

run_code_A <- function(data_name = "acc_mv", out_file = "results-code-A.txt") {

  # Open output file and make sure the sink is closed even if an error occurs.
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

    cat("Code A MANOVA and ANOVA results\n")
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

    # Treat method and design_block as categorical variables.
    acc_mv_model$method <- as.factor(acc_mv_model$method)
    acc_mv_model$design_block <- as.factor(acc_mv_model$design_block)

    response_cols <- c("rmse", "L1", "Linf")
    non_numeric_responses <- response_cols[!vapply(acc_mv_model[response_cols], is.numeric, logical(1))]
    if (length(non_numeric_responses) > 0) {
      stop("The following response column(s) must be numeric: ",
           paste(non_numeric_responses, collapse = ", "))
    }

    cat("Input checks\n")
    cat("------------\n")
    cat("Rows before complete-case filtering:", n_before, "\n")
    cat("Rows after complete-case filtering: ", n_after, "\n")
    cat("Number of methods:              ", nlevels(acc_mv_model$method), "\n")
    cat("Number of design blocks:        ", nlevels(acc_mv_model$design_block), "\n\n")

    # MANOVA uses rmse and L1 only because L1 and Linf may be collinear for
    # three-proportion accuracy errors. Linf is analyzed separately below.
    formula_manova <- cbind(rmse, L1) ~ method * design_block +
      method * G +
      method * theta_true +
      method * theta_CV +
      method * sigma +
      method * mean_nsamp +
      method * p1 + method * p2

    formula_linf <- Linf ~ method * design_block +
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

    cat("Analysis completed successfully.\n")

    invisible(list(
      manova = fit_manova,
      linf = fit_linf,
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
