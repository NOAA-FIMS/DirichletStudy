# code-B-RE-simplex.R
# Multivariate mixed model accuracy analysis using simplex_id effect 

library(lme4)
library(lmerTest)
library(emmeans)
emm_options(rg.limit = 11000)

run_code_B <- function(data_name = "acc_long", out_file = "results-code-B-RE-simplex.txt") {

  # Open output file and ensure the sink is closed even if an error occurs.
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

    cat("Code B Linear Mixed-Effects Model results\n")
    cat("Output file:", out_file, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

    # Check that the data object exists
    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' was not found.")
    }

    acc_long_in <- get(data_name, envir = parent.frame(), inherits = TRUE)

    # Required columns for the long-format factorial design
    required_cols <- c(
      "accuracy", "method", "metric", "design_block", "G",
      "theta_true", "theta_CV", "sigma", "mean_nsamp",
      "p1", "p2", "example_id", "simplex_id", "dataset_id"
    )

    missing_cols <- setdiff(required_cols, names(acc_long_in))
    if (length(missing_cols) > 0) {
      stop("Missing columns in '", data_name, "': ", paste(missing_cols, collapse = ", "))
    }

    # Data Cleaning and Factor Conversion
    acc_long_model <- acc_long_in[, required_cols, drop = FALSE]
    acc_long_model <- acc_long_model[stats::complete.cases(acc_long_model), , drop = FALSE]
    
    acc_long_model$method       <- as.factor(acc_long_model$method)
    acc_long_model$metric       <- as.factor(acc_long_model$metric)
    acc_long_model$design_block <- as.factor(acc_long_model$design_block)

	# Transform nonnegative scalar accuracy metric.
	if (any(acc_long_model$accuracy < 0, na.rm = TRUE)) {
	stop("Accuracy must be nonnegative before log1p transformation.")
	}

	acc_long_model$log_accuracy <- log1p(acc_long_model$accuracy)

    cat("Input checks\n")
    cat("------------\n")
    cat("Total rows analyzed:       ", nrow(acc_long_model), "\n")
    cat("Levels of Method:          ", nlevels(acc_long_model$method), "\n")
    cat("Levels of Metric:          ", nlevels(acc_long_model$metric), "\n\n")

    # Define Formula
    # Using log1p as per original code-B requirement
formula_lmm <- log_accuracy ~ method * metric +
  method * design_block +
  method * G +
  method * theta_true +
  method * theta_CV +
  method * sigma +
  method * mean_nsamp +
  method * p1 + method * p2 +
  (1 | simplex_id)

    cat("LMM Model Formula\n")
    cat("-----------------\n")
    print(formula_lmm)
    cat("\n")

    # Fit the Model
    cat("Fitting Linear Mixed-Effects Model (lmer)...\n")
    fit_mv_lmm <- lmer(formula_lmm, data = acc_long_model)
    cat("Model fitting complete.\n\n")

    # ANOVA Results
    cat("Type III Analysis of Variance Table with Satterthwaite's method\n")
    cat("--------------------------------------------------------------\n")
    print(anova(fit_mv_lmm))
    cat("\n")

    # Estimated Marginal Means
    cat("Post-hoc: Estimated Marginal Means (Pairwise Comparisons by Metric)\n")
    cat("------------------------------------------------------------------\n")
    emm_metric <- emmeans(fit_mv_lmm, pairwise ~ method | metric, adjust = "holm")
    print(emm_metric)
    cat("\n")

    cat("Post-hoc: Marginal Means for Methods (Averaged over Metrics)\n")
    cat("-----------------------------------------------------------\n")
    emm_method <- emmeans(fit_mv_lmm, pairwise ~ method, adjust = "holm")
    print(emm_method)
    cat("\n")
    
    cat("How much variance was explained by the random effect")
    print(VarCorr(fit_mv_lmm), comp = c("Variance", "Std.Dev."))
    

    cat("Analysis completed successfully.\n")

    invisible(list(
      model = fit_mv_lmm,
      emmeans_metric = emm_metric,
      emmeans_method = emm_method,
      output_file = out_file
    ))

  }, error = function(e) {
    cat("ERROR\n")
    cat("-----\n")
    cat(conditionMessage(e), "\n")
    stop(e)
  })
}

# Execute the analysis
run_code_B()