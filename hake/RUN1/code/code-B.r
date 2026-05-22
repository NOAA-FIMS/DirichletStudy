# code-B.R
#
# Best practical fix for the LINPACK / too-large-model-matrix issue in code-B.R.
#
# Key changes:
#   1. Use acc_mv instead of acc_long.
#   2. Fit one mixed model per accuracy metric: rmse, L1, and Linf.
#   3. Use simplex_id as the only random effect: (1 | simplex_id).
#   4. Avoid lmerTest and default emmeans calculations, which can trigger very
#      large dense matrix operations on this data size.
#   5. Reconstruct simplex_id from example_id and mesh_id if it is not already
#      present in acc_mv.
#
# Usage after metadata preparation:
#   acc_long <- load_metadata_table("acc_long")
#   acc_mv   <- make_acc_mv_from_long(acc_long)
#   source("code-B.R")
#
# Optional controls before source():
#   Sys.setenv(HAKE_CODE_B_RUN_EMMEANS = "true")       # default false
#   Sys.setenv(HAKE_CODE_B_SAVE_MODELS = "true")       # default false
#   Sys.setenv(HAKE_CODE_B_MODEL_DIR = "code_B_models")
#   Sys.setenv(HAKE_CODE_B_KEEP_MODELS = "false")      # default false

suppressPackageStartupMessages({
  library(lme4)
})

truthy <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "yes", "y", "1")
}

safe_sink_close <- function(start_output_sink, con) {
  while (sink.number(type = "output") > start_output_sink) {
    sink(type = "output")
  }
  if (!is.null(con) && isOpen(con)) close(con)
}

print_section <- function(title, underline = "-") {
  cat("\n", title, "\n", sep = "")
  cat(paste(rep(underline, nchar(title)), collapse = ""), "\n", sep = "")
}

prepare_acc_mv_for_lmm <- function(acc_mv_in) {
  required_metric_cols <- c("rmse", "L1", "Linf")
  required_model_cols <- c(
    "method", "design_block", "G", "theta_true", "theta_CV",
    "sigma", "mean_nsamp", "p1", "p2"
  )

  missing_cols <- setdiff(c(required_metric_cols, required_model_cols), names(acc_mv_in))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in acc_mv: ", paste(missing_cols, collapse = ", "))
  }

  # Reconstruct simplex_id if make_acc_mv_from_long() did not retain it.
  if (!"simplex_id" %in% names(acc_mv_in)) {
    if (all(c("example_id", "mesh_id") %in% names(acc_mv_in))) {
      acc_mv_in$simplex_id <- paste(acc_mv_in$example_id, acc_mv_in$mesh_id, sep = "_")
    } else {
      stop(
        "acc_mv does not contain simplex_id, and it cannot be reconstructed because ",
        "example_id and/or mesh_id are missing."
      )
    }
  }

  keep_cols <- unique(c(required_metric_cols, required_model_cols, "simplex_id"))
  acc_mv_model <- acc_mv_in[, keep_cols, drop = FALSE]
  acc_mv_model <- acc_mv_model[stats::complete.cases(acc_mv_model), , drop = FALSE]

  # Treat the design variables as factorial controls, consistent with the
  # 3-level factorial design matrix.
  acc_mv_model$method <- factor(
    acc_mv_model$method,
    levels = c("i", "ii", "iii", "iv", "v")
  )
  acc_mv_model$design_block <- as.factor(acc_mv_model$design_block)
  acc_mv_model$G <- as.factor(acc_mv_model$G)
  acc_mv_model$theta_true <- as.factor(acc_mv_model$theta_true)
  acc_mv_model$theta_CV <- as.factor(acc_mv_model$theta_CV)
  acc_mv_model$sigma <- as.factor(acc_mv_model$sigma)
  acc_mv_model$mean_nsamp <- as.factor(acc_mv_model$mean_nsamp)
  acc_mv_model$simplex_id <- as.factor(acc_mv_model$simplex_id)

  if (any(is.na(acc_mv_model$method))) {
    bad_methods <- unique(as.character(acc_mv_in$method[is.na(acc_mv_model$method)]))
    stop(
      "Some method values were not recognized as i, ii, iii, iv, or v. ",
      "Unrecognized values: ", paste(bad_methods, collapse = ", ")
    )
  }

  acc_mv_model
}

fit_one_metric_lmm <- function(metric_name, data, run_emmeans = FALSE) {
  if (!metric_name %in% c("rmse", "L1", "Linf")) {
    stop("Unknown metric_name: ", metric_name)
  }

  if (any(data[[metric_name]] < 0, na.rm = TRUE)) {
    stop("Metric ", metric_name, " contains negative values; log1p() is not appropriate.")
  }

  rhs <- paste(
    "method * design_block",
    "method * G",
    "method * theta_true",
    "method * theta_CV",
    "method * sigma",
    "method * mean_nsamp",
    "method * p1",
    "method * p2",
    "(1 | simplex_id)",
    sep = " + "
  )

  formula_lmm <- stats::as.formula(paste0("log1p(", metric_name, ") ~ ", rhs))

  print_section(paste("Mixed model for", metric_name), "=")
  cat("Formula:\n")
  print(formula_lmm)
  cat("\n")

  cat("Fitting lme4::lmer model using ML and calc.derivs = FALSE...\n")
  fit <- lme4::lmer(
    formula_lmm,
    data = data,
    REML = FALSE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      calc.derivs = FALSE,
      check.rankX = "message+drop.cols",
      check.conv.singular = "message"
    )
  )
  cat("Model fitting complete.\n")

  print_section(paste("Model summary for", metric_name))
  print(summary(fit))

  print_section(paste("Sequential ANOVA table for", metric_name))
  cat(
    "Note: This is the lme4 sequential ANOVA table for fixed effects. ",
    "It avoids lmerTest Satterthwaite/Kenward-Roger calculations, which can ",
    "trigger large dense matrix operations for this dataset.\n\n",
    sep = ""
  )
  print(stats::anova(fit))

  print_section(paste("Random-effect variance for", metric_name))
  print(lme4::VarCorr(fit), comp = c("Variance", "Std.Dev."))

  emm <- NULL
  if (isTRUE(run_emmeans)) {
    print_section(paste("Estimated marginal means for", metric_name))
    if (!requireNamespace("emmeans", quietly = TRUE)) {
      cat("The emmeans package is not installed; skipping emmeans.\n")
    } else {
      cat(
        "Running emmeans::emmeans(fit, pairwise ~ method). ",
        "For very large models this may still be computationally expensive.\n\n",
        sep = ""
      )
      emm <- emmeans::emmeans(fit, pairwise ~ method, adjust = "holm")
      print(emm)
    }
  }

  list(fit = fit, emmeans = emm)
}

run_code_B <- function(data_name = "acc_mv",
                       out_file = "results-code-B.txt",
                       metrics = c("rmse", "L1", "Linf"),
                       run_emmeans = truthy(Sys.getenv("HAKE_CODE_B_RUN_EMMEANS", unset = "false")),
                       save_models = truthy(Sys.getenv("HAKE_CODE_B_SAVE_MODELS", unset = "false")),
                       model_dir = Sys.getenv("HAKE_CODE_B_MODEL_DIR", unset = "code_B_models"),
                       keep_models = truthy(Sys.getenv("HAKE_CODE_B_KEEP_MODELS", unset = "false"))) {

  con <- file(out_file, open = "wt")
  sink_start <- sink.number(type = "output")
  sink(con, type = "output")

  on.exit({
    safe_sink_close(sink_start, con)
  }, add = TRUE)

  results <- list()
  model_files <- character(0)

  tryCatch({
    cat("Code B revised mixed-model results\n")
    cat("Output file:", out_file, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    cat("Input data object:", data_name, "\n")
    cat("Metrics:", paste(metrics, collapse = ", "), "\n")
    cat("Run emmeans:", run_emmeans, "\n")
    cat("Save models:", save_models, "\n")
    cat("Keep models in returned object:", keep_models, "\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' was not found.")
    }

    acc_mv_in <- get(data_name, envir = parent.frame(), inherits = TRUE)
    acc_mv_model <- prepare_acc_mv_for_lmm(acc_mv_in)

    print_section("Input checks")
    cat("Rows analyzed:", nrow(acc_mv_model), "\n")
    cat("Method counts:\n")
    print(table(acc_mv_model$method, useNA = "ifany"))
    cat("Design block counts:\n")
    print(table(acc_mv_model$design_block, useNA = "ifany"))
    cat("Number of simplex_id levels:", nlevels(acc_mv_model$simplex_id), "\n")
    cat("Object size of model data:", format(utils::object.size(acc_mv_model), units = "GB"), "\n")

    if (isTRUE(save_models)) {
      if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    }

    for (metric_name in metrics) {
      metric_result <- tryCatch({
        fit_result <- fit_one_metric_lmm(
          metric_name = metric_name,
          data = acc_mv_model,
          run_emmeans = run_emmeans
        )

        if (isTRUE(save_models)) {
          model_file <- file.path(model_dir, paste0("code_B_lmm_", metric_name, ".rds"))
          saveRDS(fit_result$fit, model_file, compress = "gzip")
          model_files <- c(model_files, model_file)
          cat("\nSaved model to:", model_file, "\n")
        }

        if (isTRUE(keep_models)) {
          fit_result
        } else {
          # Return compact metadata only. This avoids keeping three very large
          # fitted lmer objects in memory at the same time.
          list(
            formula = formula(fit_result$fit),
            logLik = as.numeric(stats::logLik(fit_result$fit)),
            AIC = stats::AIC(fit_result$fit),
            BIC = stats::BIC(fit_result$fit),
            sigma = stats::sigma(fit_result$fit),
            is_singular = lme4::isSingular(fit_result$fit),
            emmeans = fit_result$emmeans
          )
        }
      }, error = function(e) {
        print_section(paste("ERROR while fitting", metric_name))
        cat(conditionMessage(e), "\n")
        list(error = conditionMessage(e))
      })

      results[[metric_name]] <- metric_result
      invisible(gc())
    }

    print_section("Analysis status")
    failed <- vapply(results, function(x) !is.null(x$error), logical(1))
    if (any(failed)) {
      cat("One or more metric models failed:\n")
      print(names(results)[failed])
    } else {
      cat("All requested metric models completed successfully.\n")
    }

    if (length(model_files)) {
      cat("Saved model files:\n")
      print(model_files)
    }

    invisible(list(
      results = results,
      output_file = out_file,
      model_files = model_files,
      run_emmeans = run_emmeans,
      save_models = save_models,
      keep_models = keep_models
    ))

  }, error = function(e) {
    cat("\nERROR\n")
    cat("-----\n")
    cat(conditionMessage(e), "\n")
    stop(e)
  })
}

# Execute the analysis when sourced.
code_B_result <- run_code_B()
