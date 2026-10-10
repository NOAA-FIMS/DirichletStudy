# code-B.R
# 
# source("code-B.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

truthy_code_b <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "yes", "y", "1")
}

run_code_B <- function(
    data_name = "acc_long",
    out_file = "results-code-B.txt",
    model_file = "code-B-model-results.rds",
    emmeans_metric_file = "emmeans-code-B-by-metric.csv",
    emmeans_method_file = "emmeans-code-B-method.csv",
    emmeans_reference_file = "emmeans-code-B-reference-condition.csv",
    max_simplex = as.integer(Sys.getenv("CODE_B_MAX_SIMPLEX", "20000")),
    sample_seed = as.integer(Sys.getenv("CODE_B_SAMPLE_SEED", "11307")),
    save_model = truthy_code_b(Sys.getenv("CODE_B_SAVE_MODEL", "true"))) {

  if (is.na(max_simplex) || max_simplex < 100L) max_simplex <- 20000L
  if (is.na(sample_seed)) sample_seed <- 20260730L

  caller_env <- parent.frame()

  con <- file(out_file, open = "wt")
  sink_start <- sink.number(type = "output")
  sink(con, type = "output")

  on.exit({
    while (sink.number(type = "output") > sink_start) sink(type = "output")
    close(con)
  }, add = TRUE)

  tryCatch({
    cat("Code B: memory-safe simplex random-effects analysis\n")
    cat("===================================================\n")
    cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
    cat("Output file: ", out_file, "\n\n", sep = "")

    if (!exists(data_name, envir = caller_env, inherits = TRUE)) {
      stop("Required data object '", data_name,
           "' was not found. Source metadata_revised.R first.")
    }

    acc_long_in <- get(data_name, envir = caller_env, inherits = TRUE)

    required_cols <- c(
      "accuracy", "method", "metric", "design_block", "G",
      "theta_true", "theta_CV", "sigma", "mean_nsamp",
      "p1", "p2", "example_id", "simplex_id"
    )

    available_cols <- names(acc_long_in)
    missing_cols <- setdiff(required_cols, available_cols)
    if (length(missing_cols)) {
      stop("Missing required columns in '", data_name, "': ",
           paste(missing_cols, collapse = ", "))
    }

    is_arrow <- inherits(
      acc_long_in,
      c("Dataset", "FileSystemDataset", "UnionDataset", "InMemoryDataset",
        "ArrowTabular", "Table", "RecordBatchReader")
    )

    group_cols <- c(
      "simplex_id", "example_id", "method", "metric", "design_block", "G",
      "theta_true", "theta_CV", "sigma", "mean_nsamp", "p1", "p2"
    )

    cat("Data preparation\n")
    cat("----------------\n")

    if (is_arrow) {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("acc_long is an Arrow object, but package 'arrow' is unavailable.")
      }

      cat("Input type: disk-backed Arrow Dataset\n")
      cat("Reducing data in Arrow before collection...\n")

      # IMPORTANT: all filtering, transformation, grouping, and summarization
      # occur in the Arrow query. Only the reduced cell-level table is collected.
      acc_cells <- acc_long_in |>
        dplyr::select(dplyr::all_of(required_cols)) |>
        dplyr::filter(
          !is.na(accuracy), accuracy >= 0,
          !is.na(method), !is.na(metric), !is.na(design_block),
          !is.na(G), !is.na(theta_true), !is.na(theta_CV),
          !is.na(sigma), !is.na(mean_nsamp),
          !is.na(p1), !is.na(p2),
          !is.na(example_id), !is.na(simplex_id)
        ) |>
        dplyr::mutate(log_accuracy_row = log1p(accuracy)) |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
        dplyr::summarise(
          log_accuracy = mean(log_accuracy_row),
          n_obs = dplyr::n(),
          .groups = "drop"
        ) |>
        dplyr::collect()

    } else if (is.data.frame(acc_long_in)) {
      cat("Input type: in-memory data frame\n")
      cat("Reducing data before model fitting...\n")

      acc_cells <- acc_long_in |>
        dplyr::select(dplyr::all_of(required_cols)) |>
        dplyr::filter(
          dplyr::if_all(dplyr::everything(), ~ !is.na(.x)),
          accuracy >= 0
        ) |>
        dplyr::mutate(log_accuracy_row = log1p(accuracy)) |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
        dplyr::summarise(
          log_accuracy = mean(log_accuracy_row),
          n_obs = dplyr::n(),
          .groups = "drop"
        )
    } else {
      stop("Unsupported acc_long class: ", paste(class(acc_long_in), collapse = ", "))
    }

    acc_cells <- as.data.frame(acc_cells)
    if (!nrow(acc_cells)) stop("No complete nonnegative observations remained after filtering.")

    # Convert design variables to factors after collection. p1 and p2 remain numeric.
    factor_cols <- c(
      "simplex_id", "method", "metric", "design_block", "G",
      "theta_true", "theta_CV", "sigma", "mean_nsamp"
    )
    acc_cells[factor_cols] <- lapply(acc_cells[factor_cols], factor)
    acc_cells$example_id <- as.integer(acc_cells$example_id)
    acc_cells$p1 <- as.numeric(acc_cells$p1)
    acc_cells$p2 <- as.numeric(acc_cells$p2)
    acc_cells$n_obs <- as.numeric(acc_cells$n_obs)
    acc_cells$log_accuracy <- as.numeric(acc_cells$log_accuracy)

    n_cells_full <- nrow(acc_cells)
    simplex_ids <- levels(acc_cells$simplex_id)
    n_simplex_full <- length(simplex_ids)

    cat("Original observations represented: ",
        format(sum(acc_cells$n_obs), big.mark = ","), "\n", sep = "")
    cat("Reduced simplex-method-metric cells: ",
        format(n_cells_full, big.mark = ","), "\n", sep = "")
    cat("Unique simplex_id values: ",
        format(n_simplex_full, big.mark = ","), "\n", sep = "")

    # Adaptive safeguard for the lmer model matrix. Sampling is by complete
    # simplex_id, so all method x metric cells for a selected simplex are retained.
    sampled <- FALSE
    if (n_simplex_full > max_simplex) {
      set.seed(sample_seed)
      keep_ids <- sample(simplex_ids, size = max_simplex, replace = FALSE)
      acc_cells <- acc_cells[acc_cells$simplex_id %in% keep_ids, , drop = FALSE]
      acc_cells <- droplevels(acc_cells)
      sampled <- TRUE

      cat("\nMemory safeguard activated\n")
      cat("--------------------------\n")
      cat("The reduced data still contained more than ", max_simplex,
          " simplex groups.\n", sep = "")
      cat("A reproducible random sample of complete simplex groups was used.\n")
      cat("Sampling seed: ", sample_seed, "\n", sep = "")
      cat("Simplex groups fitted: ", nlevels(acc_cells$simplex_id), "\n", sep = "")
      cat("Cells fitted: ", format(nrow(acc_cells), big.mark = ","), "\n", sep = "")
      cat("Observations represented by fitted cells: ",
          format(sum(acc_cells$n_obs), big.mark = ","), "\n", sep = "")
    } else {
      cat("All reduced simplex groups will be fitted.\n")
    }

    invisible(gc())

    # Sum-to-zero contrasts make the Type III tests interpretable in the usual way.
    old_contrasts <- options("contrasts")
    on.exit(options(contrasts = old_contrasts$contrasts), add = TRUE)
    options(contrasts = c("contr.sum", "contr.poly"))

    formula_lmm <- log_accuracy ~ method * metric +
      method * design_block +
      method * G +
      method * theta_true +
      method * theta_CV +
      method * sigma +
      method * mean_nsamp +
      method * p1 + method * p2 +
      (1 | simplex_id)

    cat("\nModel specification\n")
    cat("-------------------\n")
    print(formula_lmm)
    cat("\nResponse: cell mean of log1p(accuracy)\n")
    cat("Prior weight: number of original observations in each cell\n")
    cat("Random effect: simplex_id intercept\n")

    cat("\nFitting model...\n")
    fit_mv_lmm <- lmerTest::lmer(
      formula_lmm,
      data = acc_cells,
      weights = n_obs,
      REML = TRUE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 200000),
        calc.derivs = FALSE,
        check.nobs.vs.rankZ = "ignore"
      )
    )
    cat("Model fitting complete.\n\n")

    cat("Model summary\n")
    cat("-------------\n")
    print(summary(fit_mv_lmm))
    cat("\n")

    cat("Type III analysis of variance (Satterthwaite)\n")
    cat("------------------------------------------------\n")
    anova_result <- anova(fit_mv_lmm, type = 3, ddf = "Satterthwaite")
    print(anova_result)
    cat("\n")

    # emmeans would otherwise construct the full Cartesian product of all
    # design-factor levels. Because each design factor interacts with method,
    # none can be declared a nuisance factor. Evaluate comparisons at one
    # documented representative design condition instead of increasing rg.limit.
    modal_level <- function(x) {
      tab <- table(x, useNA = "no")
      if (!length(tab)) stop("Cannot determine a representative level.")
      names(tab)[which.max(tab)]
    }

    emm_at <- list(
      design_block = modal_level(acc_cells$design_block),
      G             = modal_level(acc_cells$G),
      theta_true    = modal_level(acc_cells$theta_true),
      theta_CV      = modal_level(acc_cells$theta_CV),
      sigma         = modal_level(acc_cells$sigma),
      mean_nsamp    = modal_level(acc_cells$mean_nsamp),
      p1            = mean(acc_cells$p1, na.rm = TRUE),
      p2            = mean(acc_cells$p2, na.rm = TRUE)
    )

    cat("Estimated-marginal-means reference condition\n")
    cat("----------------------------------------------\n")
    cat("Factor covariates are fixed at their most frequent fitted level;\n")
    cat("p1 and p2 are fixed at their fitted-data means.\n")
    for (nm in names(emm_at)) {
      cat(sprintf("%-14s = %s\n", nm, format(emm_at[[nm]], digits = 8)))
    }
    cat("\n")

    # This grid has only method x metric rows (5 x 3 under the expected design),
    # rather than all combinations of the factorial-design predictors.
    emm_metric <- emmeans::emmeans(
      fit_mv_lmm,
      specs = ~ method | metric,
      at = emm_at,
      weights = "equal",
      rg.limit = 10000
    )
    pairs_metric <- emmeans::contrast(
      emm_metric, method = "pairwise", adjust = "holm"
    )

    cat("Estimated marginal means: method within metric\n")
    cat("------------------------------------------------\n")
    print(emm_metric)
    print(pairs_metric)
    cat("\n")

    # Construct a second small grid using the same fixed reference condition;
    # metric is averaged equally while all design factors remain fixed.
    emm_method <- emmeans::emmeans(
      fit_mv_lmm,
      specs = ~ method,
      at = emm_at,
      weights = "equal",
      rg.limit = 10000
    )
    pairs_method <- emmeans::contrast(
      emm_method, method = "pairwise", adjust = "holm"
    )

    cat("Estimated marginal means: method averaged equally over metric\n")
    cat("--------------------------------------------------------------\n")
    print(emm_method)
    print(pairs_method)
    cat("\n")

    # Write compact machine-readable outputs.
    metric_out <- as.data.frame(emm_metric)
    method_out <- as.data.frame(emm_method)
    utils::write.csv(metric_out, emmeans_metric_file, row.names = FALSE)
    utils::write.csv(method_out, emmeans_method_file, row.names = FALSE)
    reference_out <- data.frame(
      variable = names(emm_at),
      value = vapply(emm_at, function(z) paste(z, collapse = ","), character(1)),
      stringsAsFactors = FALSE
    )
    utils::write.csv(reference_out, emmeans_reference_file, row.names = FALSE)

    if (save_model) {
      saveRDS(fit_mv_lmm, model_file, compress = FALSE)
      cat("Saved fitted model: ", model_file, "\n", sep = "")
    }
    cat("Saved metric-specific EMMs: ", emmeans_metric_file, "\n", sep = "")
    cat("Saved overall method EMMs: ", emmeans_method_file, "\n", sep = "")
    cat("Saved EMM reference condition: ", emmeans_reference_file, "\n", sep = "")

    cat("\nConvergence and singularity checks\n")
    cat("----------------------------------\n")
    cat("Singular fit: ", lme4::isSingular(fit_mv_lmm, tol = 1e-4), "\n", sep = "")
    opt_messages <- fit_mv_lmm@optinfo$conv$lme4$messages
    if (is.null(opt_messages)) {
      cat("Optimizer messages: none\n")
    } else {
      cat("Optimizer messages:\n")
      print(opt_messages)
    }

    cat("\nAnalysis completed successfully.\n")
    if (sampled) {
      cat("NOTE: The fitted model used a reproducible sample of simplex groups because\n")
      cat("the complete reduced model exceeded the configured memory safeguard.\n")
    }

    invisible(list(
      model = fit_mv_lmm,
      anova = anova_result,
      emmeans_metric = emm_metric,
      pairs_metric = pairs_metric,
      emmeans_method = emm_method,
      pairs_method = pairs_method,
      sampled = sampled,
      n_simplex_available = n_simplex_full,
      n_simplex_fitted = nlevels(acc_cells$simplex_id),
      n_cells_fitted = nrow(acc_cells),
      output_file = out_file,
      model_file = if (save_model) model_file else NULL
    ))

  }, error = function(e) {
    cat("\nERROR\n")
    cat("-----\n")
    cat(conditionMessage(e), "\n")
    stop(e)
  })
}

# Execute when sourced.
code_B_result <- run_code_B()
