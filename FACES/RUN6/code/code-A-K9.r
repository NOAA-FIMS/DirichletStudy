# code-A-K9.r
#
# Randomly subsample 2% of unique simplex estimation-accuracy samples,
# retain all methods and accuracy metrics for each sampled unit, and run
# the Code A MANOVA, ANOVA, and estimated-marginal-means analyses.
#
# source("code-A-K9.r")

suppressPackageStartupMessages({
  library(data.table)
})

truthy_code_A <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "yes", "y", "1")
}

mode_value_code_A <- function(x) {
  x <- x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

find_acc_long_parts_code_A <- function() {
  out_dir <- Sys.getenv(
    "METADATA_OUTDIR",
    unset = if (exists("OUT_DIR", inherits = TRUE)) get("OUT_DIR", inherits = TRUE)
    else file.path(getwd(), "metadata_output")
  )
  part_dir <- file.path(out_dir, "acc_long_parts")
  files <- sort(list.files(
    part_dir,
    pattern = "^acc_long_part_[0-9]+\\.rds$",
    full.names = TRUE
  ))
  if (!length(files)) {
    stop(
      "No acc_long RDS partitions found in ", part_dir,
      ". Set METADATA_OUTDIR to the correct metadata_output directory."
    )
  }
  files
}

partition_to_acc_mv_code_A <- function(dat, sample_fraction, seed_offset) {
  required_long <- c(
    "accuracy", "metric", "method", "example_id", "mesh_id", "sim_id",
    "p1", "p2", "design_block", "G", "theta_true", "theta_CV",
    "sigma", "mean_nsamp"
  )
  missing_cols <- setdiff(required_long, names(dat))
  if (length(missing_cols)) {
    stop("acc_long partition missing: ", paste(missing_cols, collapse = ", "))
  }

  dt <- as.data.table(dat)[, ..required_long]
  dt <- dt[metric %chin% c("rmse", "L1_norm", "Linf_norm") & !is.na(accuracy)]

  id_cols <- setdiff(required_long, c("accuracy", "metric"))
  cast_formula <- as.formula(paste(paste(id_cols, collapse = " + "), "~ metric"))

  wide <- dcast(
    dt,
    formula = cast_formula,
    value.var = "accuracy",
    fun.aggregate = function(x) x[1L],
    fill = NA_real_
  )

  rm(dt, dat)
  invisible(gc())

  if ("L1_norm" %in% names(wide)) setnames(wide, "L1_norm", "L1")
  if ("Linf_norm" %in% names(wide)) setnames(wide, "Linf_norm", "Linf")

  required_wide <- c(
    "rmse", "L1", "Linf", "method", "design_block", "G",
    "theta_true", "theta_CV", "sigma", "mean_nsamp", "p1", "p2"
  )
  missing_wide <- setdiff(required_wide, names(wide))
  if (length(missing_wide)) {
    stop("Could not create acc_mv; missing: ", paste(missing_wide, collapse = ", "))
  }

  wide <- wide[complete.cases(wide[, ..required_wide])]

  # The sampling unit is a simulated simplex sample, not an individual method
  # result. Sampling units first preserves all method results for paired
  # comparisons on every selected simplex sample.
  sample_unit_cols <- c("example_id", "mesh_id", "sim_id")
  sample_units <- unique(wide[, ..sample_unit_cols])
  n_units_available <- nrow(sample_units)
  n_units_sampled <- if (n_units_available > 0L) {
    max(1L, as.integer(ceiling(sample_fraction * n_units_available)))
  } else {
    0L
  }

  if (n_units_sampled > 0L && n_units_sampled < n_units_available) {
    set.seed(seed_offset)
    sample_units <- sample_units[sample.int(n_units_available, n_units_sampled)]
  }

  sampled <- wide[sample_units, on = sample_unit_cols, nomatch = 0L][]
  attr(sampled, "n_units_available") <- n_units_available
  attr(sampled, "n_units_sampled") <- n_units_sampled
  sampled
}

build_acc_mv_sample_code_A <- function(sample_fraction, sample_seed) {
  files <- find_acc_long_parts_code_A()
  pieces <- vector("list", length(files))
  sampling_counts <- data.table(
    partition = basename(files),
    units_available = integer(length(files)),
    units_sampled = integer(length(files))
  )

  message(
    "Creating a memory-safe ", 100 * sample_fraction,
    "% simplex sample from ", length(files), " acc_long partition(s)."
  )

  for (i in seq_along(files)) {
    message("Processing partition ", i, " of ", length(files), ": ", basename(files[i]))
    part <- readRDS(files[i])
    piece <- partition_to_acc_mv_code_A(
      part,
      sample_fraction = sample_fraction,
      seed_offset = sample_seed + i
    )
    sampling_counts$units_available[i] <- attr(piece, "n_units_available")
    sampling_counts$units_sampled[i] <- attr(piece, "n_units_sampled")
    pieces[[i]] <- piece
    rm(part, piece)
    invisible(gc())
  }

  ans <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  rm(pieces)
  invisible(gc())

  ans <- as.data.frame(ans)
  attr(ans, "sampling_counts") <- sampling_counts
  ans
}

run_code_A <- function(out_file = "results-code-A-K9.txt") {
  sample_fraction <- suppressWarnings(as.numeric(
    Sys.getenv("CODE_A_SAMPLE_FRACTION", unset = "0.02")
  ))
  if (!is.finite(sample_fraction) || sample_fraction <= 0 || sample_fraction > 1) {
    stop("CODE_A_SAMPLE_FRACTION must be greater than 0 and no greater than 1.")
  }

  sample_seed <- suppressWarnings(as.integer(
    Sys.getenv("CODE_A_SEED", unset = "20260730")
  ))
  if (is.na(sample_seed)) sample_seed <- 20260730L

  run_posthoc <- truthy_code_A(
    Sys.getenv("RUN_CODE_A_POSTHOC", unset = "true")
  )
  adjust_method <- Sys.getenv("POSTHOC_ADJUST", unset = "holm")

  con <- file(out_file, open = "wt")
  sink_start <- sink.number(type = "output")
  sink(con, type = "output")
  on.exit({
    while (sink.number(type = "output") > sink_start) sink(type = "output")
    if (isOpen(con)) close(con)
  }, add = TRUE)

  tryCatch({
    cat("Code A simplex-subsampled MANOVA, ANOVA, and EMM results\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    cat("Simplex sampling fraction:", sample_fraction,
        "(", 100 * sample_fraction, "%)\n", sep = "")
    cat("Sampling unit: example_id x mesh_id x sim_id\n")
    cat("Sampling seed:", sample_seed, "\n\n")

    acc_mv_model <- build_acc_mv_sample_code_A(sample_fraction, sample_seed)
    sampling_counts <- attr(acc_mv_model, "sampling_counts")

    cat("Simplex sampling counts by partition\n")
    print(sampling_counts)
    cat("\nTotal simplex units available:", sum(sampling_counts$units_available), "\n")
    cat("Total simplex units sampled:", sum(sampling_counts$units_sampled), "\n")
    cat(
      "Realized simplex sampling percentage:",
      100 * sum(sampling_counts$units_sampled) /
        sum(sampling_counts$units_available),
      "\n\n"
    )

    required_cols <- c(
      "rmse", "L1", "Linf", "method", "design_block", "G",
      "theta_true", "theta_CV", "sigma", "mean_nsamp", "p1", "p2"
    )
    acc_mv_model <- acc_mv_model[, required_cols, drop = FALSE]
    acc_mv_model <- acc_mv_model[complete.cases(acc_mv_model), , drop = FALSE]

    factor_cols <- c(
      "method", "design_block", "G", "theta_true",
      "theta_CV", "sigma", "mean_nsamp"
    )
    for (fc in factor_cols) {
      acc_mv_model[[fc]] <- droplevels(as.factor(acc_mv_model[[fc]]))
    }

    if (any(acc_mv_model[c("rmse", "L1", "Linf")] < 0, na.rm = TRUE)) {
      stop("Accuracy metrics must be nonnegative.")
    }

    acc_mv_model$log_rmse <- log1p(acc_mv_model$rmse)
    acc_mv_model$log_L1 <- log1p(acc_mv_model$L1)
    acc_mv_model$log_Linf <- log1p(acc_mv_model$Linf)

    cat("Rows analyzed:", nrow(acc_mv_model), "\n\n")

    rhs <- paste(
      "method * design_block", "method * G", "method * theta_true",
      "method * theta_CV", "method * sigma", "method * mean_nsamp",
      "method * p1", "method * p2", sep = " + "
    )

    formula_manova <- as.formula(paste("cbind(log_rmse, log_L1) ~", rhs))
    formula_rmse <- as.formula(paste("log_rmse ~", rhs))
    formula_l1 <- as.formula(paste("log_L1 ~", rhs))
    formula_linf <- as.formula(paste("log_Linf ~", rhs))

    fit_manova <- manova(formula_manova, data = acc_mv_model)
    cat("MANOVA Pillai test\n")
    print(summary(fit_manova, test = "Pillai"))
    cat("\nUnivariate ANOVA results\n")
    print(summary.aov(fit_manova))

    fit_rmse <- lm(formula_rmse, data = acc_mv_model)
    fit_l1 <- lm(formula_l1, data = acc_mv_model)
    fit_linf <- lm(formula_linf, data = acc_mv_model)

    cat("\nLinf model summary\n")
    print(summary(fit_linf))
    cat("\nLinf ANOVA table\n")
    print(anova(fit_linf))

    posthoc <- NULL
    if (run_posthoc) {
      if (!requireNamespace("emmeans", quietly = TRUE)) {
        stop("Install emmeans to run post-hoc analyses.")
      }

      representative_at <- list(
        design_block = mode_value_code_A(acc_mv_model$design_block),
        G = mode_value_code_A(acc_mv_model$G),
        theta_true = mode_value_code_A(acc_mv_model$theta_true),
        theta_CV = mode_value_code_A(acc_mv_model$theta_CV),
        sigma = mode_value_code_A(acc_mv_model$sigma),
        mean_nsamp = mode_value_code_A(acc_mv_model$mean_nsamp),
        p1 = mean(acc_mv_model$p1),
        p2 = mean(acc_mv_model$p2)
      )

      cat("\nPost-hoc reference condition\n")
      print(representative_at)

      run_one_emm <- function(fit, label) {
        emm <- emmeans::emmeans(
          fit, ~ method,
          at = representative_at,
          weights = "equal",
          rg.limit = 10000
        )
        cmp <- pairs(emm, adjust = adjust_method)
        cat("\nEstimated marginal means:", label, "\n")
        print(summary(emm, infer = c(TRUE, TRUE)))
        cat("\nPairwise comparisons:", label, "\n")
        print(summary(cmp, infer = c(TRUE, TRUE)))
        list(emmeans = emm, pairs = cmp)
      }

      posthoc <- list(
        rmse = run_one_emm(fit_rmse, "log_rmse"),
        L1 = run_one_emm(fit_l1, "log_L1"),
        Linf = run_one_emm(fit_linf, "log_Linf")
      )
    }

    saveRDS(
      list(
        manova = fit_manova,
        rmse = fit_rmse,
        L1 = fit_l1,
        Linf = fit_linf,
        posthoc = posthoc,
        rows_analyzed = nrow(acc_mv_model),
        sample_fraction = sample_fraction,
        sampling_counts = sampling_counts,
        seed = sample_seed
      ),
      "code-A-K9-model-results.rds",
      compress = "gzip"
    )

    cat("\nAnalysis completed successfully.\n")
    cat("Saved: code-A-K9-model-results.rds\n")

    invisible(list(
      manova = fit_manova,
      rmse = fit_rmse,
      L1 = fit_l1,
      Linf = fit_linf,
      posthoc = posthoc,
      sampling_counts = sampling_counts,
      data = acc_mv_model
    ))
  }, error = function(e) {
    cat("\nERROR\n-----\n", conditionMessage(e), "\n", sep = "")
    stop(e)
  })
}

code_A_result <- run_code_A()
