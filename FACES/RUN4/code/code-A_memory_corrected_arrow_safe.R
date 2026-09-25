# code-A_memory_corrected_arrow_safe.R
# Memory-safe Analysis A for metadata_revised.R.
#
# Do not run make_acc_mv_from_long(acc_long) before sourcing this file.
# This script converts acc_long RDS partitions to acc_mv one partition at a time
# and retains only a reproducible sample before fitting the models.

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
    "HAKE_METADATA_OUTDIR",
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
      ". Set HAKE_METADATA_OUTDIR to the correct metadata_output directory."
    )
  }
  files
}

partition_to_acc_mv_code_A <- function(dat, quota, seed_offset) {
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

  if (nrow(wide) > quota) {
    set.seed(seed_offset)
    wide <- wide[sample.int(nrow(wide), quota)]
  }
  wide[]
}

build_acc_mv_sample_code_A <- function(max_rows, sample_seed) {
  files <- find_acc_long_parts_code_A()
  quota <- ceiling(max_rows / length(files) * 1.10)
  pieces <- vector("list", length(files))

  message(
    "Creating a memory-safe acc_mv sample from ", length(files),
    " acc_long partition(s)."
  )

  for (i in seq_along(files)) {
    message("Processing partition ", i, " of ", length(files), ": ", basename(files[i]))
    part <- readRDS(files[i])
    pieces[[i]] <- partition_to_acc_mv_code_A(
      part, quota = quota, seed_offset = sample_seed + i
    )
    rm(part)
    invisible(gc())
  }

  ans <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  rm(pieces)
  invisible(gc())

  if (nrow(ans) > max_rows) {
    set.seed(sample_seed)
    ans <- ans[sample.int(nrow(ans), max_rows)]
  }
  as.data.frame(ans)
}

run_code_A <- function(out_file = "results-code-A.txt") {
  max_rows <- suppressWarnings(as.integer(
    Sys.getenv("HAKE_CODE_A_MAX_ROWS", unset = "500000")
  ))
  if (is.na(max_rows) || max_rows < 1000L) max_rows <- 500000L

  sample_seed <- suppressWarnings(as.integer(
    Sys.getenv("HAKE_CODE_A_SEED", unset = "20260730")
  ))
  if (is.na(sample_seed)) sample_seed <- 20260730L

  run_posthoc <- truthy_code_A(
    Sys.getenv("HAKE_RUN_CODE_A_POSTHOC", unset = "true")
  )
  adjust_method <- Sys.getenv("HAKE_POSTHOC_ADJUST", unset = "holm")

  con <- file(out_file, open = "wt")
  sink_start <- sink.number(type = "output")
  sink(con, type = "output")
  on.exit({
    while (sink.number(type = "output") > sink_start) sink(type = "output")
    if (isOpen(con)) close(con)
  }, add = TRUE)

  tryCatch({
    cat("Code A memory-safe MANOVA, ANOVA, and EMM results\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    cat("Maximum sampled acc_mv rows:", max_rows, "\n")
    cat("Sampling seed:", sample_seed, "\n\n")

    acc_mv_model <- build_acc_mv_sample_code_A(max_rows, sample_seed)

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
        seed = sample_seed
      ),
      "code-A-model-results.rds",
      compress = "gzip"
    )

    cat("\nAnalysis completed successfully.\n")
    cat("Saved: code-A-model-results.rds\n")

    invisible(list(
      manova = fit_manova,
      rmse = fit_rmse,
      L1 = fit_l1,
      Linf = fit_linf,
      posthoc = posthoc,
      data = acc_mv_model
    ))
  }, error = function(e) {
    cat("\nERROR\n-----\n", conditionMessage(e), "\n", sep = "")
    stop(e)
  })
}

code_A_result <- run_code_A()
