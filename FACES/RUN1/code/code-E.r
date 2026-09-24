# code-E.R
#
# source("code-E.R")
#
# Optional controls:
#   Sys.setenv(CODE_E_MAX_ROWS = "500000")
#   Sys.setenv(CODE_E_MAX_DIFF_ROWS = "300000")
#   Sys.setenv(CODE_E_SEED = "20731")

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(data.table)
})

find_acc_long_parts_E <- function() {
  out_dir <- Sys.getenv(
    "METADATA_OUTDIR",
    unset = if (exists("OUT_DIR", inherits = TRUE)) get("OUT_DIR", inherits = TRUE)
    else file.path(getwd(), "metadata_output")
  )
  part_dir <- file.path(out_dir, "acc_long_parts")
  files <- sort(list.files(part_dir,
                           pattern = "^acc_long_part_[0-9]+\\.rds$",
                           full.names = TRUE))
  if (!length(files)) {
    stop("No acc_long RDS partitions found in: ", part_dir,
         "\nRun metadata_revised.R first or set METADATA_OUTDIR correctly.")
  }
  files
}

safe_label_E <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

standardize_metric_E <- function(x) {
  z <- tolower(as.character(x))
  fifelse(z %chin% c("rmse", "rmse_norm"), "rmse",
    fifelse(z %chin% c("l1", "l1_norm", "l1 norm"), "L1",
      fifelse(z %chin% c("linf", "l_inf", "linf_norm", "l_inf_norm", "linf norm"),
              "Linf", as.character(x))))
}

prepare_global_partition_E <- function(dat, quota, seed) {
  required_cols <- c("example_id", "design_block", "G", "theta_true", "theta_CV",
                     "sigma", "mean_nsamp", "p1", "p2", "method", "metric", "accuracy")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

  dt <- as.data.table(dat)[, ..required_cols]
  dt[, method := as.character(method)]
  dt[, metric := standardize_metric_E(metric)]
  dt[, accuracy := as.numeric(accuracy)]
  for (nm in c("G", "theta_true", "theta_CV", "sigma", "mean_nsamp", "p1", "p2")) {
    set(dt, j = nm, value = as.numeric(dt[[nm]]))
  }
  dt <- dt[is.finite(accuracy) & accuracy >= 0 & is.finite(p1) & is.finite(p2) &
             metric %chin% c("rmse", "L1", "Linf") &
             method %chin% c("i", "ii", "iii", "iv", "v")]
  if (nrow(dt) > quota) {
    set.seed(seed)
    dt <- dt[sample.int(nrow(dt), quota, replace = FALSE)]
  }
  dt[]
}

prepare_diff_partition_E <- function(dat, quota, seed) {
  id_vars <- c("example_id", "design_block", "G", "theta_true", "theta_CV",
               "sigma", "mean_nsamp", "p1", "p2", "metric")
  required_cols <- c(id_vars, "method", "accuracy")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

  dt <- as.data.table(dat)[, ..required_cols]
  dt[, method := as.character(method)]
  dt <- dt[method %chin% c("v", "iii")]
  if (!nrow(dt)) return(data.table())

  dt[, metric := standardize_metric_E(metric)]
  dt[, accuracy := as.numeric(accuracy)]
  for (nm in c("G", "theta_true", "theta_CV", "sigma", "mean_nsamp", "p1", "p2")) {
    set(dt, j = nm, value = as.numeric(dt[[nm]]))
  }
  dt <- dt[is.finite(accuracy) & is.finite(p1) & is.finite(p2) &
             metric %chin% c("rmse", "L1", "Linf")]
  if (!nrow(dt)) return(data.table())

  dt <- dt[, .(accuracy = mean(accuracy, na.rm = TRUE)), by = c(id_vars, "method")]
  cast_formula <- as.formula(paste(paste(id_vars, collapse = " + "), "~ method"))
  wide <- dcast(dt, formula = cast_formula, value.var = "accuracy", fill = NA_real_)
  rm(dt); invisible(gc())
  if (!all(c("v", "iii") %in% names(wide))) return(data.table())
  wide <- wide[!is.na(v) & !is.na(iii)]
  wide[, diff := v - iii]
  if (nrow(wide) > quota) {
    set.seed(seed)
    wide <- wide[sample.int(nrow(wide), quota, replace = FALSE)]
  }
  wide[]
}

build_partition_sample_E <- function(files, max_rows, seed, diff = FALSE) {
  quota <- ceiling(max_rows / length(files) * 1.10)
  pieces <- vector("list", length(files))
  for (i in seq_along(files)) {
    message(if (diff) "Difference GAM" else "Global GAM",
            ": partition ", i, " of ", length(files))
    part <- readRDS(files[[i]])
    pieces[[i]] <- if (diff) {
      prepare_diff_partition_E(part, quota, seed + i)
    } else {
      prepare_global_partition_E(part, quota, seed + i)
    }
    rm(part); invisible(gc())
  }
  out <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  rm(pieces); invisible(gc())
  if (nrow(out) > max_rows) {
    set.seed(seed)
    out <- out[sample.int(nrow(out), max_rows, replace = FALSE)]
  }
  as.data.frame(out)
}

run_code_E <- function(data_name = "acc_long",
                       out_file = "results-code-E.txt",
                       fig_dir = "figures-code-E") {
  max_rows <- suppressWarnings(as.integer(Sys.getenv("CODE_E_MAX_ROWS", "500000")))
  if (is.na(max_rows) || max_rows < 10000L) max_rows <- 500000L
  max_diff_rows <- suppressWarnings(as.integer(Sys.getenv("CODE_E_MAX_DIFF_ROWS", "300000")))
  if (is.na(max_diff_rows) || max_diff_rows < 10000L) max_diff_rows <- 300000L
  seed <- suppressWarnings(as.integer(Sys.getenv("CODE_E_SEED", "20731")))
  if (is.na(seed)) seed <- 20260731L

  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  con <- file(out_file, open = "wt")
  sink_start <- sink.number(type = "output")
  sink(con, type = "output")
  on.exit({
    while (sink.number(type = "output") > sink_start) sink(type = "output")
    if (isOpen(con)) close(con)
  }, add = TRUE)

  tryCatch({
    cat("Code E: Generalized Additive Model Analysis\n")
    cat("Maximum global-GAM rows:", max_rows, "\n")
    cat("Maximum difference-GAM rows:", max_diff_rows, "\n")
    cat("Sampling seed:", seed, "\n\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' not found.")
    }
    acc_long_in <- get(data_name, envir = parent.frame(), inherits = TRUE)
    is_arrow <- inherits(acc_long_in,
      c("Dataset", "FileSystemDataset", "UnionDataset", "arrow_dplyr_query", "RecordBatchReader"))

    if (is_arrow) {
      files <- find_acc_long_parts_E()
      cat("Input type: disk-backed Arrow dataset\n")
      cat("Number of partitions:", length(files), "\n\n")
      acc_model <- build_partition_sample_E(files, max_rows, seed, diff = FALSE)
      diff_v_iii <- build_partition_sample_E(files, max_diff_rows, seed + 100000L, diff = TRUE)
    } else if (is.data.frame(acc_long_in) || is.data.table(acc_long_in)) {
      acc_model <- as.data.frame(prepare_global_partition_E(acc_long_in, max_rows, seed))
      diff_v_iii <- as.data.frame(prepare_diff_partition_E(acc_long_in, max_diff_rows, seed + 100000L))
    } else {
      stop("Unsupported acc_long class: ", paste(class(acc_long_in), collapse = ", "))
    }

    if (!nrow(acc_model)) stop("No usable rows remained for the global GAM.")
    if (!nrow(diff_v_iii)) stop("No complete Method v versus Method iii rows remained.")

    acc_model <- acc_model %>% mutate(
      method = droplevels(factor(method)),
      metric = droplevels(factor(metric)),
      design_block = droplevels(factor(design_block)),
      across(c(G, theta_true, theta_CV, sigma, mean_nsamp, p1, p2, accuracy), as.numeric)
    ) %>% filter(complete.cases(.), accuracy >= 0)

    diff_v_iii <- diff_v_iii %>% mutate(
      design_block = droplevels(factor(design_block)),
      metric = droplevels(factor(metric)),
      across(c(G, theta_true, theta_CV, sigma, mean_nsamp, p1, p2, diff), as.numeric)
    ) %>% filter(complete.cases(.), is.finite(diff))

    cat("Global GAM rows analyzed:", nrow(acc_model), "\n")
    cat("Difference GAM rows analyzed:", nrow(diff_v_iii), "\n\n")

    fit_gam <- bam(
      log1p(accuracy) ~ metric * method +
        s(p1, p2, by = method, k = 30) +
        design_block + G + theta_true + theta_CV + sigma + mean_nsamp,
      data = acc_model, method = "fREML", discrete = TRUE, nthreads = 1, gc.level = 1
    )
    cat("Global GAM Summary:\n"); print(summary(fit_gam))
    cat("\nGAM ANOVA Table:\n"); print(anova(fit_gam)); cat("\n")

    png(file.path(fig_dir, "code_E_global_gam_check.png"), 2400, 1800, res = 300)
    old_par <- par(no.readonly = TRUE); par(mfrow = c(2, 2)); gam.check(fit_gam); par(old_par); dev.off()

    methods_list <- levels(acc_model$method)
    png(file.path(fig_dir, "code_E_global_method_surfaces_panel.png"), 3600, 2400, res = 300)
    old_par <- par(no.readonly = TRUE); par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))
    for (m in methods_list) vis.gam(fit_gam, view = c("p1", "p2"), cond = list(method = m),
                                    plot.type = "contour", main = paste("Method:", m), color = "topo")
    par(old_par); dev.off()

    for (m in methods_list) {
      png(file.path(fig_dir, paste0("code_E_global_surface_method_", safe_label_E(m), ".png")),
          2400, 1800, res = 300)
      vis.gam(fit_gam, view = c("p1", "p2"), cond = list(method = m),
              plot.type = "contour", main = paste("Method:", m), color = "topo")
      dev.off()
    }

    fit_diff_gam <- bam(
      log1p(abs(diff)) ~ metric +
        s(p1, p2, by = metric, k = 30) +
        design_block + G + theta_true + theta_CV + sigma + mean_nsamp,
      data = diff_v_iii, method = "fREML", discrete = TRUE, nthreads = 1, gc.level = 1
    )
    cat("Difference GAM Summary (Method v - iii):\n"); print(summary(fit_diff_gam)); cat("\n")

    png(file.path(fig_dir, "code_E_diff_v_iii_gam_check.png"), 2400, 1800, res = 300)
    old_par <- par(no.readonly = TRUE); par(mfrow = c(2, 2)); gam.check(fit_diff_gam); par(old_par); dev.off()

    metrics_list <- levels(diff_v_iii$metric)
    png(file.path(fig_dir, "code_E_diff_v_iii_metric_surfaces_panel.png"), 3600, 2400, res = 300)
    old_par <- par(no.readonly = TRUE); par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    for (mm in metrics_list) vis.gam(fit_diff_gam, view = c("p1", "p2"), cond = list(metric = mm),
                                     plot.type = "contour", main = paste("Difference GAM metric:", mm), color = "topo")
    par(old_par); dev.off()

    for (mm in metrics_list) {
      png(file.path(fig_dir, paste0("code_E_diff_v_iii_surface_metric_", safe_label_E(mm), ".png")),
          2400, 1800, res = 300)
      vis.gam(fit_diff_gam, view = c("p1", "p2"), cond = list(metric = mm),
              plot.type = "contour", main = paste("Difference GAM metric:", mm), color = "topo")
      dev.off()
    }

    if (requireNamespace("gratia", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)) {
      p1 <- gratia::draw(fit_gam, data = acc_model, select = "s(p1,p2)", partial_match = TRUE)
      ggplot2::ggsave(file.path(fig_dir, "code_E_global_gratia_smooths.png"), p1,
                      width = 12, height = 8, dpi = 300, bg = "white")
      p2 <- gratia::appraise(fit_gam)
      ggplot2::ggsave(file.path(fig_dir, "code_E_global_gratia_appraise.png"), p2,
                      width = 12, height = 8, dpi = 300, bg = "white")
      p3 <- gratia::draw(fit_diff_gam, data = diff_v_iii, select = "s(p1,p2)", partial_match = TRUE)
      ggplot2::ggsave(file.path(fig_dir, "code_E_diff_v_iii_gratia_smooths.png"), p3,
                      width = 12, height = 8, dpi = 300, bg = "white")
      p4 <- gratia::appraise(fit_diff_gam)
      ggplot2::ggsave(file.path(fig_dir, "code_E_diff_v_iii_gratia_appraise.png"), p4,
                      width = 12, height = 8, dpi = 300, bg = "white")
    } else {
      cat("gratia/ggplot2 not installed; gratia plots were skipped.\n")
    }

    saveRDS(list(fit_gam = fit_gam, fit_diff_gam = fit_diff_gam,
                 global_rows = nrow(acc_model), difference_rows = nrow(diff_v_iii),
                 seed = seed, figure_directory = fig_dir),
            "code-E-model-results.rds", compress = "gzip")

    cat("\nAnalysis completed successfully.\n")
    invisible(list(fit_gam = fit_gam, fit_diff_gam = fit_diff_gam,
                   global_data = acc_model, difference_data = diff_v_iii,
                   figure_directory = fig_dir))
  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n", sep = "")
    stop(e)
  })
}

code_E_result <- run_code_E()
