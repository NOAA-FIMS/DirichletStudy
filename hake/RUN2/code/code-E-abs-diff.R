# code-E.R
# Revised 1-June-2026
# Revised to save all output graphs to PNG files.
# Revised to handle memory errors in GAM fitting.

library(mgcv)
library(dplyr)

run_code_E <- function(
  data_name = "acc_long",
  out_file = "results-code-E.txt",
  fig_dir = "figures-code-E"
) {

  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE)
  }

  safe_label <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    tolower(x)
  }

  save_gratia_png <- function(plot_obj, filename, width = 12, height = 8, dpi = 300) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot_obj,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }

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
    cat("Code E: Generalized Additive Model (GAM) Analysis\n")
    cat("Output file:", out_file, "\n")
    cat("Figure output directory:", fig_dir, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' not found.")
    }

    acc_long <- get(data_name, envir = parent.frame(), inherits = TRUE)

    acc_long <- acc_long %>%
      mutate(
        method = factor(method),
        metric = factor(metric),
        design_block = factor(design_block),
        G = as.numeric(G),
        theta_true = as.numeric(theta_true),
        theta_CV = as.numeric(theta_CV),
        sigma = as.numeric(sigma),
        mean_nsamp = as.numeric(mean_nsamp),
        p1 = as.numeric(p1),
        p2 = as.numeric(p2),
        accuracy = as.numeric(accuracy)
      )
    
    # Reduce memory use by aggregating repeated design/grid cells.
    # The GAM estimates the conditional mean of log1p(accuracy), so repeated
    # rows with the same covariates can be represented by their mean response
    # and a frequency weight.
    
    acc_gam <- acc_long %>%
      mutate(
        y = log1p(accuracy)
      ) %>%
      group_by(
        metric, method, design_block,
        G, theta_true, theta_CV, sigma, mean_nsamp,
        p1, p2
      ) %>%
      summarise(
        y = mean(y, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      )
    
    cat("Rows in original acc_long:", nrow(acc_long), "\n")
    cat("Rows in aggregated acc_gam:", nrow(acc_gam), "\n")
    
    # Global GAM: Spatial error surface by method
    cat("Fitting Global GAM with spatial splines s(p1, p2) by method.\n")
    fit_gam <- bam(
      y ~ metric * method +
        s(p1, p2, by = method, k = 20) +
        design_block + G + theta_true + theta_CV + sigma + mean_nsamp,
      data = acc_gam,
      weights = n,
      method = "fREML",
      discrete = TRUE,
      nthreads = max(1, parallel::detectCores() - 1),
      gc.level = 2
    )

    cat("Global GAM Summary:\n")
    print(summary(fit_gam))
    cat("\nGAM ANOVA Table:\n")
    print(anova(fit_gam))
    cat("\n")

    # Technical diagnostics
    # Save the gam.check diagnostic panel to a PNG file.
    graphics.off()
    gam_check_file <- file.path(fig_dir, "code_E_global_gam_check.png")
    png(filename = gam_check_file, width = 2400, height = 1800, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(2, 2))
    gam.check(fit_gam)
    par(old_par)
    dev.off()
    cat("Saved GAM diagnostic plot:", gam_check_file, "\n")

    # Visualizing Method-Specific Surfaces using base mgcv
    methods_list <- levels(acc_long$method)

    # Save one combined panel of all method-specific p1/p2 surfaces.
    surface_panel_file <- file.path(fig_dir, "code_E_global_method_surfaces_panel.png")
    png(filename = surface_panel_file, width = 3600, height = 2400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))
    for (m in methods_list) {
      vis.gam(
        fit_gam,
        view = c("p1", "p2"),
        cond = list(method = m),
        plot.type = "contour",
        main = paste("Method:", m),
        color = "topo"
      )
    }
    par(old_par)
    dev.off()
    cat("Saved method surface panel:", surface_panel_file, "\n")

    # Save separate p1/p2 surface PNG files by method.
    for (m in methods_list) {
      method_file <- file.path(
        fig_dir,
        paste0("code_E_global_surface_method_", safe_label(m), ".png")
      )

      png(filename = method_file, width = 2400, height = 1800, res = 300)
      vis.gam(
        fit_gam,
        view = c("p1", "p2"),
        cond = list(method = m),
        plot.type = "contour",
        main = paste("Method:", m),
        color = "topo"
      )
      dev.off()
      cat("Saved method surface plot:", method_file, "\n")
    }

    # Modern Interpretation via gratia
    if (!requireNamespace("gratia", quietly = TRUE)) {
      install.packages("gratia")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      install.packages("ggplot2")
    }
    library(gratia)
    library(ggplot2)

    # Check the exact smooth labels generated by mgcv
    smooth_names <- gratia::smooths(fit_gam)

    cat("\nSmooth terms in fit_gam:\n")
    print(smooth_names)

    # Use explicit data with the same variable classes used to fit the GAM.
    plot_data <- acc_gam %>%
      mutate(
        method = factor(method, levels = levels(acc_long$method)),
        metric = factor(metric, levels = levels(acc_long$metric)),
        design_block = factor(design_block, levels = levels(acc_long$design_block)),
        G = as.numeric(G),
        theta_true = as.numeric(theta_true),
        theta_CV = as.numeric(theta_CV),
        sigma = as.numeric(sigma),
        mean_nsamp = as.numeric(mean_nsamp),
        p1 = as.numeric(p1),
        p2 = as.numeric(p2)
      )

    cat("\nDrawing method-specific smooth terms with gratia.\n")

    p_smooth <- gratia::draw(
      fit_gam,
      data = plot_data,
      select = "s(p1,p2)",
      partial_match = TRUE
    )

    smooth_file <- file.path(fig_dir, "code_E_global_gratia_smooths.png")
    save_gratia_png(p_smooth, smooth_file, width = 12, height = 8)
    print(p_smooth)
    cat("Saved gratia smooth plot:", smooth_file, "\n")

    # Comprehensive model appraisal
    cat("\nRunning gratia model appraisal.\n")
    p_appraise <- gratia::appraise(fit_gam)
    appraise_file <- file.path(fig_dir, "code_E_global_gratia_appraise.png")
    save_gratia_png(p_appraise, appraise_file, width = 12, height = 8)
    print(p_appraise)
    cat("Saved gratia appraisal plot:", appraise_file, "\n")

    # Partial effects of categorical parametric terms only
    # Avoid drawing all numeric parametric terms automatically because this
    # is where type mismatches can occur in the generated prediction grid.
    cat("\nDrawing categorical parametric effects with gratia.\n")

    gg_tensor <- gratia::draw(
      fit_gam,
      data = plot_data,
      parametric = TRUE,
      terms = c("metric", "method", "design_block")
    )

    tensor_file <- file.path(fig_dir, "code_E_global_parametric_terms.png")
    save_gratia_png(gg_tensor, tensor_file, width = 12, height = 8)
    print(gg_tensor)
    cat("Saved parametric-term plot:", tensor_file, "\n")

    ############################################################################################

    # Difference GAM: Comparing Method V vs Method III
    cat("\nAnalyzing difference between Method 'v' and Method 'iii'.\n")

    diff_v_iii <- acc_long %>%
      filter(method %in% c("v", "iii")) %>%
      group_by(
        example_id, design_block, G, theta_true, theta_CV, sigma,
        mean_nsamp, p1, p2, metric, method
      ) %>%
      summarise(
        accuracy = mean(accuracy, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = method,
        values_from = accuracy
      ) %>%
      filter(
        !is.na(v),
        !is.na(iii)
      ) %>%
      mutate(
        diff = v - iii
      )

    diff_v_iii <- diff_v_iii %>%
      mutate(
        design_block = factor(design_block, levels = levels(acc_long$design_block)),
        metric = factor(metric, levels = levels(acc_long$metric)),
        G = as.numeric(G),
        theta_true = as.numeric(theta_true),
        theta_CV = as.numeric(theta_CV),
        sigma = as.numeric(sigma),
        mean_nsamp = as.numeric(mean_nsamp),
        p1 = as.numeric(p1),
        p2 = as.numeric(p2),
        diff = as.numeric(diff)
      )

    diff_gam <- diff_v_iii %>%
      mutate(
        y = log1p(abs(diff))
      ) %>%
      group_by(
        metric, design_block,
        G, theta_true, theta_CV, sigma, mean_nsamp,
        p1, p2
      ) %>%
      summarise(
        y = mean(y, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      )
    
    cat("Rows in original diff_v_iii:", nrow(diff_v_iii), "\n")
    cat("Rows in aggregated diff_gam:", nrow(diff_gam), "\n")
    
    fit_diff_gam <- bam(
      y ~ metric +
        s(p1, p2, by = metric, k = 20) +
        design_block + G + theta_true + theta_CV + sigma + mean_nsamp,
      data = diff_gam,
      weights = n,
      method = "fREML",
      discrete = TRUE,
      nthreads = max(1, parallel::detectCores() - 1),
      gc.level = 2
    )

    cat("Difference GAM Summary (Method v - iii):\n")
    print(summary(fit_diff_gam))
    cat("\n")

    # Save diagnostics for the difference GAM.
    diff_check_file <- file.path(fig_dir, "code_E_diff_v_iii_gam_check.png")
    png(filename = diff_check_file, width = 2400, height = 1800, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(2, 2))
    gam.check(fit_diff_gam)
    par(old_par)
    dev.off()
    cat("Saved difference GAM diagnostic plot:", diff_check_file, "\n")

    # Save difference-GAM p1/p2 surfaces by metric.
    metrics_list <- levels(diff_v_iii$metric)

    diff_surface_panel_file <- file.path(fig_dir, "code_E_diff_v_iii_metric_surfaces_panel.png")
    png(filename = diff_surface_panel_file, width = 3600, height = 2400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    for (mm in metrics_list) {
      vis.gam(
        fit_diff_gam,
        view = c("p1", "p2"),
        cond = list(metric = mm),
        plot.type = "contour",
        main = paste("Difference GAM metric:", mm),
        color = "topo"
      )
    }
    par(old_par)
    dev.off()
    cat("Saved difference-GAM surface panel:", diff_surface_panel_file, "\n")

    for (mm in metrics_list) {
      metric_file <- file.path(
        fig_dir,
        paste0("code_E_diff_v_iii_surface_metric_", safe_label(mm), ".png")
      )

      png(filename = metric_file, width = 2400, height = 1800, res = 300)
      vis.gam(
        fit_diff_gam,
        view = c("p1", "p2"),
        cond = list(metric = mm),
        plot.type = "contour",
        main = paste("Difference GAM metric:", mm),
        color = "topo"
      )
      dev.off()
      cat("Saved difference-GAM metric surface plot:", metric_file, "\n")
    }

    diff_plot_data <- diff_gam %>%
      mutate(
        design_block = factor(design_block, levels = levels(acc_long$design_block)),
        metric = factor(metric, levels = levels(acc_long$metric)),
        G = as.numeric(G),
        theta_true = as.numeric(theta_true),
        theta_CV = as.numeric(theta_CV),
        sigma = as.numeric(sigma),
        mean_nsamp = as.numeric(mean_nsamp),
        p1 = as.numeric(p1),
        p2 = as.numeric(p2)
      )
    
    p_diff_smooth <- gratia::draw(
      fit_diff_gam,
      data = diff_plot_data,
      select = "s(p1,p2)",
      partial_match = TRUE
    )

    diff_smooth_file <- file.path(fig_dir, "code_E_diff_v_iii_gratia_smooths.png")
    save_gratia_png(p_diff_smooth, diff_smooth_file, width = 12, height = 8)
    print(p_diff_smooth)
    cat("Saved difference-GAM gratia smooth plot:", diff_smooth_file, "\n")

    p_diff_appraise <- gratia::appraise(fit_diff_gam)
    diff_appraise_file <- file.path(fig_dir, "code_E_diff_v_iii_gratia_appraise.png")
    save_gratia_png(p_diff_appraise, diff_appraise_file, width = 12, height = 8)
    print(p_diff_appraise)
    cat("Saved difference-GAM gratia appraisal plot:", diff_appraise_file, "\n")

    cat("\nAnalysis completed successfully.\n")

    invisible(list(
      fit_gam = fit_gam,
      fit_diff_gam = fit_diff_gam,
      figure_directory = fig_dir
    ))

  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n")
    stop(e)
  })
}

run_code_E()
