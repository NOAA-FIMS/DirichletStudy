# code-E.R
# Revised 1-June-2026
# Revised to save all output graphs to PNG files
# Revised to handle memory errors in GAM fitting
# Revised 2-June-2026 to add signed Method v versus Method iii 
# log-ratio GAM with directional plots

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



    ############################################################################################

    # Directional signed log-ratio GAM: Method V relative to Method III
    # Positive fitted values mean Method V has larger error than Method III.
    # Negative fitted values mean Method V has smaller error than Method III.
    cat("\nAnalyzing signed directional log-ratio: Method 'v' relative to Method 'iii'.\n")

    ratio_values <- c(diff_v_iii$v, diff_v_iii$iii)
    min_positive_accuracy <- suppressWarnings(min(
      ratio_values[is.finite(ratio_values) & ratio_values > 0],
      na.rm = TRUE
    ))

    if (!is.finite(min_positive_accuracy)) {
      eps_ratio <- 1e-12
    } else {
      eps_ratio <- max(1e-12, 0.5 * min_positive_accuracy)
    }

    cat("Directional comparison epsilon:", format(eps_ratio, scientific = TRUE), "\n")
    cat("Response definition: y = log((v + eps) / (iii + eps)).\n")
    cat("Sign convention: y > 0 means Method v has larger error than Method iii; y < 0 means Method v has smaller error than Method iii.\n")

    signed_v_iii <- diff_v_iii %>%
      mutate(
        v_safe = pmax(v, 0),
        iii_safe = pmax(iii, 0),
        signed_log_ratio = log((v_safe + eps_ratio) / (iii_safe + eps_ratio))
      )

    signed_gam <- signed_v_iii %>%
      group_by(
        metric, design_block,
        G, theta_true, theta_CV, sigma, mean_nsamp,
        p1, p2
      ) %>%
      summarise(
        y = mean(signed_log_ratio, na.rm = TRUE),
        n = dplyr::n(),
        v_mean = mean(v, na.rm = TRUE),
        iii_mean = mean(iii, na.rm = TRUE),
        .groups = "drop"
      ) %>%
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
        y = as.numeric(y)
      )

    cat("Rows in signed_v_iii:", nrow(signed_v_iii), "\n")
    cat("Rows in aggregated signed_gam:", nrow(signed_gam), "\n")

    fit_signed_gam <- bam(
      y ~ metric +
        s(p1, p2, by = metric, k = 20) +
        design_block + G + theta_true + theta_CV + sigma + mean_nsamp,
      data = signed_gam,
      weights = n,
      method = "fREML",
      discrete = TRUE,
      nthreads = max(1, parallel::detectCores() - 1),
      gc.level = 2
    )

    cat("Signed log-ratio GAM Summary (Method v / Method iii):\n")
    print(summary(fit_signed_gam))
    cat("\nSigned log-ratio GAM ANOVA Table:\n")
    print(anova(fit_signed_gam))
    cat("\n")

    signed_check_file <- file.path(fig_dir, "code_E_signed_v_over_iii_gam_check.png")
    png(filename = signed_check_file, width = 2400, height = 1800, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(2, 2))
    gam.check(fit_signed_gam)
    par(old_par)
    dev.off()
    cat("Saved signed log-ratio GAM diagnostic plot:", signed_check_file, "\n")

    signed_plot_data <- signed_gam %>%
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

    p_signed_smooth <- gratia::draw(
      fit_signed_gam,
      data = signed_plot_data,
      select = "s(p1,p2)",
      partial_match = TRUE
    )

    signed_smooth_file <- file.path(fig_dir, "code_E_signed_v_over_iii_gratia_smooths.png")
    save_gratia_png(p_signed_smooth, signed_smooth_file, width = 12, height = 8)
    print(p_signed_smooth)
    cat("Saved signed log-ratio gratia smooth plot:", signed_smooth_file, "\n")

    p_signed_appraise <- gratia::appraise(fit_signed_gam)
    signed_appraise_file <- file.path(fig_dir, "code_E_signed_v_over_iii_gratia_appraise.png")
    save_gratia_png(p_signed_appraise, signed_appraise_file, width = 12, height = 8)
    print(p_signed_appraise)
    cat("Saved signed log-ratio gratia appraisal plot:", signed_appraise_file, "\n")

    # Publication-quality directional surfaces.
    # These predictions hold design factors at representative values and show
    # the signed Method v / Method iii error ratio across the feasible simplex.
    most_common_level <- function(x) {
      tab <- table(x, useNA = "no")
      names(tab)[which.max(tab)]
    }

    representative_design_block <- most_common_level(signed_gam$design_block)

    representative_values <- signed_gam %>%
      summarise(
        G = median(G, na.rm = TRUE),
        theta_true = median(theta_true, na.rm = TRUE),
        theta_CV = median(theta_CV, na.rm = TRUE),
        sigma = median(sigma, na.rm = TRUE),
        mean_nsamp = median(mean_nsamp, na.rm = TRUE),
        .groups = "drop"
      )

    p1_seq <- seq(min(signed_gam$p1, na.rm = TRUE), max(signed_gam$p1, na.rm = TRUE), length.out = 175)
    p2_seq <- seq(min(signed_gam$p2, na.rm = TRUE), max(signed_gam$p2, na.rm = TRUE), length.out = 175)

    make_signed_prediction_grid <- function(metric_name) {
      grid <- expand.grid(
        p1 = p1_seq,
        p2 = p2_seq
      ) %>%
        dplyr::filter(
          is.finite(p1),
          is.finite(p2),
          p1 >= 0,
          p2 >= 0,
          p1 + p2 <= 1
        ) %>%
        mutate(
          metric = factor(metric_name, levels = levels(signed_gam$metric)),
          design_block = factor(representative_design_block, levels = levels(signed_gam$design_block)),
          G = representative_values$G,
          theta_true = representative_values$theta_true,
          theta_CV = representative_values$theta_CV,
          sigma = representative_values$sigma,
          mean_nsamp = representative_values$mean_nsamp
        )

      pred <- predict(
        fit_signed_gam,
        newdata = grid,
        type = "response",
        se.fit = TRUE
      )

      grid %>%
        mutate(
          fit_log_ratio = as.numeric(pred$fit),
          se_log_ratio = as.numeric(pred$se.fit),
          ratio_v_over_iii = exp(fit_log_ratio),
          pct_error_difference = 100 * (ratio_v_over_iii - 1),
          positive_region = fit_log_ratio > 0,
          negative_region = fit_log_ratio < 0,
          interpretation = dplyr::case_when(
            fit_log_ratio > 0 ~ "Method v higher error than Method iii",
            fit_log_ratio < 0 ~ "Method v lower error than Method iii",
            TRUE ~ "No directional difference"
          )
        )
    }

    signed_pred_grid <- dplyr::bind_rows(
      lapply(metrics_list, make_signed_prediction_grid)
    )

    signed_region_summary <- signed_pred_grid %>%
      group_by(metric) %>%
      summarise(
        prediction_grid_cells = dplyr::n(),
        percent_positive_cells = 100 * mean(positive_region, na.rm = TRUE),
        percent_negative_cells = 100 * mean(negative_region, na.rm = TRUE),
        min_log_ratio = min(fit_log_ratio, na.rm = TRUE),
        max_log_ratio = max(fit_log_ratio, na.rm = TRUE),
        min_percent_error_difference = min(pct_error_difference, na.rm = TRUE),
        max_percent_error_difference = max(pct_error_difference, na.rm = TRUE),
        .groups = "drop"
      )

    cat("\nSigned directional prediction-grid summary.\n")
    print(signed_region_summary)

    signed_grid_file <- file.path(fig_dir, "code_E_signed_v_over_iii_prediction_grid.csv")
    utils::write.csv(signed_pred_grid, signed_grid_file, row.names = FALSE)
    cat("Saved signed directional prediction grid:", signed_grid_file, "\n")

    plot_caption <- paste0(
      "Predictions use ", representative_design_block,
      " design block and median numeric design covariates. ",
      "Positive values indicate Method v has larger error; negative values indicate Method v has smaller error."
    )

    p_signed_full <- ggplot2::ggplot(
      signed_pred_grid,
      ggplot2::aes(x = p1, y = p2, fill = pct_error_difference)
    ) +
      ggplot2::geom_raster() +
      ggplot2::geom_contour(
        ggplot2::aes(z = pct_error_difference),
        breaks = 0,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::facet_wrap(~ metric, ncol = 2) +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::scale_fill_gradient2(
        name = "% error difference\nv relative to iii",
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0
      ) +
      ggplot2::labs(
        title = "Signed directional comparison: Method v relative to Method iii",
        subtitle = "Black contour is equality: log((v + eps) / (iii + eps)) = 0",
        x = "p1",
        y = "p2",
        caption = plot_caption
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold"),
        legend.position = "right"
      )

    signed_full_file <- file.path(fig_dir, "code_E_signed_v_over_iii_surface_panel.png")
    save_gratia_png(p_signed_full, signed_full_file, width = 12, height = 8)
    print(p_signed_full)
    cat("Saved signed directional full surface panel:", signed_full_file, "\n")

    signed_positive_grid <- signed_pred_grid %>%
      mutate(
        positive_percent_higher_error = dplyr::if_else(
          positive_region,
          pct_error_difference,
          NA_real_
        )
      )

    p_signed_positive <- ggplot2::ggplot(
      signed_positive_grid,
      ggplot2::aes(x = p1, y = p2, fill = positive_percent_higher_error)
    ) +
      ggplot2::geom_raster(na.rm = FALSE) +
      ggplot2::geom_contour(
        data = signed_pred_grid,
        ggplot2::aes(x = p1, y = p2, z = pct_error_difference),
        breaks = 0,
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::facet_wrap(~ metric, ncol = 2) +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::scale_fill_gradient(
        name = "% higher error\nMethod v",
        low = "white",
        high = "#B2182B",
        na.value = "grey95"
      ) +
      ggplot2::labs(
        title = "Positive directional regions: Method v has larger error than Method iii",
        subtitle = "Colored cells show where log((v + eps) / (iii + eps)) > 0; black contour is equality.",
        x = "p1",
        y = "p2",
        caption = plot_caption
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold"),
        legend.position = "right"
      )

    signed_positive_file <- file.path(fig_dir, "code_E_signed_v_over_iii_positive_regions_panel.png")
    save_gratia_png(p_signed_positive, signed_positive_file, width = 12, height = 8)
    print(p_signed_positive)
    cat("Saved signed positive-region panel:", signed_positive_file, "\n")

    signed_negative_grid <- signed_pred_grid %>%
      mutate(
        negative_percent_lower_error = dplyr::if_else(
          negative_region,
          -pct_error_difference,
          NA_real_
        )
      )

    p_signed_negative <- ggplot2::ggplot(
      signed_negative_grid,
      ggplot2::aes(x = p1, y = p2, fill = negative_percent_lower_error)
    ) +
      ggplot2::geom_raster(na.rm = FALSE) +
      ggplot2::geom_contour(
        data = signed_pred_grid,
        ggplot2::aes(x = p1, y = p2, z = pct_error_difference),
        breaks = 0,
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::facet_wrap(~ metric, ncol = 2) +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::scale_fill_gradient(
        name = "% lower error\nMethod v",
        low = "white",
        high = "#2166AC",
        na.value = "grey95"
      ) +
      ggplot2::labs(
        title = "Negative directional regions: Method v has smaller error than Method iii",
        subtitle = "Colored cells show where log((v + eps) / (iii + eps)) < 0; black contour is equality.",
        x = "p1",
        y = "p2",
        caption = plot_caption
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold"),
        legend.position = "right"
      )

    signed_negative_file <- file.path(fig_dir, "code_E_signed_v_over_iii_negative_regions_panel.png")
    save_gratia_png(p_signed_negative, signed_negative_file, width = 12, height = 8)
    print(p_signed_negative)
    cat("Saved signed negative-region panel:", signed_negative_file, "\n")

    for (mm in metrics_list) {
      metric_safe <- safe_label(mm)
      metric_all_data <- dplyr::filter(signed_pred_grid, metric == mm)
      metric_positive_data <- dplyr::filter(signed_positive_grid, metric == mm)
      metric_negative_data <- dplyr::filter(signed_negative_grid, metric == mm)

      p_metric_positive <- ggplot2::ggplot(
        metric_positive_data,
        ggplot2::aes(x = p1, y = p2, fill = positive_percent_higher_error)
      ) +
        ggplot2::geom_raster(na.rm = FALSE) +
        ggplot2::geom_contour(
          data = metric_all_data,
          ggplot2::aes(x = p1, y = p2, z = pct_error_difference),
          breaks = 0,
          inherit.aes = FALSE,
          color = "black",
          linewidth = 0.35
        ) +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::scale_fill_gradient(
          name = "% higher error\nMethod v",
          low = "white",
          high = "#B2182B",
          na.value = "grey95"
        ) +
        ggplot2::labs(
          title = paste("Positive directional region:", mm),
          subtitle = "Method v has larger error than Method iii; black contour is equality.",
          x = "p1",
          y = "p2",
          caption = plot_caption
        ) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          legend.position = "right"
        )

      metric_positive_file <- file.path(
        fig_dir,
        paste0("code_E_signed_v_over_iii_positive_region_metric_", metric_safe, ".png")
      )
      save_gratia_png(p_metric_positive, metric_positive_file, width = 8, height = 6)
      print(p_metric_positive)
      cat("Saved signed positive-region metric plot:", metric_positive_file, "\n")

      p_metric_negative <- ggplot2::ggplot(
        metric_negative_data,
        ggplot2::aes(x = p1, y = p2, fill = negative_percent_lower_error)
      ) +
        ggplot2::geom_raster(na.rm = FALSE) +
        ggplot2::geom_contour(
          data = metric_all_data,
          ggplot2::aes(x = p1, y = p2, z = pct_error_difference),
          breaks = 0,
          inherit.aes = FALSE,
          color = "black",
          linewidth = 0.35
        ) +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::scale_fill_gradient(
          name = "% lower error\nMethod v",
          low = "white",
          high = "#2166AC",
          na.value = "grey95"
        ) +
        ggplot2::labs(
          title = paste("Negative directional region:", mm),
          subtitle = "Method v has smaller error than Method iii; black contour is equality.",
          x = "p1",
          y = "p2",
          caption = plot_caption
        ) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          legend.position = "right"
        )

      metric_negative_file <- file.path(
        fig_dir,
        paste0("code_E_signed_v_over_iii_negative_region_metric_", metric_safe, ".png")
      )
      save_gratia_png(p_metric_negative, metric_negative_file, width = 8, height = 6)
      print(p_metric_negative)
      cat("Saved signed negative-region metric plot:", metric_negative_file, "\n")
    }
    cat("\nAnalysis completed successfully.\n")

    invisible(list(
      fit_gam = fit_gam,
      fit_diff_gam = fit_diff_gam,
      fit_signed_gam = fit_signed_gam,
      signed_prediction_grid = signed_pred_grid,
      figure_directory = fig_dir
    ))

  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n")
    stop(e)
  })
}

run_code_E()
