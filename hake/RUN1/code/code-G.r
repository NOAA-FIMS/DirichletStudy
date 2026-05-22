# code-G.R
# Revised 20-May-2026
# Revised again to replace the unexported namespace call for pairwise contrasts
# with emmeans::contrast(..., method = "pairwise").
# Purpose:
#   Rank-based accuracy analysis for compositional likelihood simulations.
#   This revision avoids the emmeans reference-grid error and removes the
#   singular mixed-model random-effect structure in the original script.
#
# Key changes from the original code-G.R:
#   1. Converts ordered design factors that are used as covariates to numeric
#      before modeling. This prevents emmeans from building a large factorial
#      reference grid over G, theta_true, theta_CV, sigma, and mean_nsamp.
#   2. Removes random intercepts for example_id and dataset_id. Because ranks
#      are computed within dataset_id and metric, the average rank is fixed
#      within each ranking block, which makes random intercepts redundant and
#      can produce boundary/singular fits.
#   3. Adds a metric-adjusted composite rank analysis that averages ranks over
#      RMSE, L1, and Linf. This reduces pseudo-replication from highly collinear
#      accuracy metrics.
#   4. Uses direct paired rank contrasts by dataset_id and method. Lower ranks
#      indicate better accuracy.
#   5. Saves plots to PNG files for reproducibility.

library(dplyr)
library(tidyr)
library(emmeans)
library(ggplot2)

run_code_G <- function(
  data_name = "acc_long",
  out_file = "results-code-G-revised.txt",
  fig_dir = "figures-code-G"
) {

  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE)
  }

  safe_numeric <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    suppressWarnings(as.numeric(x))
  }

  finite_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
  }

  save_plot <- function(plot_obj, filename, width = 10, height = 7, dpi = 300) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot_obj,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }

  make_paired_contrasts <- function(dat, response, block_vars, group_var = NULL) {
    # Direct paired rank contrasts. For a contrast A - B, a negative estimate
    # means method A has a lower, and therefore better, rank than method B.

    if (!response %in% names(dat)) {
      stop("Response column '", response, "' not found.")
    }

    method_levels <- levels(dat$method)
    if (is.null(method_levels)) {
      method_levels <- sort(unique(as.character(dat$method)))
    }

    one_group <- function(dd, group_value = NULL) {
      select_vars <- unique(c(block_vars, "method", response))

      wide <- dd %>%
        dplyr::select(dplyr::all_of(select_vars)) %>%
        dplyr::distinct() %>%
        tidyr::pivot_wider(
          names_from = method,
          values_from = dplyr::all_of(response)
        )

      methods_here <- method_levels[method_levels %in% names(wide)]
      if (length(methods_here) < 2L) {
        return(data.frame())
      }

      out <- lapply(utils::combn(methods_here, 2L, simplify = FALSE), function(mm) {
        diff_vec <- wide[[mm[1L]]] - wide[[mm[2L]]]
        diff_vec <- diff_vec[is.finite(diff_vec)]
        n_used <- length(diff_vec)

        if (n_used == 0L) {
          mean_diff <- NA_real_
          se_diff <- NA_real_
          t_value <- NA_real_
          p_value <- NA_real_
        } else {
          mean_diff <- mean(diff_vec)
          se_diff <- stats::sd(diff_vec) / sqrt(n_used)
          if (!is.finite(se_diff) || se_diff == 0) {
            t_value <- NA_real_
            p_value <- NA_real_
          } else {
            t_value <- mean_diff / se_diff
            p_value <- 2 * stats::pt(abs(t_value), df = n_used - 1L, lower.tail = FALSE)
          }
        }

        data.frame(
          contrast = paste(mm[1L], "-", mm[2L]),
          estimate = mean_diff,
          SE = se_diff,
          df = ifelse(n_used > 1L, n_used - 1L, NA_integer_),
          t_value = t_value,
          p_value = p_value,
          n_pairs = n_used,
          interpretation = dplyr::case_when(
            is.na(mean_diff) ~ "not estimable",
            mean_diff < 0 ~ paste(mm[1L], "has lower mean rank than", mm[2L]),
            mean_diff > 0 ~ paste(mm[2L], "has lower mean rank than", mm[1L]),
            TRUE ~ "equal mean rank"
          ),
          stringsAsFactors = FALSE
        )
      })

      ans <- dplyr::bind_rows(out)
      if (nrow(ans) > 0L) {
        ans$p_value_holm <- stats::p.adjust(ans$p_value, method = "holm")
        if (!is.null(group_var)) {
          ans[[group_var]] <- group_value
          ans <- ans %>% dplyr::relocate(dplyr::all_of(group_var))
        }
      }
      ans
    }

    if (is.null(group_var)) {
      one_group(dat)
    } else {
      split_dat <- split(dat, dat[[group_var]], drop = TRUE)
      dplyr::bind_rows(lapply(names(split_dat), function(g) one_group(split_dat[[g]], g)))
    }
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
    cat("Code G: Revised Rank-Based Accuracy Analysis\n")
    cat("Output file:", out_file, "\n")
    cat("Figure output directory:", fig_dir, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' not found.")
    }
    acc_long <- get(data_name, envir = parent.frame(), inherits = TRUE)

    required_cols <- c(
      "dataset_id", "example_id", "method", "metric", "accuracy",
      "design_block", "G", "theta_true", "theta_CV", "sigma",
      "mean_nsamp", "p1", "p2"
    )
    missing_cols <- setdiff(required_cols, names(acc_long))
    if (length(missing_cols) > 0L) {
      stop("Input data are missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    method_order <- c("i", "ii", "iii", "iv", "v")
    observed_methods <- unique(as.character(acc_long$method))
    method_levels <- c(
      method_order[method_order %in% observed_methods],
      sort(setdiff(observed_methods, method_order))
    )

    metric_order <- c("rmse", "L1_norm", "Linf_norm", "L1", "Linf")
    observed_metrics <- unique(as.character(acc_long$metric))
    metric_levels <- c(
      metric_order[metric_order %in% observed_metrics],
      sort(setdiff(observed_metrics, metric_order))
    )

    acc_long <- acc_long %>%
      dplyr::mutate(
        dataset_id = factor(dataset_id),
        example_id = factor(example_id),
        method = factor(as.character(method), levels = method_levels),
        metric = factor(as.character(metric), levels = metric_levels),
        design_block = factor(design_block),
        G = safe_numeric(G),
        theta_true = safe_numeric(theta_true),
        theta_CV = safe_numeric(theta_CV),
        sigma = safe_numeric(sigma),
        mean_nsamp = safe_numeric(mean_nsamp),
        p1 = safe_numeric(p1),
        p2 = safe_numeric(p2),
        accuracy = as.numeric(accuracy)
      ) %>%
      dplyr::filter(
        !is.na(dataset_id), !is.na(method), !is.na(metric),
        is.finite(accuracy)
      )

    cat("Computing ranks within dataset_id and metric.\n")
    rank_dat <- acc_long %>%
      dplyr::group_by(dataset_id, metric) %>%
      dplyr::mutate(
        rank = rank(accuracy, ties.method = "average"),
        best_accuracy = min(accuracy, na.rm = TRUE),
        is_best = accuracy == best_accuracy
      ) %>%
      dplyr::ungroup()

    cat("\nRank interpretation:\n")
    cat("  Rank 1 is the lowest accuracy error and therefore the best method.\n")
    cat("  Pairwise rank contrasts use A - B; negative estimates favor A.\n\n")

    cat("Descriptive Statistics: Mean Rank and Probability of Being Best\n")
    cat("----------------------------------------------------------------\n")
    rank_summary <- rank_dat %>%
      dplyr::group_by(metric, method) %>%
      dplyr::summarise(
        mean_rank = mean(rank, na.rm = TRUE),
        prob_best = mean(is_best, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      )
    print(as.data.frame(rank_summary))
    cat("\n")

    cat("Checking rank concordance among metrics.\n")
    metric_rank_wide <- rank_dat %>%
      dplyr::select(dataset_id, method, metric, rank) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = metric, values_from = rank)

    metric_cols <- setdiff(names(metric_rank_wide), c("dataset_id", "method"))
    if (length(metric_cols) >= 2L) {
      metric_rank_cor <- stats::cor(
        metric_rank_wide[, metric_cols, drop = FALSE],
        use = "pairwise.complete.obs"
      )
      print(round(metric_rank_cor, 4))
      cat("\n")
    }

    cat("Metric-specific paired rank contrasts.\n")
    cat("--------------------------------------\n")
    paired_by_metric <- make_paired_contrasts(
      dat = rank_dat,
      response = "rank",
      block_vars = c("dataset_id"),
      group_var = "metric"
    )
    print(as.data.frame(paired_by_metric))
    cat("\n")

    cat("Creating metric-adjusted composite rank by averaging ranks over metrics.\n")
    composite_dat <- rank_dat %>%
      dplyr::group_by(
        dataset_id, example_id, design_block, G, theta_true, theta_CV,
        sigma, mean_nsamp, p1, p2, method
      ) %>%
      dplyr::summarise(
        mean_rank = mean(rank, na.rm = TRUE),
        prob_best_across_metrics = mean(is_best, na.rm = TRUE),
        n_metrics = dplyr::n_distinct(metric),
        .groups = "drop"
      )

    cat("\nComposite Rank Summary by Method\n")
    cat("--------------------------------\n")
    composite_summary <- composite_dat %>%
      dplyr::group_by(method) %>%
      dplyr::summarise(
        mean_composite_rank = mean(mean_rank, na.rm = TRUE),
        sd_composite_rank = stats::sd(mean_rank, na.rm = TRUE),
        mean_probability_best = mean(prob_best_across_metrics, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      )
    print(as.data.frame(composite_summary))
    cat("\n")

    cat("Composite paired rank contrasts.\n")
    cat("--------------------------------\n")
    composite_contrasts <- make_paired_contrasts(
      dat = composite_dat,
      response = "mean_rank",
      block_vars = c("dataset_id"),
      group_var = NULL
    )
    print(as.data.frame(composite_contrasts))
    cat("\n")

    cat("Fitting metric-adjusted fixed-effect rank model.\n")
    cat("Random intercepts are intentionally not used because ranks are fixed within ranking blocks.\n")
    fit_comp <- stats::lm(
      mean_rank ~ method * (design_block + G + theta_true + theta_CV + sigma + mean_nsamp + p1 + p2),
      data = composite_dat
    )

    cat("\nComposite Rank Model ANOVA Table\n")
    cat("--------------------------------\n")
    print(stats::anova(fit_comp))
    cat("\n")

    at_comp <- list(
      G = finite_mean(composite_dat$G),
      theta_true = finite_mean(composite_dat$theta_true),
      theta_CV = finite_mean(composite_dat$theta_CV),
      sigma = finite_mean(composite_dat$sigma),
      mean_nsamp = finite_mean(composite_dat$mean_nsamp),
      p1 = finite_mean(composite_dat$p1),
      p2 = finite_mean(composite_dat$p2)
    )

    cat("Estimated marginal mean composite ranks by method.\n")
    emm_comp <- emmeans::emmeans(
      fit_comp,
      specs = ~ method,
      at = at_comp,
      weights = "proportional",
      rg.limit = 50000
    )
    print(emm_comp)
    cat("\nPairwise EMM contrasts for composite ranks.\n")
    print(emmeans::contrast(emm_comp, method = "pairwise", adjust = "holm"))
    cat("\n")

    cat("Metric-specific fixed-effect rank models and EMM contrasts.\n")
    cat("----------------------------------------------------------\n")
    metric_model_results <- list()
    for (mm in levels(rank_dat$metric)) {
      dd <- rank_dat %>% dplyr::filter(metric == mm)
      if (nrow(dd) == 0L) next

      cat("\nMetric:", mm, "\n")
      fit_metric <- stats::lm(
        rank ~ method * (design_block + G + theta_true + theta_CV + sigma + mean_nsamp + p1 + p2),
        data = dd
      )
      print(stats::anova(fit_metric))

      at_metric <- list(
        G = finite_mean(dd$G),
        theta_true = finite_mean(dd$theta_true),
        theta_CV = finite_mean(dd$theta_CV),
        sigma = finite_mean(dd$sigma),
        mean_nsamp = finite_mean(dd$mean_nsamp),
        p1 = finite_mean(dd$p1),
        p2 = finite_mean(dd$p2)
      )

      emm_metric <- emmeans::emmeans(
        fit_metric,
        specs = ~ method,
        at = at_metric,
        weights = "proportional",
        rg.limit = 50000
      )
      print(emm_metric)
      print(emmeans::contrast(emm_metric, method = "pairwise", adjust = "holm"))

      metric_model_results[[as.character(mm)]] <- list(
        fit = fit_metric,
        emmeans = emm_metric
      )
    }

    cat("\nGenerating visualization plots.\n")

    p_best <- ggplot(rank_summary, aes(x = method, y = prob_best, fill = method)) +
      geom_col(color = "black", alpha = 0.8) +
      facet_wrap(~ metric) +
      theme_minimal() +
      labs(
        title = "Probability of Method Being Best",
        subtitle = "Rank 1 is assigned to the lowest accuracy error within each dataset and metric",
        y = "Probability of best rank",
        x = "Estimation method"
      ) +
      guides(fill = "none")
    best_file <- file.path(fig_dir, "code_G_probability_best_by_metric.png")
    save_plot(p_best, best_file)
    cat("Saved plot:", best_file, "\n")

    p_comp <- ggplot(composite_summary, aes(x = method, y = mean_composite_rank)) +
      geom_point(size = 3) +
      geom_errorbar(
        aes(
          ymin = mean_composite_rank - 1.96 * sd_composite_rank / sqrt(n),
          ymax = mean_composite_rank + 1.96 * sd_composite_rank / sqrt(n)
        ),
        width = 0.15
      ) +
      theme_minimal() +
      labs(
        title = "Mean Composite Rank by Method",
        subtitle = "Ranks are averaged over RMSE, L1, and Linf; lower values are better",
        y = "Mean composite rank",
        x = "Estimation method"
      )
    comp_file <- file.path(fig_dir, "code_G_composite_mean_rank.png")
    save_plot(p_comp, comp_file)
    cat("Saved plot:", comp_file, "\n")

    rank_dist_dat <- rank_dat %>%
      dplyr::group_by(method, metric, rank) %>%
      dplyr::tally(name = "n") %>%
      dplyr::group_by(method, metric) %>%
      dplyr::mutate(percentage = n / sum(n)) %>%
      dplyr::ungroup()

    p_dist <- ggplot(rank_dist_dat, aes(x = method, y = percentage, fill = factor(rank))) +
      geom_col(color = "white", linewidth = 0.2) +
      facet_wrap(~ metric) +
      theme_minimal() +
      labs(
        title = "Complete Rank Distribution Profile",
        subtitle = "Lower ranks indicate lower error",
        y = "Proportion of datasets",
        x = "Estimation method",
        fill = "Rank"
      )
    dist_file <- file.path(fig_dir, "code_G_rank_distribution_by_metric.png")
    save_plot(p_dist, dist_file)
    cat("Saved plot:", dist_file, "\n")

    diag_dat <- data.frame(
      fitted = stats::fitted(fit_comp),
      resid = stats::resid(fit_comp)
    )
    p_res <- ggplot(diag_dat, aes(x = fitted, y = resid)) +
      geom_point(alpha = 0.15) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_minimal() +
      labs(
        title = "Composite Rank Model Diagnostics",
        x = "Fitted composite rank",
        y = "Residual"
      )
    res_file <- file.path(fig_dir, "code_G_composite_rank_residuals.png")
    save_plot(p_res, res_file)
    cat("Saved plot:", res_file, "\n")

    cat("\nAnalysis completed successfully.\n")

    invisible(list(
      rank_summary = rank_summary,
      paired_by_metric = paired_by_metric,
      composite_summary = composite_summary,
      composite_contrasts = composite_contrasts,
      fit_comp = fit_comp,
      emm_comp = emm_comp,
      metric_model_results = metric_model_results,
      figure_directory = fig_dir
    ))

  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n")
    stop(e)
  })
}

run_code_G()
