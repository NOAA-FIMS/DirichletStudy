# code-C1.R
# Revised to match the output style of code-A.R

library(dplyr)
library(tidyr)

run_code_C1 <- function(data_name = "acc_long", out_file = "results-code-C1.txt") {
  
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
    cat("Code C1: Pairwise Hotelling's T^2 (MANOVA) results\n")
    cat("Output file:", out_file, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' not found.")
    }
    acc_long <- get(data_name, envir = parent.frame(), inherits = TRUE)

    # Helper function for data transformation

make_pairwise_diff <- function(dat, m1, m2, eps = 1e-6) {

  id_vars <- c(
    "example_id", "design_block", "G",
    "theta_true", "theta_CV", "sigma", "mean_nsamp",
    "p1", "p2"
  )

  dat %>%
    select(
      all_of(id_vars),
      method, metric, accuracy
    ) %>%
    filter(method %in% c(m1, m2)) %>%
    mutate(
      accuracy = as.numeric(accuracy),

      # Standardize metric names so MANOVA can find them later.
      metric = case_when(
        tolower(as.character(metric)) %in% c("rmse", "rmse_norm") ~ "rmse",
        tolower(as.character(metric)) %in% c("l1", "l1_norm", "l1 norm") ~ "L1",
        tolower(as.character(metric)) %in% c("linf", "l_inf", "linf_norm", "l_inf_norm", "linf norm") ~ "Linf",
        TRUE ~ as.character(metric)
      ),

      accuracy_log = log(accuracy + eps)
    ) %>%

    # Collapse duplicate cells before widening.
    summarise(
      accuracy_log = mean(accuracy_log, na.rm = TRUE),
      .by = c(all_of(id_vars), metric, method)
    ) %>%

    pivot_wider(
      id_cols = c(all_of(id_vars), metric),
      names_from = method,
      values_from = accuracy_log
    ) %>%

    filter(
      !is.na(.data[[m1]]),
      !is.na(.data[[m2]])
    ) %>%

    mutate(
      diff = .data[[m1]] - .data[[m2]]
    ) %>%

    select(
      all_of(id_vars),
      metric,
      diff
    ) %>%

    pivot_wider(
      id_cols = all_of(id_vars),
      names_from = metric,
      values_from = diff
    )
}

    # Define method pairs for comparison
    methods <- c("i", "ii", "iii", "iv", "v")
    pairs <- combn(methods, 2, simplify = FALSE)

    for (p in pairs) {
      m_target <- p[2]
      m_ref <- p[1]
      comp_name <- paste0(m_target, "_vs_", m_ref)
      
      cat("==================================================================\n")
      cat("Comparison:", comp_name, "\n")
      cat("==================================================================\n")
      
      # Prepare multivariate difference data
      diff_mv <- make_pairwise_diff(acc_long, m_target, m_ref)
      
      # Hotelling's T^2 equivalent (Intercept only MANOVA)
      fit_hotelling <- stats::manova(cbind(rmse, L1, Linf) ~ 1, data = diff_mv)
      
      cat("MANOVA results: Pillai test\n")
      print(summary(fit_hotelling, test = "Pillai"))
      cat("\n")
    }

    cat("Analysis completed successfully.\n")
  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n")
    stop(e)
  })
}

run_code_C1()