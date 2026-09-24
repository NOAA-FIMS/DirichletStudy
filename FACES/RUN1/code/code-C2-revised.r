# code-C2.R
# Revised to match the output style of code-A.R

library(dplyr)
library(tidyr)

run_code_C2 <- function(data_name = "acc_long", out_file = "results-code-C2.txt") {
  
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
    cat("Code C2: Factorial MANOVA on Pairwise Differences\n")
    cat("Output file:", out_file, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' not found.")
    }
    acc_long <- get(data_name, envir = parent.frame(), inherits = TRUE)

    make_pairwise_diff <- function(dat, m1, m2) {
      dat %>%
        filter(method %in% c(m1, m2)) %>%
        pivot_wider(names_from = method, values_from = accuracy) %>%
        mutate(diff = .data[[m1]] - .data[[m2]]) %>%
        select(-all_of(c(m1, m2))) %>%
        pivot_wider(names_from = metric, values_from = diff)
    }

    methods <- c("i", "ii", "iii", "iv", "v")
    pairs <- combn(methods, 2, simplify = FALSE)

    for (p in pairs) {
      m_target <- p[2]
      m_ref <- p[1]
      comp_name <- paste0(m_target, "_vs_", m_ref)
      
      cat("==================================================================\n")
      cat("Factorial Analysis of Difference:", comp_name, "\n")
      cat("==================================================================\n")
      
      diff_mv <- make_pairwise_diff(acc_long, m_target, m_ref)
      
      # Full Factorial formula for the difference vector
      formula_diff <- cbind(rmse, L1, Linf) ~ design_block + G + theta_true +
                      theta_CV + sigma + mean_nsamp + p1 + p2

      fit_diff_manova <- stats::manova(formula_diff, data = diff_mv)

      cat("MANOVA Table (Pillai):\n")
      print(summary(fit_diff_manova, test = "Pillai"))
      cat("\n")
      
      cat("Univariate Breakdown for each Error Metric:\n")
      print(summary.aov(fit_diff_manova))
      cat("\n")
    }

    cat("Analysis completed successfully.\n")
  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n")
    stop(e)
  })
}

run_code_C2()