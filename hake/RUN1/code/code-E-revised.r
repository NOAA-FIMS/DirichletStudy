# code-E.R
# Revised to match the output style of code-A.R

library(mgcv)
library(dplyr)

run_code_E <- function(data_name = "acc_long", out_file = "results-code-E.txt") {
  
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
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' not found.")
    }
    acc_long <- get(data_name, envir = parent.frame(), inherits = TRUE)

    # Global GAM: Spatial error surface by method
    cat("Fitting Global GAM with spatial splines s(p1, p2) by method...\n")
    fit_gam <- gam(
      log1p(accuracy) ~ metric * method +
        s(p1, p2, by = method, k = 30) +
        design_block + G + theta_true + theta_CV + sigma + mean_nsamp,
      data = acc_long,
      method = "REML"
    )

    cat("Global GAM Summary:\n")
    print(summary(fit_gam))
    cat("\nGAM ANOVA Table:\n")
    print(anova(fit_gam))
    cat("\n")

    # Difference GAM: Comparing Method V vs Method IV
    cat("Analyzing difference between Method 'v' and Method 'iv'...\n")
    
    # Internal helper to generate diff
    diff_v_iv <- acc_long %>%
      filter(method %in% c("v", "iv")) %>%
      tidyr::pivot_wider(names_from = method, values_from = accuracy) %>%
      mutate(diff = v - iv)

    fit_diff_gam <- gam(
      log1p(abs(diff)) ~ metric +
        s(p1, p2, by = metric, k = 30) +
        design_block + G + theta_true + theta_CV + sigma + mean_nsamp,
      data = diff_v_iv,
      method = "REML"
    )

    cat("Difference GAM Summary (Method v - iv):\n")
    print(summary(fit_diff_gam))
    cat("\n")

    cat("Analysis completed successfully.\n")
    
    invisible(list(fit_gam = fit_gam, fit_diff_gam = fit_diff_gam))

  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n")
    stop(e)
  })
}

run_code_E()