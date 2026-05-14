# code-G.R
# Revised to match the output style of code-A.R

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

run_code_G <- function(data_name = "acc_long", out_file = "results-code-G.txt") {
  
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
    cat("Code G: Rank-Based Accuracy Analysis\n")
    cat("Output file:", out_file, "\n")
    cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

    if (!exists(data_name, envir = parent.frame(), inherits = TRUE)) {
      stop("Required data object '", data_name, "' not found.")
    }
    acc_long <- get(data_name, envir = parent.frame(), inherits = TRUE)

    cat("Computing ranks within dataset_id and metric...\n")
    rank_dat <- acc_long %>%
      group_by(dataset_id, metric) %>%
      mutate(rank = rank(accuracy, ties.method = "average")) %>%
      ungroup()

    cat("Descriptive Statistics: Mean Rank and Probability of being Best (Rank 1)\n")
    cat("----------------------------------------------------------------------\n")
    rank_summary <- rank_dat %>%
      group_by(metric, method) %>%
      summarize(
        mean_rank = mean(rank, na.rm = TRUE),
        prob_best = mean(rank == 1, na.rm = TRUE),
        .groups = "drop"
      )
    print(as.data.frame(rank_summary))
    cat("\n")

    cat("Fitting Linear Mixed Model on Ranks...\n")
    fit_rank <- lmer(
      rank ~ method * metric +
        design_block + G + theta_true + theta_CV + sigma + mean_nsamp +
        p1 + p2 +
        (1 | example_id) + (1 | dataset_id),
      data = rank_dat
    )

    cat("LMM ANOVA Table (Satterthwaite approximation):\n")
    print(anova(fit_rank))
    cat("\n")

    cat("Pairwise Comparisons of Ranks (Method | Metric):\n")
    cat("------------------------------------------------\n")
    pw_comp <- emmeans(fit_rank, pairwise ~ method | metric, adjust = "holm")
    print(pw_comp)
    cat("\n")

    cat("Analysis completed successfully.\n")
    
    invisible(list(rank_summary = rank_summary, fit_rank = fit_rank))

  }, error = function(e) {
    cat("ERROR\n-----\n", conditionMessage(e), "\n")
    stop(e)
  })
}

run_code_G()