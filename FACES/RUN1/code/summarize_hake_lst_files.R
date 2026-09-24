# summarize_hake_lst_files.R
# Build a compact summary from ex1...ex972/*.lst files

extract_first <- function(x, pattern, group = 1, perl = TRUE) {
  m <- regexec(pattern, x, perl = perl)
  r <- regmatches(x, m)
  out <- sapply(r, function(z) if (length(z) >= group + 1) z[group + 1] else NA_character_)
  out
}

extract_block <- function(lines, start_pattern, end_pattern = NULL) {
  s <- grep(start_pattern, lines, ignore.case = TRUE)
  if (length(s) == 0) return(character(0))
  s <- s[1]
  if (is.null(end_pattern)) return(lines[s:length(lines)])
  e <- grep(end_pattern, lines[(s + 1):length(lines)], ignore.case = TRUE)
  if (length(e) == 0) return(lines[s:length(lines)])
  e <- s + e[1] - 1
  lines[s:e]
}

parse_emmeans_metric <- function(lines, metric_name) {
  # Look for an emmeans table that contains method + emmean for the given metric
  idx <- grep(metric_name, lines, ignore.case = TRUE)
  if (length(idx) == 0) {
    return(data.frame(method = paste0("(", 1:5, ")"), emmean = NA_real_))
  }
  
  # Search forward from first mention for rows beginning with method labels i-v or 1-5
  start <- idx[1]
  search_lines <- lines[start:min(length(lines), start + 80)]
  
  # Typical formats:
  # i   0.115 ...
  # (i) 0.115 ...
  # 1   0.115 ...
  pat_list <- c(
    "^\\s*\\(?([ivx]+)\\)?\\s+([0-9.]+(?:e[-+]?\\d+)?)\\b",
    "^\\s*\\(?([1-5])\\)?\\s+([0-9.]+(?:e[-+]?\\d+)?)\\b"
  )
  
  vals <- data.frame(method = character(), emmean = numeric())
  
  for (pat in pat_list) {
    for (ln in search_lines) {
      m <- regexec(pat, ln, perl = TRUE, ignore.case = TRUE)
      rr <- regmatches(ln, m)[[1]]
      if (length(rr) >= 3) {
        meth_raw <- tolower(rr[2])
        meth <- switch(meth_raw,
                       "i" = "(i)", "ii" = "(ii)", "iii" = "(iii)",
                       "iv" = "(iv)", "v" = "(v)",
                       "1" = "(i)", "2" = "(ii)", "3" = "(iii)",
                       "4" = "(iv)", "5" = "(v)",
                       NA_character_)
        if (!is.na(meth)) {
          vals <- rbind(vals, data.frame(method = meth, emmean = as.numeric(rr[3])))
        }
      }
    }
    vals <- vals[!duplicated(vals$method), , drop = FALSE]
    if (nrow(vals) >= 5) break
  }
  
  target_methods <- paste0("(", c("i","ii","iii","iv","v"), ")")
  out <- merge(data.frame(method = target_methods), vals, by = "method", all.x = TRUE, sort = FALSE)
  out
}

parse_test_stat <- function(lines, metric_name, test_name) {
  idx <- grep(metric_name, lines, ignore.case = TRUE)
  if (length(idx) == 0) return(data.frame(stat = NA_real_, df = NA_real_, p = NA_character_))
  search_lines <- lines[idx[1]:min(length(lines), idx[1] + 120)]
  
  txt <- paste(search_lines, collapse = "\n")
  
  if (tolower(test_name) == "friedman") {
    stat <- extract_first(txt, "Friedman[^\\n]*?(?:chi-squared|chi\\^?2)\\s*=\\s*([0-9.]+)")
    if (all(is.na(stat))) stat <- extract_first(txt, "Friedman[^\\n]*?=\\s*([0-9.]+)")
    df   <- extract_first(txt, "Friedman[^\\n]*?df\\s*=\\s*([0-9.]+)")
    p    <- extract_first(txt, "Friedman[^\\n]*?p\\s*[<=]\\s*([^\\s,\\n]+)")
    return(data.frame(stat = as.numeric(stat), df = as.numeric(df), p = p))
  }
  
  if (tolower(test_name) == "kendall") {
    stat <- extract_first(txt, "Kendall'?s\\s+W\\s*=\\s*([0-9.]+)")
    return(data.frame(stat = as.numeric(stat), df = NA_real_, p = NA_character_))
  }
  
  if (tolower(test_name) == "anova") {
    stat <- extract_first(txt, "(?:ANOVA|Analysis of Variance)[^\\n]*?F\\s*=\\s*([0-9.]+)")
    if (all(is.na(stat))) stat <- extract_first(txt, "\\bF\\s*=\\s*([0-9.]+)")
    p    <- extract_first(txt, "(?:ANOVA|Analysis of Variance)[^\\n]*?p\\s*[<=]\\s*([^\\s,\\n]+)")
    return(data.frame(stat = as.numeric(stat), df = NA_real_, p = p))
  }
  
  data.frame(stat = NA_real_, df = NA_real_, p = NA_character_)
}

detect_distribution <- function(lines, path) {
  txt <- paste(lines[1:min(length(lines), 120)], collapse = "\n")
  low <- tolower(txt)
  
  if (grepl("poisson", low)) return("Poisson")
  if (grepl("lognormal|log-normal", low)) return("Lognormal")
  if (grepl("negative binomial|neg binomial|nb\\b", low)) return("Negative binomial")
  if (grepl("log-uniform|log uniform", low)) return("Log-uniform")
  
  # fallback from file number if experiment ordering is 4-way repeating
  ex_num <- as.integer(sub(".*ex([0-9]+).*", "\\1", path, ignore.case = TRUE))
  if (!is.na(ex_num)) {
    r <- ((ex_num - 1) %% 4) + 1
    return(c("Poisson", "Lognormal", "Negative binomial", "Log-uniform")[r])
  }
  
  NA_character_
}

parse_one_lst <- function(f) {
  lines <- tryCatch(readLines(f, warn = FALSE, encoding = "UTF-8"),
                    error = function(e) readLines(f, warn = FALSE))
  dist_type <- detect_distribution(lines, f)
  
  rmse_tab <- parse_emmeans_metric(lines, "RMSE")
  l1_tab   <- parse_emmeans_metric(lines, "L1")
  linf_tab <- parse_emmeans_metric(lines, "Linf")
  
  rmse_best <- rmse_tab$method[which.min(rmse_tab$emmean)][1]
  l1_best   <- l1_tab$method[which.min(l1_tab$emmean)][1]
  linf_best <- linf_tab$method[which.min(linf_tab$emmean)][1]
  
  rmse_fried <- parse_test_stat(lines, "RMSE", "friedman")
  rmse_kw    <- parse_test_stat(lines, "RMSE", "kendall")
  rmse_an    <- parse_test_stat(lines, "RMSE", "anova")
  
  l1_fried <- parse_test_stat(lines, "L1", "friedman")
  l1_kw    <- parse_test_stat(lines, "L1", "kendall")
  l1_an    <- parse_test_stat(lines, "L1", "anova")
  
  linf_fried <- parse_test_stat(lines, "Linf", "friedman")
  linf_kw    <- parse_test_stat(lines, "Linf", "kendall")
  linf_an    <- parse_test_stat(lines, "Linf", "anova")
  
  get_metric_vals <- function(tab, prefix) {
    x <- setNames(tab$emmean, tab$method)
    out <- c(
      x["(i)"], x["(ii)"], x["(iii)"], x["(iv)"], x["(v)"]
    )
    names(out) <- paste0(prefix, "_", c("i","ii","iii","iv","v"))
    out
  }
  
  ex_num <- as.integer(sub(".*ex([0-9]+).*", "\\1", f, ignore.case = TRUE))
  
  data.frame(
    experiment = ex_num,
    file = basename(f),
    distribution = dist_type,
    rmse_best = rmse_best,
    l1_best = l1_best,
    linf_best = linf_best,
    
    t(get_metric_vals(rmse_tab, "rmse")),
    t(get_metric_vals(l1_tab,   "l1")),
    t(get_metric_vals(linf_tab, "linf")),
    
    rmse_friedman_chisq = rmse_fried$stat,
    rmse_friedman_df    = rmse_fried$df,
    rmse_friedman_p     = rmse_fried$p,
    rmse_kendall_w      = rmse_kw$stat,
    rmse_anova_f        = rmse_an$stat,
    rmse_anova_p        = rmse_an$p,
    
    l1_friedman_chisq   = l1_fried$stat,
    l1_friedman_df      = l1_fried$df,
    l1_friedman_p       = l1_fried$p,
    l1_kendall_w        = l1_kw$stat,
    l1_anova_f          = l1_an$stat,
    l1_anova_p          = l1_an$p,
    
    linf_friedman_chisq = linf_fried$stat,
    linf_friedman_df    = linf_fried$df,
    linf_friedman_p     = linf_fried$p,
    linf_kendall_w      = linf_kw$stat,
    linf_anova_f        = linf_an$stat,
    linf_anova_p        = linf_an$p,
    
    stringsAsFactors = FALSE
  )
}

# Find all .lst files under ex1...ex972
lst_files <- list.files(
  path = ".",
  pattern = "\\.lst$",
  recursive = TRUE,
  full.names = TRUE
)

lst_files <- lst_files[grepl("ex[0-9]+", lst_files, ignore.case = TRUE)]

cat("Found", length(lst_files), ".lst files\n")

summary_df <- do.call(
  rbind,
  lapply(lst_files, function(f) {
    cat("Parsing:", f, "\n")
    tryCatch(parse_one_lst(f),
             error = function(e) {
               message("Failed on ", f, ": ", conditionMessage(e))
               NULL
             })
  })
)

summary_df <- summary_df[order(summary_df$experiment), ]

write.csv(summary_df, "hake_972_experiment_summary.csv", row.names = FALSE)

# Also write frequency tables by distribution
make_freq <- function(df, metric_col) {
  as.data.frame.matrix(table(df$distribution, df[[metric_col]]))
}

write.csv(make_freq(summary_df, "rmse_best"), "hake_rmse_best_frequency_by_distribution.csv")
write.csv(make_freq(summary_df, "l1_best"),   "hake_l1_best_frequency_by_distribution.csv")
write.csv(make_freq(summary_df, "linf_best"), "hake_linf_best_frequency_by_distribution.csv")

cat("Wrote:\n")
cat("  hake_972_experiment_summary.csv\n")
cat("  hake_rmse_best_frequency_by_distribution.csv\n")
cat("  hake_l1_best_frequency_by_distribution.csv\n")
cat("  hake_linf_best_frequency_by_distribution.csv\n")