## file = run_hake_functional_analysis.R
##
## Run DirichletStudy functional analysis for the linear and saturated
## Dirichlet-multinomial variants using the hake simulation inputs.
##
## Default inputs:
##   hake/hake.inp - simulation settings, including K, h, theta_true, bounds
##   hake/hake.csv - generated hake simulation output with Xg_k count columns
##
## Default outputs:
##   hake/hake_functional_analysis_results.rds
##   hake/hake_functional_analysis_linear_metrics.csv
##   hake/hake_functional_analysis_saturated_metrics.csv
##   hake/hake_functional_analysis_summary.csv

suppressPackageStartupMessages(library(DirichletStudy))

args <- commandArgs(trailingOnly = TRUE)

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(script_path) == 1) dirname(normalizePath(script_path)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."))

arg_or_default <- function(i, default) {
  if (length(args) >= i && nzchar(args[[i]])) args[[i]] else default
}

inp_path <- normalizePath(arg_or_default(1, file.path(script_dir, "hake.inp")), mustWork = TRUE)
hake_csv_path <- normalizePath(arg_or_default(2, file.path(script_dir, "hake.csv")), mustWork = TRUE)
output_dir <- normalizePath(arg_or_default(3, script_dir), mustWork = TRUE)

parse_hake_inp <- function(path) {
  tbl <- read.table(
    file = path,
    sep = "=",
    comment.char = "#",
    strip.white = TRUE,
    blank.lines.skip = TRUE,
    col.names = c("key", "value"),
    colClasses = c("character", "character")
  )
  setNames(tbl$value, tbl$key)
}

parse_numeric_vector <- function(value, n, default) {
  if (is.null(value) || is.na(value) || !nzchar(trimws(value))) {
    return(rep(default, n))
  }
  pieces <- unlist(strsplit(value, "[,[:space:]]+"))
  pieces <- pieces[nzchar(pieces)]
  out <- as.numeric(pieces)
  if (length(out) != n || any(!is.finite(out))) {
    stop(sprintf("Expected %d numeric values, got '%s'.", n, value))
  }
  out
}

nexcom_step <- function(N, K, P, MTC, I, J) {
  if (MTC == FALSE) {
    P[1] <- N
    I <- N
    J <- 0
    if (K != 1) {
      for (ii in 2:K) P[ii] <- 0
      MTC <- (P[K] != N)
      return(list(P = P, MTC = MTC, I = I, J = J))
    }
  }
  if (I > 1) J <- 0
  J <- J + 1
  I <- P[J]
  P[J] <- 0
  P[1] <- I - 1
  P[J + 1] <- P[J + 1] + 1
  MTC <- (P[K] != N)
  list(P = P, MTC = MTC, I = I, J = J)
}

nexcom <- function(N, K) {
  rn <- nexcom_step(N, K, P = integer(K), MTC = FALSE, I = 0, J = 0)
  out <- data.frame(P = rbind(rn$P))
  while (rn$MTC == TRUE) {
    rn <- nexcom_step(N, K, P = rn$P, MTC = rn$MTC, I = rn$I, J = rn$J)
    out <- rbind(out, data.frame(P = rbind(rn$P)))
  }
  out
}

make_open_simplex <- function(K, h, p_lbound, p_ubound) {
  if (!is.numeric(K) || length(K) != 1 || K != as.integer(K) || K < 2) {
    stop("K must be a single integer >= 2.")
  }
  if (!is.numeric(h) || length(h) != 1 || !(h > 0) || !(h < 1 / K)) {
    stop("h must be a single numeric with 0 < h < 1/K.")
  }

  N <- as.integer(ceiling(1 / h))
  simplex <- as.matrix((nexcom(N, K) + 1) / (N + K))
  storage.mode(simplex) <- "double"

  keep <- apply(sweep(simplex, 2, p_ubound, `<=`), 1, all) &
    apply(sweep(simplex, 2, p_lbound, `>=`), 1, all)
  simplex <- simplex[keep, , drop = FALSE]
  colnames(simplex) <- paste0("p", seq_len(K))
  simplex
}

aggregate_hake_counts <- function(path, K) {
  dat <- read.csv(path, check.names = FALSE)
  count_cols <- grep("^X[[:digit:]]+_[[:digit:]]+$", names(dat), value = TRUE)
  if (length(count_cols) == 0) {
    stop("No hake count columns matching X<group>_<category> were found.")
  }

  category <- as.integer(sub("^X[[:digit:]]+_", "", count_cols))
  if (!all(seq_len(K) %in% category)) {
    stop("The hake CSV does not contain count columns for every category.")
  }

  counts <- vapply(seq_len(K), function(k) {
    sum(dat[count_cols[category == k]], na.rm = TRUE)
  }, numeric(1))
  counts <- as.integer(round(counts))
  if (sum(counts) <= 0) {
    stop("Aggregated hake counts sum to zero.")
  }
  counts
}

require_interface_method <- function(object, method, interface_name) {
  method_fun <- tryCatch(do.call("$", list(object, method)), error = function(e) NULL)
  if (!is.function(method_fun)) {
    stop(sprintf(
      "%s does not expose %s(). Install the current DirichletStudy branch before running this script.",
      interface_name, method
    ))
  }
}

write_variant_metrics <- function(result, path) {
  metrics <- result$paper_metrics
  parameter_sets <- as.data.frame(result$parameter_sets)
  names(parameter_sets) <- paste0("p", seq_len(ncol(parameter_sets)))

  out <- cbind(
    parameter_sets,
    data.frame(
      loglik = result$values,
      delta_loglik = metrics$delta_loglik,
      relative_likelihood = metrics$relative_likelihood,
      gradient_norm = metrics$gradient_norm,
      max_abs_gradient = metrics$max_abs_gradient,
      boundary_distance = metrics$boundary_distance,
      l1_distance_from_observed = metrics$l1_distance_from_observed,
      l2_distance_from_observed = metrics$l2_distance_from_observed
    )
  )
  write.csv(out, path, row.names = FALSE)
}

params <- parse_hake_inp(inp_path)

K <- as.integer(params[["K"]])
h <- as.numeric(params[["h"]])
theta <- as.numeric(params[["theta_true"]])
beta <- as.numeric(params[["theta_true"]])
p_lbound <- parse_numeric_vector(params[["p_lbound"]], K, 0)
p_ubound <- parse_numeric_vector(params[["p_ubound"]], K, 1)

simplex <- make_open_simplex(K, h, p_lbound, p_ubound)
counts <- aggregate_hake_counts(hake_csv_path, K)

linear <- new(DirichletLinearInterface)
require_interface_method(linear, "setCounts", "DirichletLinearInterface")
linear$setCounts(counts)
linear$theta <- theta
linear$setSimplexData(simplex)

saturated <- new(DirichletSaturatedInterface)
require_interface_method(saturated, "setCounts", "DirichletSaturatedInterface")
saturated$setCounts(counts)
saturated$beta <- beta
saturated$setSimplexData(simplex)

study <- new(DirichletStudyInterface)
study$addStudy(linear$getId())
study$addStudy(saturated$getId())

ok <- study$runAnalysis()
if (!isTRUE(ok)) stop("Functional analysis failed.")

results <- study$getResults()
variants <- vapply(results, function(x) x$variant, character(1))
linear_result <- results[[which(variants == "linear")[[1]]]]
saturated_result <- results[[which(variants == "saturated")[[1]]]]

artifact <- list(
  input = list(
    inp_path = inp_path,
    hake_csv_path = hake_csv_path,
    K = K,
    h = h,
    p_lbound = p_lbound,
    p_ubound = p_ubound,
    counts = counts,
    observed_proportions = counts / sum(counts),
    theta = theta,
    beta = beta
  ),
  simplex = simplex,
  results = results
)

rds_path <- file.path(output_dir, "hake_functional_analysis_results.rds")
linear_csv_path <- file.path(output_dir, "hake_functional_analysis_linear_metrics.csv")
saturated_csv_path <- file.path(output_dir, "hake_functional_analysis_saturated_metrics.csv")
summary_csv_path <- file.path(output_dir, "hake_functional_analysis_summary.csv")

saveRDS(artifact, rds_path)
write_variant_metrics(linear_result, linear_csv_path)
write_variant_metrics(saturated_result, saturated_csv_path)

summary <- data.frame(
  variant = c("linear", "saturated"),
  dispersion = c(linear_result$dispersion, saturated_result$dispersion),
  effective_sample_size = c(linear_result$effective_sample_size, saturated_result$effective_sample_size),
  max_loglik = c(linear_result$max_value, saturated_result$max_value),
  min_loglik = c(linear_result$min_value, saturated_result$min_value),
  boundary_gradient_ratio = c(
    linear_result$paper_metrics$boundary_gradient_ratio,
    saturated_result$paper_metrics$boundary_gradient_ratio
  ),
  mean_gradient_norm = c(
    mean(linear_result$paper_metrics$gradient_norm, na.rm = TRUE),
    mean(saturated_result$paper_metrics$gradient_norm, na.rm = TRUE)
  ),
  mean_abs_derivative_step = c(
    mean(linear_result$paper_metrics$mean_abs_derivative_step, na.rm = TRUE),
    mean(saturated_result$paper_metrics$mean_abs_derivative_step, na.rm = TRUE)
  )
)
write.csv(summary, summary_csv_path, row.names = FALSE)

cat("Functional analysis complete.\n")
cat("Input file: ", inp_path, "\n", sep = "")
cat("Hake CSV: ", hake_csv_path, "\n", sep = "")
cat("Simplex rows: ", nrow(simplex), "\n", sep = "")
cat("Counts: ", paste(counts, collapse = ", "), "\n", sep = "")
cat("Results RDS: ", rds_path, "\n", sep = "")
cat("Linear metrics CSV: ", linear_csv_path, "\n", sep = "")
cat("Saturated metrics CSV: ", saturated_csv_path, "\n", sep = "")
cat("Summary CSV: ", summary_csv_path, "\n", sep = "")
