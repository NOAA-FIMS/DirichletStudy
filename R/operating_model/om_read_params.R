# om_read_params.R
# Reads an input CSV of parameters (parameter,value,comment), creates R objects,
# prints them for validation, and returns a named list.
#
# Usage:
#   source("om_read_params.R")
#   p <- create_and_read_parameters("om_parameters.csv", assign_to_env = TRUE)
#   # This creates variables like Y, A, Linf, K, etc. in the parent frame.

create_and_read_parameters <- function(csv_path = "om_parameters.csv",
                                       assign_to_env = TRUE,
                                       envir = parent.frame(),
                                       verbose = TRUE) {
  if (!file.exists(csv_path)) {
    stop(sprintf("File not found: %s", csv_path))
  }
  # Read CSV (expects columns: parameter, value, comment)
  df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  required_cols <- c("parameter", "value", "comment")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("CSV is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Coerce types: integers for Y, A; numeric for others
  df$value <- suppressWarnings(as.numeric(df$value))
  # Enforce integer for Y and A if present
  if ("Y" %in% df$parameter) df$value[df$parameter == "Y"] <- as.integer(df$value[df$parameter == "Y"])
  if ("A" %in% df$parameter) df$value[df$parameter == "A"] <- as.integer(df$value[df$parameter == "A"])

  # Build named list
  params <- setNames(as.list(df$value), df$parameter)

  # Optionally assign to environment (creates variables Y, A, Linf, ...)
  if (isTRUE(assign_to_env)) {
    list2env(params, envir = envir)
  }

  # Print for validation
  if (isTRUE(verbose)) {
    cat("Parameter values loaded from", csv_path, "\n")
    print(df[, c("parameter", "value", "comment")], row.names = FALSE)
  }

  invisible(params)
}
