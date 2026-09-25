# ==========================================================
# Delete hake_exN.csv files from folders ex1 to ex729
# ==========================================================

for (i in 1:729) {
  
  file_to_delete <- file.path(
    paste0("ex", i),
    paste0("hake_ex", i, ".csv")
  )
  
  if (file.exists(file_to_delete)) {
    file.remove(file_to_delete)
    cat("Deleted:", file_to_delete, "\n")
  } else {
    cat("Not found:", file_to_delete, "\n")
  }
}

cat("Finished.\n")