# file = Exclude the Log-uniform distribution.R
# First, builds the full 4-block matrix by distribution
# Second, removes experiments from the Log-uniform 
# distribution block (4) for randomly composition sizes

source("build_hake_factorial_design_matrix.R")

# Build the full 4-block matrix first
design_full <- build_hake_factorial_design_matrix_default(
  output_csv = "hake_factorial_design_matrix_full_4blocks.csv",
  legend_csv = "hake_factorial_design_legend_full_4blocks.csv",
  run_verification = FALSE
)

# Exclude Log-uniform
design_3block <- subset(design_full, design_block != "Log-uniform")

# Renumber example IDs and input-file names so downstream scripts expect hake_ex1 ... hake_ex729
design_3block$example_id <- seq_len(nrow(design_3block))
design_3block$inp_file <- sprintf("hake_ex%d.inp", design_3block$example_id)

# Write the filtered matrix
write.csv(
  design_3block,
  "hake_factorial_design_matrix.csv",
  row.names = FALSE,
  na = ""
)

# Check the result
table(design_3block$design_block)
nrow(design_3block)