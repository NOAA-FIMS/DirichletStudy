# file = test_build_write_G123.R

source("build_hake_factorial_design_matrix_revised_G123.R")

design <- build_hake_factorial_design_matrix_default(
  output_csv = "hake_factorial_design_matrix.csv",
  legend_csv = "hake_factorial_design_legend.csv",
  run_verification = FALSE
)

source("write_hake_inp_files_from_matrix_unique_seed_revised_G123.R")

manifest <- write_hake_inp_files_from_matrix(
  matrix_file = "hake_factorial_design_matrix.csv",
  output_dir = ".",
  validate_common_G = TRUE,
  random_seed_generator_seed = 11131
)

table(manifest$G)
length(unique(manifest$random.seed))
nrow(manifest)