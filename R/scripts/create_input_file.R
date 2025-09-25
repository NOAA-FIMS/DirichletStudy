# file = create_input_file.R
# Create DM_rmse_All_3_simplex.inp with defaults (run once or keep in script)
inp_path <- "DM_rmse_All_3_simplex.inp"
writeLines(c(
  "# Input parameters for DM_rmse_All_3_simplex",
  "K = 3",
  "G = 10",
  "h = 0.025",
  "theta_true = 0.30",
  "Nmin = 50",
  "Nmax = 500",
  "nsamples = 1",
  "random.seed = 1097"
), con = inp_path)
