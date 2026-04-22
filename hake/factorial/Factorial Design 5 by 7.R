grid <- expand.grid(
  A=1:5, B=1:5, C=1:5, D=1:5, E=1:5, F=1:5, G=1:5
)
# expand.grid makes A vary fastest by default
grid$Run <- seq_len(nrow(grid))
grid <- grid[, c("Run","A","B","C","D","E","F","G")]
