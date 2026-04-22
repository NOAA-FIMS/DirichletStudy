grid <- expand.grid(
  A=1:3, B=1:3, C=1:3, D=1:3, E=1:3, F=1:3, G=1:3
)
# expand.grid makes A vary fastest by default
grid$Run <- seq_len(nrow(grid))
grid <- grid[, c("Run","A","B","C","D","E","F","G")]

