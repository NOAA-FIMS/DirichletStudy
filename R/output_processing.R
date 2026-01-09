
#' Compute lagged differences of a numeric vector
lagged_diff <- function(x, lag = 1) {
  if (length(x) <= lag) return(numeric(0))
  x[(lag + 1):length(x)] - x[1:(length(x) - lag)]
}

#' Compute nearest neighbor smoothness of a function over points in space
nn_smoothness <- function(X, f, scale = NULL, eps = 1e-8) {
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(scale)) {
    scale <- apply(X, 2, function(z) diff(range(z)))
    scale[scale == 0] <- 1
  }

  Xs <- sweep(X, 2, scale, "/")

  slopes <- numeric(n)

  for (k in seq_len(n)) {
    d <- sqrt(rowSums((Xs[k, ] - Xs)^2))
    d[k] <- Inf
    j <- which.min(d)

    slopes[k] <- (f[k] - f[j]) / d[j]
  }

  m <- mean(slopes)
  denom <- if (abs(m) < eps) mean(abs(slopes)) else abs(m)

  list(
    slopes = slopes,
    smoothness_cv = sd(slopes) / denom
  )
}


#' Compute nearest neighbor smoothness for each column of a gradient matrix
nn_gradient_smoothness <- function(X, grad, scale = NULL) {
  p <- ncol(grad)
  out <- vector("list", p)

  for (i in seq_len(p)) {
    out[[i]] <- nn_smoothness(X, grad[, i], scale)
  }

  names(out) <- colnames(grad)
  out
}

#' Assemble output data frame from evaluation results
assemble_output <- function(eval) {
  df <- as.data.frame(eval$X)
  colnames(df) <- paste0("x", seq_len(ncol(df)))
  df$f <- eval$f

  if (!is.null(eval$grad)) {
    df <- cbind(df, as.data.frame(eval$grad))
  }

  df
}