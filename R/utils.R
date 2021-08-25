VAR <- function(x, p) {
  T <- nrow(x)
  Y <- x[(p + 1L):T, ]
  X <- NULL
  for (i in 1:p) X <- cbind(X, x[(p + 1L - i):(T - i), ])
  A <- qr.coef(qr(X), Y) # solve(t(X) %*% X) %*% t(X) %*% Y
  res <- Y - X %*% A

  return(list(Y = Y, X = X, A = A, res = res))
}

ginv <- MASS::ginv
