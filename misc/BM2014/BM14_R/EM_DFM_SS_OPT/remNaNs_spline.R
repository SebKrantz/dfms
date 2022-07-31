
remNaN1o <- function(X, k) {
  d <- dim(X)
  T <- d[1L]
  N <- d[2L]
  k2 <- 2L * k + 1L
  indNaN <- is.na(X)
  wna <- which(indNaN)
  X[wna] <- fmedian(X, TRA = 1L)[wna]
  X_MA <- roll_mean(rbind(X[rep(1L, k), ], X, X[rep(T, k), ]), k2)
  X[wna] <- X_MA[-(1:(k2-1L)), ][wna]
  return(list(X = X, indNaN = indNaN))
}


remNaN2 <- function(X, k, threshold = 0.8) {
  d <- dim(X)
  T <- d[1L]
  N <- d[2L]
  k2 = 2L * k + 1L
  indNaN = is.na(X)
  rem1 = rowSums(indNaN) > N * threshold
  nanLead = cumsum(rem1) == 1:T
  nanEnd = rev(cumsum(rev(rem1)) == 1:T);
  nanLE = nanLead | nanEnd
  X = X[!nanLE, ]
  indNaN = is.na(X)

  for (i in 1:N) {
    x = X[, i]
    nnai = which(!is.na(x))
    ln = length(nnai)
    t1 = nnai[1L]
    t2 = nnai[ln]
    # Cubic spline to interpolate any internal missing values...
    if(ln != t2-t1+1L) x[t1:t2] = spline(nnai, x[nnai], xout = t1:t2)$y
    isnanx = which(is.na(x))
    x[isnanx] = fmedian.default(x)
    x_MA = filter(c(rep(x[1L], k), x, rep(x[T], k)), rep(1, k2)/k2, sides = 1L)
    x[isnanx] = x_MA[-(1:(k2-1L))][isnanx]
    X[, i] = x
  }
  return(list(X = X, indNaN = indNaN))
}

# only remove rows with leading and closing zeros
remNaN3 <- function(X) {
  d <- dim(X)
  T <- d[1L]
  N <- d[2L]
  indNaN = is.na(X)
  rem1 = rowSums(indNaN) == N
  nanLead = cumsum(rem1) == 1:T
  nanEnd = rev(cumsum(rev(rem1)) == 1:T)
  nanLE = nanLead | nanEnd
  X = X[!nanLE, ]
  return(list(X = X, indNaN = is.na(X)))
}

remNaNs <- function(X, options) { # output: [X, indNaN]
  switch(options$method,
         # case 1 % replace all the missing values
         remNaN1o(X, options$k),
         # case 2 % replace missing values after removing leading and closing zeros (= mostly missing rows)
         remNaN2(X, options$k),
         # case 3 % only remove rows with leading and closing zeros
         remNaN3(X),
         stop("remNaN option not implemented yet"))
}


