
lagnam <- function(nam, p) list(nam, as.vector(t(outer(paste0("L", seq_len(p)), nam, paste, sep = "."))))

msum <- function(x) {
  stats <- qsu(x)
  med <- fmedian(x)
  res <- if(is.matrix(x)) cbind(stats[, 1:2], Median = med, stats[, -(1:2)]) else
    c(stats[1:2], Median = med, stats[-(1:2)])
  class(res) <- "qsu"
  res
}

AC1 <- function(res, anymiss) {
  ACF <- cov(res[-1L, ], res[-nrow(res), ],
             use = if(anymiss) "pairwise.complete.obs" else "everything")
  diag(ACF) / fvar(res)
}

unscale <- function(x, stats) TRA.matrix(TRA.matrix(x, stats[, "SD"], "*"), stats[, "Mean"], "+", set = TRUE)

ftail <- function(x, p) {n <- dim(x)[1L]; x[(n-p+1L):n, , drop = FALSE]}

#' (Fast) Barebones Vector-Autoregression
#'
#' Quickly estimate an VAR(p) model using Armadillo's inverse function.
#'
#' @param x data matrix with time series in columns - without missing values.
#' @param p integer. The lag order of the VAR.
#'
#' @returns A list containing matrices \code{Y = x[-(1:p), ]}, \code{X} which contains lags 1 - p of \code{x} combined column-wise,
#' \code{A} which is the np x n transition matrix, where n is the number of series in \code{x}, and the VAR residual matrix \code{res = Y - X \%*\% A}.
#'
#' @export
.VAR <- function(x, p = 1L) {
  T <- dim(x)[1L]
  Y <- x[(p + 1L):T, ]
  X <- do.call(cbind, lapply(1:p, function(i) x[(p + 1L - i):(T - i), ]))
  # A <- qr.coef(qr(X), Y) # solve(t(X) %*% X) %*% t(X) %*% Y
  A <- ainv(crossprod(X)) %*% crossprod(X, Y) # Faster !!!

  return(list(Y = Y, X = X, A = A, res = Y - X %*% A))
}

# ginv <- MASS::ginv # use apinv


#' Convergence Test for EM-Algorithm
#'
#' @param loglik vurrent value of the log-likelihood function.
#' @param previous_loglik value of the log-likelihood function at the previous iteration.
#' @param tol the numerical tolerance of the test. If |LL(t) - LL(t-1)| / avg < tol, where avg = (|LL(t)| + |LL(t-1)|)/2, then algorithm has converged.
#' @param check.increased logical. Check if likelihood has increased.
#' @return A logical statement indicating whether EM algorithm has converged. if \code{check.increased = TRUE}, a vector with 2 elements indicating the convergence status and whether the likelihood has decreased.
#' @export
em_converged <- function(loglik, previous_loglik, tol = 1e-4, check.increased = TRUE) { # [converged, decrease]
  # EM_CONVERGED Has EM converged?
  # [converged, decrease] = em_converged(loglik, previous_loglik, tol)
  #
  # We have converged if the slope of the log-likelihood function falls below 'tol',
  # i.e., |f(t) - f(t-1)| / avg < tol,
  # where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
  # 'tol' defaults to 1e-4.
  #
  # This stopping criterion is from Numerical Recipes in C p423
  #
  # If we are doing MAP estimation (using priors), the likelihood can decrase,
  # even though the mode of the posterior is increasing.

 delta_loglik = abs(loglik - previous_loglik)
 avg_loglik = (abs(loglik) + abs(previous_loglik) + .Machine$double.eps)/2
 test = delta_loglik / avg_loglik < tol
 converged = is.finite(test) && test

 if(check.increased) {
    test <- loglik - previous_loglik < -1e-3
    if(is.finite(test) && test) { # allow for a little imprecision
      sprintf('******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik)
      decrease = TRUE
    } else decrease = FALSE
    return(c(converged = converged, decrease = decrease))
 }

 converged
}


findNA_LE <- function(W, threshold) {
  d <- dim(W)
  st <- 1:d[1L]
  rem1 <- rowSums(W) > d[2L] * threshold
  nanLead <- cumsum(rem1) == st
  nanEnd <- rev(cumsum(rev(rem1)) == st)
  return(nanLead | nanEnd)
}


impNA_MA <- function(X, W, k) {
  d <- dim(X)
  T <- d[1L]
  k2 <- 2L * k + 1L
  weights <- rep(1, k2)/k2
  sss <- -(1:(k2-1L))
  for (i in 1:d[2L]) {
    x <- X[, i]
    nai <- which(W[, i])
    x[nai] <- fmedian.default(x)
    x_MA <- filter(c(rep(x[1L], k), x, rep(x[T], k)), weights, sides = 1L)
    x[nai] <- x_MA[sss][nai]
    X[, i] <- x
  }
  return(X)
}


impNA_spline <- function(X, W, k) {
  d <- dim(X)
  T <- d[1L]
  k2 <- 2L * k + 1L
  weights <- rep(1, k2)/k2
  sss <- -(1:(k2-1L))

  for (i in 1:d[2L]) {
    x = X[, i]
    nnai <- which(!W[, i])
    ln = length(nnai)
    t1 = nnai[1L]
    t2 = nnai[ln]
    # Cubic spline to interpolate any internal missing values...
    if(ln != t2-t1+1L) x[t1:t2] = spline(nnai, x[nnai], xout = t1:t2)$y
    isnanx = which(is.na(x))
    x[isnanx] = fmedian.default(x)
    x_MA = filter(c(rep(x[1L], k), x, rep(x[T], k)), weights, sides = 1L)
    x[isnanx] = x_MA[sss][isnanx]
    X[, i] = x
  }
  return(X)
}

#' Remove and Impute Missing Values in a Multivariate Time Series
#'
#' This function imputes missing (and infinite) values in a stationary multivariate time series using various
#' methods, and removes cases with too many missing values.
#'
#' @param X a matrix or multivariate time series where each column is a series.
#' @inheritParams DFM
#'
#' @returns A list with the imputed matrix \code{X_imp}, a missingness matrix \code{W} matching the dimensions of \code{X_imp},
#' and a vector or cases \code{na.rm} indicating cases with too many missing values that were removed beforehand.
#' @export
tsremimpNA <- function(X,
                       max.missing = 0.5,
                       na.rm.method = c("LE", "all"),
                       na.impute = c("median", "rnorm", "median.ma", "median.ma.spline"),
                       ma.terms = 3L) {
  W <- !is.finite(X) # is.na(X)
  n <- dim(X)[2L]
  na.rm <- NULL
  if(max.missing < 1) {
    thresh <- switch(na.rm.method[1L],
                     LE = findNA_LE(W, max.missing),
                     all = rowSums(W) > max.missing * n,
                     stop("Unknown na.rm.method:", na.rm.method[1L]))
    if(any(thresh)) {
      na.rm <- which(thresh)
      X <- X[-na.rm, ]
      W <- W[-na.rm, ]
    }
  }
  list(X_imp = switch(na.impute[1L],
                  median = replace(X, W, fmedian(X, TRA = 1L)[W]),
                  rnorm = replace(X, W, rnorm(sum(W))),
                  median.ma = impNA_MA(X, W, ma.terms),
                  median.ma.spline = impNA_spline(X, W, ma.terms),
                  stop("Unknown na.impute option:", na.impute[1L])),
       W = W,
       na.rm = na.rm)
}

