setCN <- function(x, nam) `dimnames<-`(x, list(NULL, nam))

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

unscale <- function(x, stats) TRA.matrix(TRA.matrix(x, stats[, "SD"], "*"), stats[, "Mean"], "+")

ftail <- function(x, p) {n <- dim(x)[1L]; x[(n-p+1L):n, , drop = FALSE]}

#' Fast Vector-Autoregression
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
fVAR <- function(x, p = 1L) {
  T <- dim(x)[1L]
  Y <- x[(p + 1L):T, ]
  X <- do.call(cbind, lapply(1:p, function(i) x[(p + 1L - i):(T - i), ]))
  # A <- qr.coef(qr(X), Y) # solve(t(X) %*% X) %*% t(X) %*% Y
  A <- ainv(crossprod(X)) %*% crossprod(X, Y) # Faster !!!

  return(list(Y = Y, X = X, A = A, res = Y - X %*% A))
}

# ginv <- MASS::ginv # use apinv


#' Convergence test for EM-algorithm.
#'
#' @param loglik Current value of the log-likelihood function
#' @param previous_loglik Value of the log-likelihood function at the previous
#  iteration
#' @param threshold If difference is less than threshold, then algorithm has
#' converged
#' @return A logical statement indicating whether EM algorithm has converged
#' according to slope convergence test
em_converged <- function(loglik, previous_loglik, threshold = 1e-4) { # , check_increased = TRUE

  # converged <- FALSE
  # decrease <- 0
  # if (check_increased) {
  #   if (loglik - previous_loglik < -0.001) {
  #     #            message("*** Likelihood decreased from ", previous_loglik, " to ", loglik, "\n")
  #     decrease <- 1
  #   }
  # }
  if(is.finite(loglik) && is.finite(previous_loglik)) {
    delta_loglik <- abs(loglik - previous_loglik)
    avg_loglik <- (abs(loglik) + abs(previous_loglik) + .Machine$double.eps)/2
    if(delta_loglik/avg_loglik < threshold) return(TRUE)
  }
  return(FALSE)
  # return(list('converged'=converged, 'decrease'=decrease))
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
                       na.impute = c("median", "rnrom", "median.ma", "median.ma.spline"),
                       na.impute.MA = 3L) {
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
                  rnrom = replace(X, W, rnorm(sum(W))),
                  median.ma = impNA_MA(X, W, na.impute.MA),
                  median.ma.spline = impNA_spline(X, W, na.impute.MA),
                  stop("Unknown na.impute option:", na.impute[1L])),
       W = W,
       na.rm = na.rm)
}

