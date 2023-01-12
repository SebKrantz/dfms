#' Generates Dimension names for the Factor Transition Matrix A
#' @param nam a vector of names.
#' @param p integer. The lag order.
#' @noRd
lagnam <- function(nam, p) list(nam, as.vector(t(outer(paste0("L", seq_len(p)), nam, paste, sep = "."))))

#' Computes Summary Statistics: collapse::qsu() + median
#' @param x a numeric vector, matrix or data frame like object.
#' @noRd
msum <- function(x) {
  stats <- qsu(x)
  med <- fmedian(x)
  res <- if(is.matrix(stats))
    cbind(stats[, 1:2, drop = FALSE], Median = med, stats[, -(1:2), drop = FALSE]) else
    c(stats[1:2], Median = med, stats[-(1:2)])
  class(res) <- "qsu"
  res
}

#' Computes AR(1) Autocorrelation Coefficient (of DFM residuals)
#' @param res a numeric matrix.
#' @param anymiss logical indicating whether \code{res} contains any missing values.
#' @noRd
AC1 <- function(res, anymiss) {
  ACF <- cov(res[-1L, ], res[-nrow(res), ],
             use = if(anymiss) "pairwise.complete.obs" else "everything")
  diag(ACF) / fvar(res)
}

#' Undo Standardizing
#' @param x a vector, matrix or data frame like object.
#' @param stats a statistics matrix of class 'qsu'.
#' @noRd
unscale <- function(x, stats) TRA.matrix(TRA.matrix(x, stats[, "SD"], "*"), stats[, "Mean"], "+", set = TRUE)

#' Faster tail() Function: For use in forecast.dfm()
#' @param x a matrix.
#' @param p the number of last rows to take from \code{x}.
#' @noRd
ftail <- function(x, p) {n <- dim(x)[1L]; x[(n-p+1L):n, , drop = FALSE]}

#' (Fast) Barebones Vector-Autoregression
#'
#' Quickly estimate a VAR(p) model using Armadillo's inverse function.
#'
#' @param x data numeric matrix with time series in columns - without missing values.
#' @param p positive integer. The lag order of the VAR.
#'
#' @returns A list containing matrices \code{Y = x[-(1:p), ]}, \code{X} which contains lags 1 - p of \code{x} combined column-wise,
#' \code{A} which is the \eqn{np \times n}{np x n} transition matrix, where n is the number of series in \code{x}, and the VAR residual matrix \code{res = Y - X \%*\% A}.
#'
#' @returns A list with the following elements:
#'  \item{\code{Y}}{\code{x[-(1:p), ]}. }
#'  \item{\code{X}}{lags 1 - p of \code{x} combined column-wise. }
#'  \item{\code{A}}{\eqn{np \times n}{np x n} transition matrix, where n is the number of series in \code{x}. }
#'  \item{\code{res}}{VAR residual matrix: \code{Y - X \%*\% A}. }
#'
#' @examples
#' var = .VAR(diff(EuStockMarkets), 3)
#' str(var)
#' var$A
#' rm(var)
#'
#' @export
.VAR <- function(x, p = 1L) {
  TT <- dim(x)[1L]
  Y <- x[(p + 1L):TT, ]
  X <- do.call(cbind, lapply(1:p, function(i) x[(p + 1L - i):(TT - i), ]))
  # A <- qr.coef(qr(X), Y) # solve(t(X) %*% X) %*% t(X) %*% Y
  A <- ainv(crossprod(X)) %*% crossprod(X, Y) # Faster !!!

  return(list(Y = Y, X = X, A = A, res = Y - X %*% A))
}

# ginv <- MASS::ginv # use apinv


#' Convergence Test for EM-Algorithm
#'
#' @param loglik numeric. Current value of the log-likelihood function.
#' @param previous_loglik numeric. Value of the log-likelihood function at the previous iteration.
#' @param tol numerical. The tolerance of the test. If |LL(t) - LL(t-1)| / avg < tol, where avg = (|LL(t)| + |LL(t-1)|)/2, then algorithm has converged.
#' @param check.increased logical. Check if likelihood has increased.
#' @return A logical statement indicating whether EM algorithm has converged. if \code{check.increased = TRUE}, a vector with 2 elements indicating the convergence status and whether the likelihood has decreased.
#'
#' @examples
#' em_converged(1001, 1000)
#' em_converged(10001, 10000)
#' em_converged(10001, 10000, check = TRUE)
#' em_converged(10000, 10001, check = TRUE)
#' @export
em_converged <- function(loglik, previous_loglik, tol = 1e-4, check.increased = FALSE) { # [converged, decrease]
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

 delta_loglik <- abs(loglik - previous_loglik)
 avg_loglik <- (abs(loglik) + abs(previous_loglik) + .Machine$double.eps)/2
 test <- delta_loglik / avg_loglik < tol
 converged <- is.finite(test) && test

 if(check.increased) {
    test <- loglik - previous_loglik < -1e-3
    if(is.finite(test) && test) { # allow for a little imprecision
      sprintf('******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik)
      decrease <- TRUE
    } else decrease <- FALSE
    return(c(converged = converged, decrease = decrease))
 }

 converged
}

#' Detect Missing Values at the Beginning and End of the Sample
#' @param W a logical matrix indicating missingness: the result of \code{\link{is.na}}.
#' @param threshold the proportion of values missing in a row for the case to be considered missing.
#' @noRd
findNA_LE <- function(W, threshold) {
  d <- dim(W)
  st <- 1:d[1L]
  rem1 <- rowSums(W) > d[2L] * threshold
  nanLead <- cumsum(rem1) == st
  nanEnd <- rev(cumsum(rev(rem1)) == st)
  return(nanLead | nanEnd)
}

#' Impute Time Series with Moving Average
#' @param X numeric matrix representing a multivariate time series.
#' @param W a logical matrix indicating missingness: the result of \code{\link{is.na}}.
#' @param k integer indicating order of MA.
#' @noRd
impNA_MA <- function(X, W, k) {
  d <- dim(X)
  TT <- d[1L]
  k2 <- 2L * k + 1L
  weights <- rep(1, k2)/k2
  sss <- -(1:(k2-1L))
  for (i in 1:d[2L]) {
    x <- X[, i]
    nai <- which(W[, i])
    x[nai] <- fmedian.default(x)
    x_MA <- filter(c(rep(x[1L], k), x, rep(x[TT], k)), weights, sides = 1L)
    x[nai] <- x_MA[sss][nai]
    X[, i] <- x
  }
  return(X)
}

#' Impute Time Series with Moving Average and Cubic Spline for Internal Missing Value
#' @inheritParams impNA_MA
#' @noRd
impNA_spline <- function(X, W, k) {
  d <- dim(X)
  TT <- d[1L]
  k2 <- 2L * k + 1L
  weights <- rep(1, k2)/k2
  sss <- -(1:(k2-1L))

  for (i in 1:d[2L]) {
    x <- X[, i]
    nnai <- which(!W[, i])
    ln <- length(nnai)
    t1 <- nnai[1L]
    t2 <- nnai[ln]
    # Cubic spline to interpolate any internal missing values...
    if(ln != t2-t1+1L) x[t1:t2] <- spline(nnai, x[nnai], xout = t1:t2)$y
    isnanx <- which(is.na(x))
    x[isnanx] <- fmedian.default(x)
    x_MA <- filter(c(rep(x[1L], k), x, rep(x[TT], k)), weights, sides = 1L)
    x[isnanx] <- x_MA[sss][isnanx]
    X[, i] <- x
  }
  return(X)
}

#' Remove and Impute Missing Values in a Multivariate Time Series
#'
#' This function imputes missing values in a stationary multivariate time series using various
#' methods, and removes cases with too many missing values.
#'
#' @param X a \code{T x n} numeric data matrix (incl. ts or xts objects) or data frame of stationary time series.
#' @param max.missing numeric. Proportion of series missing for a case to be considered missing.
#' @param na.rm.method character. Method to apply concerning missing cases selected through \code{max.missing}: \code{"LE"} only removes cases at the beginning or end of the sample, whereas \code{"all"} always removes missing cases.
#' @param na.impute character. Method to impute missing values for the PCA estimates used to initialize the EM algorithm. Note that data are standardized (scaled and centered) beforehand. Available options are:
#'    \tabular{llll}{
#' \code{"median"} \tab\tab simple series-wise median imputation. \cr\cr
#' \code{"rnorm"} \tab\tab imputation with random numbers drawn from a standard normal distribution. \cr\cr
#' \code{"median.ma"} \tab\tab values are initially imputed with the median, but then a moving average is applied to smooth the estimates. \cr\cr
#' \code{"median.ma.spline"} \tab\tab "internal" missing values (not at the beginning or end of the sample) are imputed using a cubic spline, whereas missing values at the beginning and end are imputed with the median of the series and smoothed with a moving average.\cr\cr
#' }
#' @param ma.terms the order of the (2-sided) moving average applied in \code{na.impute} methods \code{"median.ma"} and \code{"median.ma.spline"}.
#'
#' @returns The imputed matrix \code{X_imp}, with attributes:
#'  \item{\code{"missing"}}{a missingness matrix \code{W} matching the dimensions of \code{X_imp}. }
#'  \item{\code{"rm.rows"}}{and a vector of indices of rows (cases) with too many missing values that were removed. }
#'
#'
#' @examples
#' library(xts)
#' str(tsnarmimp(BM14_M))
#'
#' @export
tsnarmimp <- function(X,
                       max.missing = 0.8,
                       na.rm.method = c("LE", "all"),
                       na.impute = c("median.ma.spline", "median.ma", "median", "rnorm"),
                       ma.terms = 3L) {

  if(length(max.missing) + length(ma.terms) != 2L) stop("Parameters max.missing and ma.terms must be length 1")
  if(!is.numeric(max.missing) || max.missing < 0 || max.missing > 1) stop("max.missing must be a proportion between 0 and 1")
  if(!is.numeric(ma.terms) || ma.terms <= 0) stop("ma.terms needs to be a positive integer")
  if(!is.integer(ma.terms)) ma.terms <- as.integer(ma.terms)

  W <- is.na(X)
  n <- dim(X)[2L]
  rm.rows <- NULL
  if(max.missing < 1) {
    thresh <- switch(tolower(na.rm.method[1L]),
                     le = findNA_LE(W, max.missing),
                     all = rowSums(W) > max.missing * n,
                     stop("Unknown na.rm.method:", na.rm.method[1L]))
    if(any(thresh)) {
      rm.rows <- which(thresh)
      X <- X[-rm.rows, ]
      W <- W[-rm.rows, ]
    }
  }
  X_imp <- switch(tolower(na.impute[1L]),
                  median = fmedian(X, TRA = "replace_NA"), # replace(X, W, fmedian(X, TRA = 1L)[W]),
                  rnorm = replace(X, W, rnorm(sum(W))),
                  median.ma = impNA_MA(X, W, ma.terms),
                  median.ma.spline = impNA_spline(X, W, ma.terms),
                  stop("Unknown na.impute option:", na.impute[1L]))

   attr(X_imp, "missing") <- W
   attr(X_imp, "rm.rows") <- rm.rows
   return(X_imp)
}

