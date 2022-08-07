Estep <- function(X, A, C, Q, R, F0, P0) {
  .Call(Cpp_Estep, X, A, C, Q, R, F0, P0)
}

#' Implementation of a Kalman filter
#' @param X Data matrix (T x n)
#' @param A Transition matrix (rp x rp)
#' @param C Observation matrix (n x rp)
#' @param Q State covariance (rp x rp)
#' @param R Observation covariance (n x n)
#' @param F0 Initial state vector (rp x 1)
#' @param P0 Initial state covariance (rp x rp)
#' @param loglik logical. Compute log-likelihood?
#'
#' @returns Predicted and filtered state vectors and covariances, including a prediction for period T+1.
#' \tabular{lll}{
#' F \tab\tab T x rp filtered state vectors \cr\cr
#' P \tab\tab rp x rp x T filtered state covariances \cr\cr
#' F_pred \tab\tab T+1 x rp predicted state vectors, prediction for period t = 1 is F0 \cr\cr
#' P_pred \tab\tab rp x rp x T+1 predicted state covariances, prediction for period t = 1 is P0 \cr\cr
#' loglik \tab\tab value of the log likelihood
#' }
#' @export
KalmanFilter <- function(X, A, C, Q, R, F0, P0, loglik = FALSE) {
  .Call(Cpp_KalmanFilter, X, A, C, Q, R, F0, P0, loglik)
}

#' Runs a Kalman smoother
#' @param A Transition matrix (rp x rp)
#' @param F State estimates (T x rp)
#' @param F_pred State predicted estimates (T x rp) or (T+1 x rp)
#' @param P Variance estimates (rp x rp x T)
#' @param P_pred Predicted variance estimates (rp x rp x T) or (rp x rp x T+1)
#'
#' @returns Smoothed state and covariance estimates, including initial (t = 0) values.
#' \tabular{lll}{
#' F_smooth \tab\tab T x rp smoothed state vectors, equal to the filtered state in period T \cr\cr
#' P_smooth \tab\tab rp x rp x T smoothed state covariance, equal to the filtered covariance in period T \cr\cr
#' F_smooth_0 \tab\tab 1 x rp initial smoothed state vectors, based on F0 \cr\cr
#' P_smooth_0 \tab\tab rp x rp initial smoothed state covariance, based on P0
#' }
#' @export
KalmanSmoother <- function(A, F, F_pred, P, P_pred) {
  .Call(Cpp_KalmanSmoother, A, F, F_pred, P, P_pred)
}

#' Kalman Filter and Smoother
#' @inheritParams KalmanFilter
#' @param loglik integer. 0 does not compute the likelihood, 1 computes a standard Kalman Filter likelihood, 2 computes the likelihood for Banbura and Modungo (2014).
#'
#' @returns All results from \code{\link{KalmanFilter}} and \code{\link{KalmanSmoother}}, and additionally
#' a rp x rp x T matrix \code{PPm_smooth}, which is equal to the estimate of Cov(F_smooth_t, F_smooth_t-1|T) and needed for EM iterations.
#' @export
KalmanFilterSmoother <- function(X, A, C, Q, R, F0, P0, loglik = 0L) {
  .Call(Cpp_KalmanFilterSmoother, X, A, C, Q, R, F0, P0, loglik)
}



#' @title Armadillo's Inverse Functions
#' @name ainv
#' @aliases ainv
#' @aliases apinv
#'
#' @description Matrix inverse and pseudo-inverse by the Armadillo C++ library.
#'
#' @param x a numeric matrix, must be square for \code{ainv}.
#'
#' @returns The matrix-inverse or pseudo-inverse.
#' @export
ainv <- function(x) .Call(Cpp_ainv, x)
# {
#   dn <- dimnames(x)
#   if(is.null(dn)) return(.Call(Cpp_ainv, x))
#   `dimnames<-`(.Call(Cpp_ainv, x), dn)
# }

#'
#' @rdname ainv
#' @export
apinv <- function(x) .Call(Cpp_apinv, x)
# {
#   dn <- dimnames(x)
#   if(is.null(dn)) return(.Call(Cpp_apinv, x))
#   `dimnames<-`(.Call(Cpp_apinv, x), dn)
# }
