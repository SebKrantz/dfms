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
KalmanFilter <- function(X, A, C, Q, R, F0, P0) {
  .Call(Cpp_KalmanFilter, X, A, C, Q, R, F0, P0)
}

#' Runs a Kalman smoother
#' @param A Transition matrix (rp x rp)
#' @param C Observation matrix (n x rp)
#' @param R Observation covariance (n x n)
#' @param ZTf State estimates
#' @param ZTp State predicted estimates
#' @param VTf_v Variance estimates
#' @param VTp_v Predicted variance estimates
#' @return List of smoothed estimates
KalmanSmoother <- function(A, C, R, ZTf, ZTp, VTf_v, VTp_v) {
  .Call(Cpp_KalmanSmoother, A, C, R, ZTf, ZTp, VTf_v, VTp_v)
}

#' Kalman Filter and Smoother
#' @param X Data matrix (T x n)
#' @param A Transition matrix (rp x rp)
#' @param C Observation matrix (n x rp)
#' @param Q State covariance (rp x rp)
#' @param R Observation covariance (n x n)
#' @param F0 Initial state vector (rp x 1)
#' @param P0 Initial state covariance (rp x rp)
KalmanFilterSmoother <- function(X, A, C, Q, R, F0, P0) {
  .Call(Cpp_KalmanFilterSmoother, X, A, C, Q, R, F0, P0)
}

ainv <- function(x) {
  .Call(Cpp_ainv, x)
}

apinv <- function(x) {
  .Call(Cpp_apinv, x)
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
