#' Implementation of a Kalman filter
#' @param y Data matrix (T x n)
#' @param H Observation matrix
#' @param Q State covariance
#' @param R Observation covariance
#' @param F Transition matrix
#' @param x0 Initial state vector
#' @param P0 Initial state covariance
#' @export
KalmanFilter <- function(y, H, Q, R, F, x0, P0) {
  .Call(Cpp_KalmanFilter, y, H, Q, R, F, x0, P0)
}

#' Runs a Kalman smoother
#' @param F transition matrix
#' @param H observation matrix
#' @param R Observation covariance
#' @param xfT State estimates
#' @param xpTm State predicted estimates
#' @param PfT_v Variance estimates
#' @param PpT_v Predicted variance estimates
#' @return List of smoothed estimates
#' @export
KalmanSmoother <- function(F, H, R, xfT, xpT, PfT_v, PpT_v) {
  .Call(Cpp_KalmanSmoother, F, H, R, xfT, xpT, PfT_v, PpT_v)
}


KalmanFilterSmoother <- function(y, H, Q, R, F, x0, P0) {
  .Call(Cpp_KalmanFilterSmoother, y, H, Q, R, F, x0, P0)
}

Estep <- function(y, H, Q, R, F, x0, P0) {
  .Call(Cpp_Estep, y, H, Q, R, F, x0, P0)
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
ainv <- function(x) {
  dn <- dimnames(x)
  if(is.null(dn)) return(.Call(Cpp_ainv, x))
  `dimnames<-`(.Call(Cpp_ainv, x), dn)
}

#'
#' @rdname ainv
#' @export
apinv <- function(x) {
  dn <- dimnames(x)
  if(is.null(dn)) return(.Call(Cpp_apinv, x))
  `dimnames<-`(.Call(Cpp_apinv, x), dn)
}
