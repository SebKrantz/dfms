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
