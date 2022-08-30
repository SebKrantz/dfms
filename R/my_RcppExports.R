Estep <- function(X, A, C, Q, R, F_0, P_0) {
  .Call(Cpp_Estep, X, A, C, Q, R, F_0, P_0)
}

#' Fast Kalman Filter
#'
#' @description A simple and fast C++ implementation of the Kalman Filter for stationary data with time-invariant system matrices and missing data.
#' @param X Data matrix (T x n)
#' @param A Transition matrix (rp x rp)
#' @param C Observation matrix (n x rp)
#' @param Q State covariance (rp x rp)
#' @param R Observation covariance (n x n)
#' @param F_0 Initial state vector (rp x 1)
#' @param P_0 Initial state covariance (rp x rp)
#' @param loglik logical. Compute log-likelihood?
#'
#' @details The underlying state space model is:
#'
#' \deqn{\textbf{x}_t = \textbf{C} \textbf{F}_t + \textbf{e}_t \tilde N(\textbf{0}, \textbf{R})}{x(t) = C F(t) + e(t) ~ N(0, R)}
#' \deqn{\textbf{F}_t = \textbf{A F}_{t-1} + \textbf{u}_t \tilde N(0, \textbf{Q})}{F(t) = A F(t-1) + u(t) ~ N(0, Q)}
#'
#' where \eqn{x_t}{x(t)} is \code{X[t, ]}. The filter then first performs a time update (prediction)
#'
#' \deqn{\textbf{F}_t = \textbf{A F}_{t-1}}{F(t) = A F(t-1)}
#' \deqn{\textbf{P}_t = \textbf{A P}_{t-1}\textbf{A}' + \textbf{Q}}{P(t) = A P(t-1) A' + Q}
#'
#' where \eqn{P_t = Cov(F_t)}{P(t) = Cov(F(t))}. This is followed by the measurement update (filtering)
#'
#' \deqn{\textbf{K}_t = \textbf{P}_t \textbf{C}' (\textbf{C P}_t \textbf{C}' + \textbf{R})^{-1}}{K(t) = P(t) C' inv(C P(t) C' + R)}
#' \deqn{\textbf{F}_t = \textbf{F}_t + \textbf{K}_t (\textbf{x}_t - \textbf{C F}_t)}{F(t) = F(t) + K(t) (x(t) - C F(t))}
#' \deqn{\textbf{P}_t = \textbf{P}_t - \textbf{K}_t\textbf{C P}_t}{P(t) = P(t) - K(t) C P(t)}
#'
#' If a row of the data is all missing the measurement update is skipped i.e. the prediction becomes the filtered value. The log-likelihood is
#' computed as
#' \deqn{1/2 \sum_t \log(|St|)-e_t'S_te_t-n\log(2\pi)}{1/2 sum_t[log(det(S(t))) - e(t)' S(t) e(t) - n log(2 pi)]}
#' where \eqn{S_t = (C P_t C' + R)^{-1}}{S(t) = inv(C P(t) C' + R)} and \eqn{e_t = x_t - C F_t}{e(t) = x(t) - C F(t)} is the prediction error.
#'
#' For further details see any textbook on time series such as Shumway & Stoffer (2017), which provide an analogous R implementation in \code{astsa::Kfilter0}.
#' For another fast (C-based) implementation that allows time-varying system matrices and non-stationary data see \code{FKF::fkf}.
#'
#' @references
#' Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and Its Applications: With R Examples. Springer.
#'
#' Harvey, A. C. (1990). Forecasting, structural time series models and the Kalman filter.
#'
#' Hamilton, J. D. (1994). Time Series Analysis. Princeton university press.
#'
#' @returns Predicted and filtered state vectors and covariances.
#' \tabular{lll}{
#' F \tab\tab T x rp filtered state vectors \cr\cr
#' P \tab\tab rp x rp x T filtered state covariances \cr\cr
#' F_pred \tab\tab T x rp predicted state vectors \cr\cr
#' P_pred \tab\tab rp x rp x T predicted state covariances \cr\cr
#' loglik \tab\tab value of the log likelihood
#' }
#' @export
fKF <- function(X, A, C, Q, R, F_0, P_0, loglik = FALSE) {
  .Call(Cpp_fKF, X, A, C, Q, R, F_0, P_0, loglik)
}

#' Fast Kalman Smoother
#' @param A Transition matrix (rp x rp)
#' @param F State estimates (T x rp)
#' @param F_pred State predicted estimates (T x rp)
#' @param P Variance estimates (rp x rp x T)
#' @param P_pred Predicted variance estimates (rp x rp x T)
#' @param F_0 Initial state vector (rp x 1) or empty (NULL)
#' @param P_0 Initial state covariance (rp x rp) or empty (NULL)
#'
#' @details The Kalman Smoother is given by:
#'
#' \deqn{\textbf{J}_t = \textbf{P}_t \textbf{A} + inv(\textbf{P}^{pred}_{t+1})}{J(t) = P(t) A inv(P_pred(t+1))}
#' \deqn{\textbf{F}^{smooth}_t = \textbf{F}_t + \textbf{J}_t (\textbf{F}^{smooth}_{t+1} - \textbf{F}^{pred}_{t+1})}{F_smooth(t) = F(t) + J(t) (F_smooth(t+1) - F_pred(t+1))}
#' \deqn{\textbf{P}^{smooth}_t = \textbf{P}_t + \textbf{J}_t (\textbf{P}^{smooth}_{t+1} - \textbf{P}^{pred}_{t+1}) \textbf{J}_t'}{P_smooth(t) = P(t) + J(t) (P_smooth(t+1) - P_pred(t+1)) J(t)'}
#'
#' The initial smoothed values for period t = T are set equal to the filtered values. If \code{F_0} and \code{P_0} are supplied, the smoothed initial conditions (t = 0 values) are also calculated and returned.
#'
#' @returns Smoothed state and covariance estimates, including initial (t = 0) values.
#' \tabular{lll}{
#' F_smooth \tab\tab T x rp smoothed state vectors, equal to the filtered state in period T \cr\cr
#' P_smooth \tab\tab rp x rp x T smoothed state covariance, equal to the filtered covariance in period T \cr\cr
#' F_smooth_0 \tab\tab 1 x rp initial smoothed state vectors, based on F_0 \cr\cr
#' P_smooth_0 \tab\tab rp x rp initial smoothed state covariance, based on P_0
#' }
#'
#' @references
#' Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and Its Applications: With R Examples. Springer.
#'
#' Harvey, A. C. (1990). Forecasting, structural time series models and the Kalman filter.
#'
#' @export
fKS <- function(A, F, F_pred, P, P_pred, F_0 = NULL, P_0 = NULL) {
  .Call(Cpp_fKS, A, F, F_pred, P, P_pred, F_0, P_0)
}

# @param loglik integer. 0 does not compute the likelihood, 1 computes a standard Kalman Filter likelihood, 2 computes the likelihood for Banbura and Modungo (2014).
#' Fast Kalman Filter and Smoother
#' @inheritParams fKF
#'
#' @returns All results from \code{\link{fKF}} and \code{\link{fKS}}, and additionally
#' a rp x rp x T matrix \code{PPm_smooth}, which is equal to the estimate of Cov(F_smooth_t, F_smooth_t-1|T) and needed for EM iterations.
#' @export
fKFS <- function(X, A, C, Q, R, F_0, P_0, loglik = FALSE) {
  .Call(Cpp_fKFS, X, A, C, Q, R, F_0, P_0, loglik)
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
