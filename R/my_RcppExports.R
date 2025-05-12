Estep <- function(X, A, C, Q, R, F_0, P_0) {
  .Call(Cpp_Estep, X, A, C, Q, R, F_0, P_0)
}

#' (Fast) Stationary Kalman Filter
#'
#' @description A simple and fast C++ implementation of the Kalman Filter for stationary data (or random walks - data should be mean zero and without a trend) with time-invariant system matrices and missing data.
#' @param X numeric data matrix (\eqn{T \times n}{T x n}).
#' @param A transition matrix (\eqn{rp \times rp}{rp x rp}).
#' @param C observation matrix (\eqn{n \times rp}{n x rp}).
#' @param Q state covariance (\eqn{rp \times rp}{rp x rp}).
#' @param R observation covariance (\eqn{n \times n}{n x n}).
#' @param F_0 initial state vector (\eqn{rp \times 1}{rp x 1}).
#' @param P_0 initial state covariance (\eqn{rp \times rp}{rp x rp}).
#' @param loglik logical. Compute log-likelihood?
#'
#' @details The underlying state space model is:
#'
#' \deqn{\textbf{x}_t = \textbf{C} \textbf{F}_t + \textbf{e}_t \sim N(\textbf{0}, \textbf{R})}{x(t) = C F(t) + e(t) ~ N(0, R)}
#' \deqn{\textbf{F}_t = \textbf{A F}_{t-1} + \textbf{u}_t \sim N(\textbf{0}, \textbf{Q})}{F(t) = A F(t-1) + u(t) ~ N(0, Q)}
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
#' For another fast (C-based) implementation that also allows time-varying system matrices and non-stationary data see \code{FKF::fkf}.
#'
#' @references
#' Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and Its Applications: With R Examples. Springer.
#'
#' Harvey, A. C. (1990). Forecasting, structural time series models and the Kalman filter.
#'
#' Hamilton, J. D. (1994). Time Series Analysis. Princeton university press.
#'
#' @returns Predicted and filtered state vectors and covariances.
#' \item{\code{F}}{\eqn{T \times rp}{T x rp} filtered state vectors. }
#' \item{\code{P}}{\eqn{rp \times rp \times T}{rp x rp x T} filtered state covariances. }
#' \item{\code{F_pred}}{\eqn{T \times rp}{T x rp} predicted state vectors. }
#' \item{\code{P_pred}}{\eqn{rp \times rp \times T}{rp x rp x T} predicted state covariances. }
#' \item{\code{loglik}}{value of the log likelihood. }
#'
#' @seealso \code{\link{FIS}} \code{\link{SKFS}} \link{dfms-package}
#' @examples # See ?SKFS
#'
#' @export
SKF <- function(X, A, C, Q, R, F_0, P_0, loglik = FALSE) {
  .Call(Cpp_SKF, X, A, C, Q, R, F_0, P_0, loglik)
}

#' (Fast) Fixed-Interval Smoother (Kalman Smoother)
#' @param A transition matrix (\eqn{rp \times rp}{rp x rp}).
#' @param F state estimates (\eqn{T \times rp}{T x rp}).
#' @param F_pred state predicted estimates (\eqn{T \times rp}{T x rp}).
#' @param P variance estimates (\eqn{rp \times rp \times T}{rp x rp x T}).
#' @param P_pred predicted variance estimates (\eqn{rp \times rp \times T}{rp x rp x T}).
#' @param F_0 initial state vector (\eqn{rp \times 1}{rp x 1}) or empty (\code{NULL}).
#' @param P_0 initial state covariance (\eqn{rp \times rp}{rp x rp}) or empty (\code{NULL}).
#'
#' @details The Kalman Smoother is given by:
#'
#' \deqn{\textbf{J}_t = \textbf{P}_t \textbf{A} + inv(\textbf{P}^{pred}_{t+1})}{J(t) = P(t) A inv(P_pred(t+1))}
#' \deqn{\textbf{F}^{smooth}_t = \textbf{F}_t + \textbf{J}_t (\textbf{F}^{smooth}_{t+1} - \textbf{F}^{pred}_{t+1})}{F_smooth(t) = F(t) + J(t) (F_smooth(t+1) - F_pred(t+1))}
#' \deqn{\textbf{P}^{smooth}_t = \textbf{P}_t + \textbf{J}_t (\textbf{P}^{smooth}_{t+1} - \textbf{P}^{pred}_{t+1}) \textbf{J}_t'}{P_smooth(t) = P(t) + J(t) (P_smooth(t+1) - P_pred(t+1)) J(t)'}
#'
#' The initial smoothed values for period t = T are set equal to the filtered values. If \code{F_0} and \code{P_0} are supplied, the smoothed initial conditions (t = 0 values) are also calculated and returned.
#' For further details see any textbook on time series such as Shumway & Stoffer (2017), which provide an analogous R implementation in \code{astsa::Ksmooth0}.
#'
#'
#' @returns Smoothed state and covariance estimates, including initial (t = 0) values.
#' \item{\code{F_smooth}}{\eqn{T \times rp}{T x rp} smoothed state vectors, equal to the filtered state in period \eqn{T}.}
#' \item{\code{P_smooth}}{\eqn{rp \times rp \times T}{rp x rp x T} smoothed state covariance, equal to the filtered covariance in period \eqn{T}.}
#' \item{\code{F_smooth_0}}{\eqn{1 \times rp}{1 x rp} initial smoothed state vectors, based on \code{F_0}.}
#' \item{\code{P_smooth_0}}{\eqn{rp \times rp}{rp x rp} initial smoothed state covariance, based on \code{P_0}.}
#'
#' @references
#' Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and Its Applications: With R Examples. Springer.
#'
#' Harvey, A. C. (1990). Forecasting, structural time series models and the Kalman filter.
#'
#' @seealso \code{\link{SKF}} \code{\link{SKFS}} \link{dfms-package}
#' @examples # See ?SKFS
#'
#' @export
FIS <- function(A, F, F_pred, P, P_pred, F_0 = NULL, P_0 = NULL) {
  .Call(Cpp_FIS, A, F, F_pred, P, P_pred, F_0, P_0)
}

# @param loglik integer. 0 does not compute the likelihood, 1 computes a standard Kalman Filter likelihood, 2 computes the likelihood for Banbura and Modugno (2014).
#' (Fast) Stationary Kalman Filter and Smoother
#' @inheritParams SKF
#'
#' @returns All results from \code{\link{SKF}} and \code{\link{FIS}}, and additionally
#' a \eqn{rp \times rp \times T}{rp x rp x T} matrix \code{PPm_smooth}, which is equal to the estimate of \eqn{Cov(F^{smooth}_t, F^{smooth}_{t-1} | T)}{Cov(F_smooth(t), F_smooth(t-1) | T)} and needed for EM iterations.
#' See 'Property 6.3: The Lag-One Covariance Smoother' in Shumway & Stoffer (2017).
#'
#'
#' @seealso \code{\link{SKF}} \code{\link{FIS}} \link{dfms-package}
#'
#' @references
#' Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and Its Applications: With R Examples. Springer.
#'
#' @examples
#' library(collapse)
#'
#' ## Two-Step factor estimates from monthly BM (2014) data
#' X <- fscale(diff(qM(BM14_M))) # Standardizing as KF has no intercept
#' r <- 5L # 5 Factors
#' p <- 3L # 3 Lags
#' n <- ncol(X)
#'
#' ## Initializing the Kalman Filter with PCA results
#' X_imp <- tsnarmimp(X)                 # Imputing Data
#' v <- eigen(cov(X_imp))$vectors[, 1:r] # PCA
#' F_pc <- X_imp %*% v                   # Principal component factor estimates
#' C <- cbind(v, matrix(0, n, r*p-r))    # Observation matrix
#' res <- X - tcrossprod(F_pc, v)        # Residuals from static predictions
#' R <- diag(fvar(res))                  # Observation residual covariance
#' var <- .VAR(F_pc, p)                  # VAR(p)
#' A <- rbind(t(var$A), diag(1, r*p-r, r*p))
#' Q <- matrix(0, r*p, r*p)              # VAR residual matrix
#' Q[1:r, 1:r] <- cov(var$res)
#' F_0 <- var$X[1L, ]                    # Initial factor estimate and covariance
#' P_0 <- ainv(diag((r*p)^2) - kronecker(A,A)) %*% unattrib(Q)
#' dim(P_0) <- c(r*p, r*p)
#'
#' ## Run standartized data through Kalman Filter and Smoother once
#' kfs_res <- SKFS(X, A, C, Q, R, F_0, P_0, FALSE)
#'
#' ## Two-step solution is state mean from the Kalman Smoother
#' F_kal <- kfs_res$F_smooth[, 1:r, drop = FALSE]
#' colnames(F_kal) <- paste0("f", 1:r)
#'
#' ## See that this is equal to the Two-Step estimate by DFM()
#' all.equal(F_kal, DFM(X, r, p, em.method = "none", pos.corr = FALSE)$F_2s)
#'
#' ## Same in two steps using SKF() and FIS()
#' kfs_res2 <- with(SKF(X, A, C, Q, R, F_0, P_0, FALSE), FIS(A, F, F_pred, P, P_pred))
#' F_kal2 <- kfs_res2$F_smooth[, 1:r, drop = FALSE]
#' colnames(F_kal2) <- paste0("f", 1:r)
#' all.equal(F_kal, F_kal2)
#'
#' rm(X, r, p, n, X_imp, v, F_pc, C, res, R, var, A, Q, F_0, P_0, kfs_res, F_kal, kfs_res2, F_kal2)
#'
#' @export
SKFS <- function(X, A, C, Q, R, F_0, P_0, loglik = FALSE) {
  .Call(Cpp_SKFS, X, A, C, Q, R, F_0, P_0, loglik)
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
#' @examples
#' ainv(crossprod(diff(EuStockMarkets)))
#'
#' @returns The matrix-inverse or pseudo-inverse.
#'
#' @seealso \link{dfms-package}
#'
#' @export
ainv <- function(x) .Call(Cpp_ainv, x)
# {
#   dn <- dimnames(x)
#   if(is.null(dn)) return(.Call(Cpp_ainv, x))
#   `dimnames<-`(.Call(Cpp_ainv, x), dn)
# }


#' @rdname ainv
#' @export
apinv <- function(x) .Call(Cpp_apinv, x)
# {
#   dn <- dimnames(x)
#   if(is.null(dn)) return(.Call(Cpp_apinv, x))
#   `dimnames<-`(.Call(Cpp_apinv, x), dn)
# }
