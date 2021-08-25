#' Estimate a Dynamic Factor Model
#'
#' Estimate a Dynamic Factor Model with arbitrary patterns of missing data
#' on stationary data of a single frequency, with time-invariant system matrices
#' and idiosynchratic measurement (observation) errors.
#'
#' @param X data matrix or frame.
#' @param r number of factors.
#' @param p number of lags in factor VAR.
#' @param max.inter maximum number of EM iterations.
#' @param tol EM convergence tolerance.
#' @param miss.threshold proportions of variables missing for a case to be removed.
#'
#' @details
#'
#' The State-Space Model estimated is:
#'
#' Xt = C Ft + et ~ N(0, R)\cr
#' Ft = A Ft-1 + ut ~ N(0, Q)
#'
#' Where
#' \itemize{
#' \item{n}{Number of variables}
#' \item{T}{Number of observations}
#' \item{F}{Factor vector (stacked, starting with all factors at time t, followed by t-1 etc., rp x 1)}
#' \item{C}{Observation matrix (n x rp, only the first n x r terms are non-zero (no relationship with lagged factors))}
#' \item{A}{State transition matrix (rp x rp, lower parts are known (identity matrix))}
#' \item{Q}{State covariance matrix (rp x rp, top r x r part is contemporaneous, rest is 0 (we assume all relationship between lagged factors is captured in A))}
#' \item{R}{Observation covariance matrix (n x n, diagonal by assumption)}
#' }
#' @useDynLib DFM, .registration = TRUE
#' @export

DFM <- function(X, r, p = 1L,
                rQ = c("none", "diagonal", "identity"),
                rR = c("diagonal", "identity"),
                max.iter = 100, tol = 1e-4,
                miss.threshold = 0.8) {

  rp = r * p
  X = fscale(qM(X))
  T = dim(X)[1L]
  n = dim(X)[2L]

  # Missing values
  X_imp = X
  anymiss = anyNA(X)
  if(anymiss) { # Simple median imputation
    W = is.na(X)
    X_imp[W] = fmedian(X, TRA = 1L)[W]
  }

  # Run PCA to get initial factor estimates:
  v = svd(X_imp, nu = 0L, nv = min(as.integer(r), n, T))$v
  F_ini = X_imp %*% v

  # Observation equation -------------------------------
  # Static predictions (all.equal(unattrib(HDB(X_imp, F_ini)), unattrib(F_ini %*% t(v))))
  C_ini <- cbind(v, matrix(0, n, rp - r))
  res = X_imp - F_ini %*% t(v) # residuals from static predictions
  if(anymiss) res[W] = NA
  R_ini = diag(fvar(res)) # Covariance (assumed idiosynchratic)
  # R_ini = cov(res) # unrestricted covariance estimate

  # Transition equation -------------------------------
  var = VAR(F_ini, p)
  A_ini = rbind(t(var$A), diag(1, rp-r, rp)) # var$A is rp x r matrix
  Q_ini = matrix(0, rp, rp)
  Q_ini[1:r, 1:r] = cov(var$res) # unrestricted covariance estimate

  # Initial state and state covariance (P) ------------
  x0 = var$X[1L, ] # rep(0, rp) # This should better be called f0, the factors are the state
  # Kalman gain is normally A %*% t(A) + Q, but here A is somewhat tricky...
  P0 = matrix(ginv(kronecker(A_ini, A_ini)) %*% unattrib(Q_ini), rp, rp)
  # BM2014:
  # P0 = matrix(solve(diag(rp^2) - kronecker(A_ini, A_ini)) %*% unattrib(Q_ini), rp, rp)

  # FKF::fkf(x0, P0, x0, rep(0, n), A_ini, C_ini, Q_ini, R_ini, X)

  ## Run standartized data through Kalman filter and smoother once
  kf_res <- KalmanFilter(X, C_ini, Q_ini, R_ini, A_ini, x0, P0)
  ks_res <- with(kf_res, KalmanSmoother(A_ini, C_ini, R_ini, xF, xP, Pf, Pp))

  ## Two-step solution is state mean from the Kalman smoother
  F_kal <- ks_res$xS
  return(F_kal[, 1:r])

}
