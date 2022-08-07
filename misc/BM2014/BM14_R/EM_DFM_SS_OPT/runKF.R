#--------------------------------------------------------------------------
# KALMAN FILTER
#--------------------------------------------------------------------------

runKF <- function(Y, A, C, Q, R, Z_0, V_0, S) {

  S = SKF(Y, C, R, A, Q, Z_0, V_0, S)   # Stationary Kalman Filter
  return(FIS(Y, C, R, A, Q, S))         # Fixed-Interval Smoother

}


#______________________________________________________________________
# Kalman filter for stationary systems with time-varying system matrices and missing data.
#
# The model is        Y_t   = C * Z_t + r_t ~ N(0, R)
#                     Z_t+1 = A * Z_t + q_t ~ N(0, Q)
#                     Cov(Z_t) = V; Cov(Z_t+1,Z_t) = VV
#
#______________________________________________________________________
# INPUT
#        Y         Data                                  (n x T)
#        C         Observation Matrix                    (n x rp)
#        R         Observation Covariance Matrix         (n x n)
#        A         Transition Matrix                     (rp x rp)
#        Q         Transition Covariance Matrix          (rp x rp)
#        Z_0       Factor Initialization (zeros)         (rp x 1)
#        V_0       Factor Covariance, PCA Estimates      (rp x rp)
#        S         Structure with empty results matrices
#
# OUTPUT (including initisation values in first obs)
#        S.ZT       Predicted state vector  Z_t|t-1      (rp x T)
#        S.ZT_0     Filtered  state vector  Z_t|t        (rp x T+1)
#        S.VT       Predicted covariance of Z_t|t-1      (rp x rp x T)
#        S.VT_0     Filtered  covariance of Z_t|t        (rp x rp x T+1)
#        S.loglik   Value of likelihood function
#        S.KC       Kalman Gain * observation matrix of final period, for smoothing...

SKF <- function(Y, C, R, A, Q, Z_0, V_0, S) {

  # Output structure & dimensions
  T  = ncol(Y)
  S$ZT_0[, 1L]  = Z_0
  S$VT_0[,, 1L] = V_0

  for (t in 1:T) {

    ### (1) Time update (Predict future state / factors before new measurement comes in)
    Z   = A %*% Z_0                         # A-priori prediction of factors from transition equation: Z = Z_t|t-1
    V   = A %*% tcrossprod(V_0, A) + Q      # A-priori estimate of the transition error covariance:    V = V_t|t-1
    V   =  0.5 * (V + t(V))                 # Ensure symmetry

    ### (2) handling the missing data: simply removing...
    c("y_t", "C_t", "R_t") %=% MissData(Y[, t], C, R)

    ### (3) Measurement Update (Correct)
    if(!length(y_t)) { # If no data, take predictions and move on
      Z_0 = Z
      V_0 = V
    } else {
      # Kalman Gain K = VC'(CVC' + R)^(-1)
      # is chosen to be the gain or blending factor that minimizes the a-posteriori error covariance V.
      # As the measurement error covariance R_t approaches zero, the gain K is larger and weights the measurement residual more heavily
      # If measurement error tends to zero lim R -> 0 then K -> C^(-1).
      VC  = tcrossprod(V, C_t)
      iF  = ainv(C_t %*% VC %+=% R_t)    # Inverse: needed for Kalman Gain and likelihood
      K = VC %*% iF                      # Kalman Gain
      MR  = y_t - C_t %*% Z              # Measurement Residual: also needed for likelihood
      Z_0  = Z + K %*% MR                # A-posteriori factor estimate: Kalman gain weights contribution of new measurement to prediction
      V_0  = V  - tcrossprod(K, VC)      # This now gets the a-posteriori error covariance estimate...
      # ... if measurement error tends to zero lim R -> 0 then K -> C^(-1), and V_0 = V - V = 0 (perfect measurement, zero error)
      V_0  =  0.5 * (V_0 + t(V_0))       # Ensure symmetry
      # Compute log-likelihood
      S$loglik = S$loglik + 0.5 * (log(det(iF)) - crossprod(MR, iF) %*% MR)
    }

    # Save results: Initial prediction
    S$ZT[, t]  = Z  # Z = Z_t|t-1
    S$VT[,, t] = V  # V = V_t|t-1

    # Save results: Final prediction
    S$ZT_0[, t + 1L]  = Z_0  # Z_0 = Z_t|t
    S$VT_0[,, t + 1L] = V_0  # V_0 = V_t|t
  }

  S$KC = if(!length(y_t)) matrix(0, ncol(C), ncol(C)) else K %*% C_t # Needed for smoothing...
  return(S)
}


#______________________________________________________________________
# Fixed intervall smoother (see Harvey, 1989, p. 154)
# FIS returns the smoothed state vector Zsmooth and its covar matrix Vsmooth
# Use this in conjnuction with function SKF
#______________________________________________________________________
# INPUT
#        Y, C, R, A, Q                See SKF Inputs
#        S.ZT, S.ZT_0, S.VT, S.VT_0   See SKF Outputs
#
# OUTPUT
#        S Smoothed estimates added to above
#        S.Zsmooth  : Estimates     Z_t|T                    (T x rp)
#        S.Vsmooth :  V_t|A   = Cov(Z_t|T)               (T x rp x rp)
#        S.VT_0 : Cov(Z_tZ_t-1|T)
#        where rp is the dim of state vector and t = 1 ...T is time

FIS <- function(Y, C, R, A, Q, S) {

  # Output structure & dimensions
  list2env(S, envir = environment())
  c("rp", "T") %=% dim(ZT)
  Tp = T + 1L

  # First Smoothing observation is last filtering observation
  ZsmoothT[, Tp]   = ZT_0[, Tp]
  VsmoothT[,, Tp]  = VT_0[,, Tp]
  VVsmoothT[,, T]  = (diag(rp) - KC) %*% A %*% VT_0[,, T]

  # If rp = 1, need to set dimensions...
  dm = c(rp, rp)
  tmp = VT[,, T]
  dim(tmp) = dm
  J_2 = tcrossprod(VT_0[,, T], A) %*% apinv(tmp)

  # TODO: Add comments + can optimize??
  for (t in T:1) {
    Z_0 = ZT_0[, t] # same as t-1: because has one more slot...
    V_0 = VT_0[,, t]
    J_1 = J_2

    ZsmoothT[, t]  = Z_0 + J_1 %*% (ZsmoothT[, t+1L] - A %*% Z_0)
    VsmoothT[,, t] = V_0 + J_1 %*% tcrossprod(VsmoothT[,, t+1L] - VT[,, t], J_1)

    if(t > 1L) {
      tmp = VT[,, t-1L]
      dim(tmp) = dm
      J_2 = tcrossprod(VT_0[,, t-1L], A) %*% apinv(tmp)
      VVsmoothT[,, t-1L] = tcrossprod(V_0, J_2) + J_1 %*% tcrossprod(VVsmoothT[,, t] - A %*% V_0, J_2)
    }
  }

  return(list(Zsmooth = ZsmoothT,
              Vsmooth = VsmoothT,
              VVsmooth = VVsmoothT,
              loglik = loglik))
}

#______________________________________________________________________
# PROC missdata
# PURPOSE: eliminates the rows in Y & matrices C, R that correspond to
#          missing data (NaN) in Y
# INPUT    y_t           vector of observations at time t   (n x 1)
#          C, R          KF system matrices: n x rp and n x n
#
# OUTPUT   y_t           vector of observations (reduced)   (# x 1)
#          C R           KF system matrices     (reduced)   (# x ?)
#______________________________________________________________________

MissData <- function(y_t, C, R) {
  if(!anyNA(y_t)) return(list(y_t, C, R))
  ix = whichNA(y_t, invert = TRUE)
  # if(length(ix) == length(y_t)) return(list(y_t, C, R))
  list(y_t[ix],
       C[ix,, drop = FALSE],
       R[ix, ix, drop = FALSE])
}
