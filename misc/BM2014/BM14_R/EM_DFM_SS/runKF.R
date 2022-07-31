#--------------------------------------------------------------------------
# KALMAN FILTER
#--------------------------------------------------------------------------
# y = y_est; x_0 = Z_0; Sig_0 = V_0
runKF <- function(y, A, C, Q, R, x_0, Sig_0) { # [xsmooth, Vsmooth, VVsmooth, loglik]

  S = SKF(y, C, R, A, Q, x_0, Sig_0) # Stationary Kalman Filter
  S = FIS(y, C, R, A, Q, S)          # Fixed-Interval Smoother

  list(xsmooth = S$AmT,
       Vsmooth = S$PmT,
       VVsmooth = S$PmT_1,
       loglik = S$loglik)
}


#______________________________________________________________________
# Kalman filter for stationary systems with time-varying system matrices
# and missing data.
#
# The model is        y_t   = Z * a_t + eps_t
#                     a_t+1 = T * a_t + u_t
#
#______________________________________________________________________
# INPUT
#        Y         Data                                  (nobs x n)
# OUTPUT
#        S.Am       Predicted state vector  A_t|t-1      (nobs x m)
#        S.AmU      Filtered  state vector  A_t|t        (nobs+1 x m)
#        S.Pm       Predicted covariance of A_t|t-1      (nobs x m x m)
#        S.PmU      Filtered  covariance of A_t|t        (nobs+1 x m x m)
#        S.loglik   Value of likelihood function
SKF <- function(Y, Z, R, T, Q, A_0, P_0) {

  # Output structure & dimensions
  c("n", "m") %=% dim(Z)
  nobs  = ncol(Y)
  S <- list(Am  = matrix(NA_real_, m, nobs),
            Pm  = array(NA_real_, c(m, m, nobs)),
            AmU = matrix(NA_real_, m, nobs + 1L),
            PmU = array(NA_real_, c(m, m, nobs + 1L)),
            loglik = 0)

  # ______________________________________________________________________
  Au = A_0  # A_0|0;
  Pu = P_0  # P_0|0

  S$AmU[, 1L]  = Au
  S$PmU[,, 1L] = Pu

  for (t in 1:nobs) {
    # print(t)
    # A = A_t|t-1 & P = P_t|t-1
    A   = T %*% Au
    P   = T %*% tcrossprod(Pu, T) + Q
    P   =  0.5 * (P + t(P))

    # handling the missing data
    c("y_t", "Z_t", "R_t", "L_t") %=% MissData(Y[, t], Z, R)

      if(!length(y_t)) {
        Au = A
        Pu = P
      } else {
        PZ  = tcrossprod(P, Z_t)
        iF  = ainv(Z_t %*% PZ + R_t) # solve(, tol = .Machine$double.eps/1e10)
        PZF = PZ %*% iF
        V   = y_t - Z_t %*% A
        Au  = A  + PZF %*% V
        Pu  = P  - tcrossprod(PZF, PZ)
        Pu  =  0.5 * (Pu + t(Pu))
        S$loglik = S$loglik + 0.5 * (log(det(iF))  - crossprod(V, iF) %*% V)
      }

    S$Am[, t]  = A
    S$Pm[,, t] = P

    # Au = A_t|t & Pu = P_t|t
    S$AmU[, t + 1L]  = Au
    S$PmU[,, t + 1L] = Pu
  }

  S$KZ = if(!length(y_t)) matrix(0, m, m) else PZF %*% Z_t
  return(S)
}


#______________________________________________________________________
# Fixed intervall smoother (see Harvey, 1989, p. 154)
# FIS returns the smoothed state vector AmT and its covar matrix PmT
# Use this in conjnuction with function SKF
#______________________________________________________________________
# INPUT
#        Y         Data                                 (nobs x n)
#        S Estimates from Kalman filter SKF
#          S.Am   : Estimates     a_t|t-1                  (nobs x m)
#          S.Pm   : P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
#          S.AmU  : Estimates     a_t|t                    (nobs x m)
#          S.PmU  : P_t|t   = Cov(a_t|t)               (nobs x m x m)
# OUTPUT
#        S Smoothed estimates added to above
#          S.AmT  : Estimates     a_t|T                    (nobs x m)
#          S.PmT :  P_t|T   = Cov(a_t|T)               (nobs x m x m)
#          S.PmU : Cov(a_ta_t-1|T)
#        where m is the dim of state vector and t = 1 ...T is time

# Y = y; Z = C; T = A
FIS <- function(Y, Z, R, T, Q, S) {

  c("m", "nobs")  %=% dim(S$Am)
  nobsp = nobs + 1L
  S$AmT             = matrix(0, m, nobsp)
  S$PmT             = array(0, c(m, m, nobsp))
  S$PmT_1           = array(0, c(m, m, nobs))
  S$AmT[, nobsp]    = drop(S$AmU[, nobsp])
  S$PmT[,, nobsp]   = drop(S$PmU[,, nobsp])
  S$PmT_1[,, nobs]  = (diag(m) - S$KZ) %*% T %*% drop(S$PmU[,, nobs])

  J_2 = drop(S$PmU[,, nobs]) %*% t(T) %*% apinv(drop(S$Pm[,, nobs])) # pracma::pinv

  for (t in nobs:1) {
    PmU = drop(S$PmU[,, t])
    Pm1 = drop(S$Pm[,, t])
    P_T = drop(S$PmT[,, t + 1L])
    P_T1 = drop(S$PmT_1[,, t])

    J_1 = J_2

    S$AmT[, t]  = S$AmU[, t] + J_1 %*% (S$AmT[, t + 1L] - T %*% S$AmU[, t])
    S$PmT[,, t] = PmU        + J_1 %*% (P_T - Pm1) %*% t(J_1)

    if(t > 1L) {
      J_2 = drop(S$PmU[,, t-1L]) %*% t(T) %*% apinv(drop(S$Pm[,, t-1L]))
      S$PmT_1[,, t-1L] = tcrossprod(PmU, J_2) + J_1 %*% (P_T1 - T %*% PmU) %*% t(J_2)
    }
  }
  return(S)
}

#______________________________________________________________________
# PROC missdata
# PURPOSE: eliminates the rows in y & matrices Z, G that correspond to
#          missing data (NaN) in y
# INPUT    y             vector of observations at time t  (n x 1 )
#          S             KF system matrices             (structure)
#                        must contain Z & G
# OUTPUT   y             vector of observations (reduced)   (# x 1)
#          Z G           KF system matrices     (reduced)   (# x ?)
#          L             To restore standard dimensions     (n x #)
#                        where # is the nr of available data in y
#______________________________________________________________________

MissData <- function(y, C, R) { # [y,C,R,L]
  ix = which(!is.na(y))
  e  = diag(length(y))

  list(y  =  y[ix],
       C  =  C[ix, ],
       R  =  R[ix, ix],
       L  = e[, ix])
}
