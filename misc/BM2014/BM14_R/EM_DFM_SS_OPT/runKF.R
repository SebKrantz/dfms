#--------------------------------------------------------------------------
# KALMAN FILTER
#--------------------------------------------------------------------------
# y = y_est; x_0 = Z_0; Sig_0 = V_0
runKF <- function(y, A, C, Q, R, x_0, Sig_0, S) { # [xsmooth, Vsmooth, VVsmooth, loglik]

  S = SKF(y, C, R, A, Q, x_0, Sig_0, S) # Stationary Kalman Filter
  S = FIS(y, C, R, A, Q, S)             # Fixed-Interval Smoother

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
SKF <- function(Y, Z, R, T, Q, A_0, P_0, S) {

  # Output structure & dimensions
  c("n", "m") %=% dim(Z) # Note: m = rp
  nobs  = ncol(Y)

  # ______________________________________________________________________
  Au = A_0  # A_0|0;
  Pu = P_0  # P_0|0

  S$AmU[, 1L]  = Au
  S$PmU[,, 1L] = Pu
  e = diag(n)     # Used to be in MissData

  for (t in 1:nobs) {
    # print(t)
    # A = A_t|t-1 & P = P_t|t-1
    A   = T %*% Au
    P   = T %*% tcrossprod(Pu, T) %+=% Q
    P   =  0.5 * (P + t(P))

    # handling the missing data
    c("y_t", "Z_t", "R_t", "L_t") %=% MissData(Y[, t], Z, R, e)

      if(!length(y_t)) {
        Au = A
        Pu = P
      } else {
        PZ  = tcrossprod(P, Z_t)
        iF  = cinv(Z_t %*% PZ %+=% R_t)
        PZF = PZ %*% iF
        V   = y_t - Z_t %*% A
        Au  = A  + PZF %*% V
        Pu  = P  - tcrossprod(PZF, PZ)
        Pu  =  0.5 * (Pu + t(Pu))
        S$loglik = S$loglik + 0.5 * (log(det(iF)) - crossprod(V, iF) %*% V)
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

  c("m", "nobs") %=% dim(S$Am)
  nobsp = nobs + 1L
  S$AmT             = matrix(0, m, nobsp)
  S$PmT             = array(0, c(m, m, nobsp))
  S$PmT_1           = array(0, c(m, m, nobs))
  S$AmT[, nobsp]    = S$AmU[, nobsp]
  S$PmT[,, nobsp]   = S$PmU[,, nobsp]
  S$PmT_1[,, nobs]  = (diag(m) - S$KZ) %*% T %*% S$PmU[,, nobs]

  dm = c(m, m)
  tmp = S$Pm[,, nobs]
  dim(tmp) = dm
  J_2 = tcrossprod(S$PmU[,, nobs], T) %*% apinv(tmp)

  for (t in nobs:1) {
    PmU = S$PmU[,, t]
    Pm1 = S$Pm[,, t]
    P_T = S$PmT[,, t + 1L]
    P_T1 = S$PmT_1[,, t]
    J_1 = J_2

    S$AmT[, t]  = S$AmU[, t] + J_1 %*% (S$AmT[, t + 1L] - T %*% S$AmU[, t])
    S$PmT[,, t] = PmU        + J_1 %*% tcrossprod(P_T - Pm1, J_1)

    if(t > 1L) {
      tmp = S$Pm[,, t-1L]
      dim(tmp) = dm
      J_2 = tcrossprod(S$PmU[,, t-1L], T) %*% apinv(tmp)
      S$PmT_1[,, t-1L] = tcrossprod(PmU, J_2) + J_1 %*% tcrossprod(P_T1 - T %*% PmU, J_2)
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

MissData <- function(y, C, R, e) { # [y,C,R,L]
  if(!anyNA(y)) return(list(y, C, R, e))
  ix = whichNA(y, invert = TRUE)
  # if(length(ix) == length(y)) return(list(y, C, R, e))
  list(y[ix],
       C[ix,, drop = FALSE],
       R[ix, ix, drop = FALSE],
       e[, ix, drop = FALSE])
}
