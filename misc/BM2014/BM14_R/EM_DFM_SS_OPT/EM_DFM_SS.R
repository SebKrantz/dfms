library(collapse) # 1.8+
library(roll)
library(DFM)

# Note: This is a decently optimized version of the original BM 2014 code

setwd("~/Documents/R/DFM/misc/BM2014/BM14_R/EM_DFM_SS_OPT")
source("remNaNs_spline.R")
source("runKF.R")
source("Procedures.R")

# TODO: make it work with just 1 series?
# Also make sure cinv solution is robust!!
# Also try to use FKF package...
# FKF::fkf

EM_DFM_SS_OPT <- function(X, r, p = 1L, max_iter = 100L, thresh = 1e-4) { # Res

  # --------------------------------------------------------------------------
  # Preparation of the data
  #--------------------------------------------------------------------------

  X = qM(X)
  c("T", "N") %=% dim(X)

  # Standardise x
  Mx = fmean(X)
  Wx = fsd(X)
  xNaN = fscale(X)

  # --------------------------------------------------------------------------
  # Initial Conditions
  #--------------------------------------------------------------------------

  # Removing missing values (for initial estimators)
  optNaN = list()
  optNaN$method = 2 # Remove leading and closing zeros
  optNaN$k = 3      # order of the moving average for replacing the missing observations

  # return(InitCond(xNaN, r, p, optNaN))
  c("A", "C", "Q", "R", "Z_0", "V_0") %=% InitCond(xNaN, r, p, optNaN)

  # some auxiliary variables for the iterations
  previous_loglik = -Inf
  num_iter = 0L
  LL = -Inf
  converged = FALSE

  # y for the estimation is WITH missing data
  y = t(xNaN)

  #--------------------------------------------------------------------------
  # THE EM LOOP
  #--------------------------------------------------------------------------

  # The model can be written as
  # y = C*Z + e
  # Z = A*Z(-1) + v
  # where y is NxT, Z is (pr)xT, etc.

  # remove the leading and ending nans for the estimation
  optNaN$method = 3
  y_est = t(remNaNs(xNaN, optNaN)$X)
  dnkron = matrix(1, r, r) %x% diag(nrow(y_est)) # Used to be inside EMstep
  dnkron_ind = whichv(dnkron, 1)
  nobs = ncol(y_est)
  rp = r*p
  # Matrices for Kalman Filter, used to be in SKF
  S = list(Am  = matrix(NA_real_, rp, nobs),
            Pm  = array(NA_real_, c(rp, rp, nobs)),
            AmU = matrix(NA_real_, rp, nobs + 1L),
            PmU = array(NA_real_, c(rp, rp, nobs + 1L)),
            loglik = 0)

  # return(EMstep_OPT(y_est, A, C, Q, R, Z_0, V_0, r, dnkron, dnkron_ind, S))

  while(num_iter < max_iter && !converged) {

    c("C_new", "R", "A", "Q", "Z_0", "V_0", "loglik") %=% EMstep_OPT(y_est, A, C, Q, R, Z_0, V_0, r, dnkron, dnkron_ind, S)
    C[, 1:r] = C_new # C = C_new

    # Checking convergence
    c("converged", "decrease") %=% em_converged(loglik, previous_loglik, thresh, TRUE)

    LL = c(LL, loglik)
    previous_loglik = loglik
    num_iter =  num_iter + 1L

  }

  # final run of the Kalman filter
  Zsmooth = t(runKF(y, A, C, Q, R, Z_0, V_0, S)$xsmooth)
  F = Zsmooth[-1L, 1:r, drop = FALSE]
  Res = list()
  Res$x_sm = tcrossprod(Zsmooth[-1L, , drop = FALSE], C)
  Res$X_sm = setop(Res$x_sm %r*% Wx, "+", Mx, rowwise = TRUE)
  Res$F = F

  #--------------------------------------------------------------------------
  #   Loading the structure with the results
  #--------------------------------------------------------------------------
  Res$C = C
  Res$R = R
  Res$A = A
  Res$Q = Q
  Res$Z_0 = Z_0
  Res$V_0 = V_0
  Res$r = r
  Res$p = p
  Res$Mx = Mx
  Res$Wx = Wx

  return(Res)
}




