library(collapse) # 1.8+
library(roll)
library(DFM)

setwd("~/Documents/R/DFM/misc/BM2014/BM14_R/EM_DFM_SS")
source("remNaNs_spline.R")
source("runKF.R")
source("Procedures.R")

EM_DFM_SS <- function(X, r, p = 1L, max_iter = 100L, thresh = 1e-4, na.method = 2L, ma.k = 3L) {

  X = qM(X) # Matlab works with matrices...

  # thresh = 1e-4
  # r = P$r # number of factors
  # p = P$p # number of lags in the factor VAR
  # max_iter = P$max_iter

  # --------------------------------------------------------------------------
  # Preparation of the data
  #--------------------------------------------------------------------------

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
  optNaN$method = na.method # See remNaNs_spline.R
  optNaN$k = ma.k           # order of the moving average for replacing the missing observations

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

  while(num_iter < max_iter && !converged) {

    c("C_new", "R", "A", "Q", "Z_0", "V_0", "loglik") %=% EMstep(y_est, A, C, Q, R, Z_0, V_0, r)
    C[, 1:r] = C_new # C = C_new

    # Checking convergence
    c("converged", "decrease") %=% em_converged(loglik, previous_loglik, thresh, TRUE)

    LL = c(LL, loglik)
    previous_loglik = loglik
    num_iter =  num_iter + 1L

  }

  # final run of the Kalman filter
  Zsmooth = t(runKF(y, A, C, Q, R, Z_0, V_0)$xsmooth)
  F = Zsmooth[-1L, 1:r]
  Res = list()
  Res$x_sm = tcrossprod(Zsmooth[-1L, ], C)
  Res$X_sm = (Res$x_sm %r*% Wx) %r+% Mx
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




