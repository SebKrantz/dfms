library(collapse) # 1.8+
library(roll)
library(DFM)

################################################################################
# This is an R-optimized version of the original BM2014 code,
# that also introduces clear and consistent notation.
#
# The model is        Y_t   = C * Z_t + r_t ~ N(0, R)
#                     Z_t+1 = A * Z_t + q_t ~ N(0, Q)
#                     Cov(Z_t) = V; Cov(Z_t+1,Z_t) = VV
#
# It is noteworthy that this code is significantly faster than BM2014's
# Matlab code, particularly because of the reduced use of kronecker products
################################################################################

setwd("~/Documents/R/DFM")
source("misc/BM2014/BM14_R/EM_DFM_SS_OPT/remNaNs_spline.R")
source("misc/BM2014/BM14_R/EM_DFM_SS_OPT/runKF.R")
source("misc/BM2014/BM14_R/EM_DFM_SS_OPT/Procedures.R")

# TODO: make it work with just 1 series?
# Also make sure cinv solution is robust!!
# Also try to use FKF package...
# FKF::fkf

EM_DFM_SS_OPT <- function(X, r, p = 1L, max_iter = 100L, thresh = 1e-4, na.method = 2L, ma.k = 3L) {

  # --------------------------------------------------------------------------
  # Preparation of the data
  #--------------------------------------------------------------------------

  X = qM(X)
  c("T", "N") %=% dim(X)

  # Standardise data
  Mx = fmean(X)
  Wx = fsd(X)
  X_STD = fscale(X)

  # --------------------------------------------------------------------------
  # Initial Conditions
  #--------------------------------------------------------------------------

  # Removing missing values (for initial estimators)
  optNaN = list()
  optNaN$method = na.method # See remNaNs_spline.R
  optNaN$k = ma.k           # order of the moving average for replacing the missing observations

  # return(InitCond(X_STD, r, p, optNaN))
  c("A", "C", "Q", "R", "Z_0", "V_0") %=% InitCond(X_STD, r, p, optNaN)

  # some auxiliary variables for the iterations
  previous_loglik = -Inf
  num_iter = 0L
  LL = -Inf
  converged = FALSE

  # Y for the final estimation WITH missing data
  Y = t(X_STD)

  #--------------------------------------------------------------------------
  # THE EM LOOP
  #--------------------------------------------------------------------------

  # Remove the leading and ending nans for the estimation
  optNaN$method = 3
  Y_narm = t(remNaNs(X_STD, optNaN)$X)
  dimnames(Y_narm) = NULL
  dnkron = matrix(1, r, r) %x% diag(nrow(Y_narm)) # Used to be inside EMstep
  dnkron_ind = whichv(dnkron, 1)
  T = ncol(Y_narm)
  rp = r*p

  # Matrices for Kalman Filter and Smoother
  S = list(ZT  = matrix(NA_real_, rp, T),
           VT  = array(NA_real_, c(rp, rp, T)),
           ZT_0 = matrix(NA_real_, rp, T + 1L),
           VT_0 = array(NA_real_, c(rp, rp, T + 1L)),
           loglik = 0,
           ZsmoothT  = matrix(0, rp, T + 1L),
           VsmoothT  = array(0, c(rp, rp, T + 1L)),
           VVsmoothT = array(0, c(rp, rp, T)))

  while(num_iter < max_iter && !converged) {

    c("A", "C", "Q", "R", "Z_0", "V_0", "loglik") %=% EMstep_OPT(Y_narm, A, C, Q, R, Z_0, V_0, r, dnkron, dnkron_ind, S)

    # Checking convergence
    c("converged", "decrease") %=% em_converged(loglik, previous_loglik, thresh, TRUE)

    LL = c(LL, loglik)
    previous_loglik = loglik
    num_iter =  num_iter + 1L
  }

  # Final run of the Kalman filter, on data with leaning and ending nans
  Zsmooth = t(runKF(Y, A, C, Q, R, Z_0, V_0, S)$Zsmooth)

  #--------------------------------------------------------------------------
  #   Loading the structure with the results
  #--------------------------------------------------------------------------
  Res = list()
  Res$x_sm = tcrossprod(Zsmooth[-1L, , drop = FALSE], C)      # Data Prediction (standardized)
  Res$X_sm = setop(Res$x_sm %r*% Wx, "+", Mx, rowwise = TRUE) # Data Prediction at original scale
  Res$F = Zsmooth[-1L, 1:r, drop = FALSE]                     # Smoothed Factor estimates
  Res$C = C
  Res$R = R
  Res$A = A
  Res$Q = Q
  Res$Z_0 = Z_0 # Smoothed factor estimate before final KF run
  Res$V_0 = V_0 # Smoothed factor covariance estimate before final KF run
  Res$r = r
  Res$p = p
  Res$Mx = Mx
  Res$Wx = Wx
  Res$LL = LL  # Sequence of log-likelihoods

  return(Res)
}




