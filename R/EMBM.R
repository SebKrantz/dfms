#' Expectation Maximization Algorithm following Banbura and Modungo (2014)
#' @inheritParams SKF
#' @inheritParams impNA_MA
#' @param XWO \code{X} with missing values set to 0.
#' @param rQi,rRi restrictions on Q and R passed to DFM(), and turned into integers such that identity = 0L, diagonal = 1L, none = 2L.
#' The other parameters are matrix dimensions, which do not need to be recalculated every iteration, and inputs to avoid repeated computations of
#' kronecker products in the EM loop. see DFM() code where these objects are generated.
#'
#' where
#' @noRd
EMstepBMOPT <- function(X, A, C, Q, R, F_0, P_0, XW0, W, n, r, sr, TT, dgind, dnkron, dnkron_ind, rQi, rRi) {

  kfs_res = SKFS(X, A, C, Q, R, F_0, P_0, TRUE)
  # kfs_res = tryCatch(SKFS(X, A, C, Q, R, F_0, P_0), error = function(e) return(list(NULL, X, A, C, Q, R, F_0, P_0)))
  # if(is.null(kfs_res[[1]])) return(kfs_res)
  Zsmooth = kfs_res$F_smooth
  Vsmooth = kfs_res$P_smooth
  VVsmooth = kfs_res$PPm_smooth
  Zsmooth0 = kfs_res$F_smooth_0
  Vsmooth0 = kfs_res$P_smooth_0
  # return(kfs_res[.c(F_smooth, P_smooth, F_smooth_0, P_smooth_0, PPm_smooth)])

  tmp = rbind(Zsmooth0, Zsmooth[-TT,, drop = FALSE])
  tmp2 = sum3(Vsmooth[,, -TT, drop = FALSE])
  EZZ = crossprod(Zsmooth) %+=% (tmp2 + Vsmooth[,, TT]) # E(Z'Z)
  EZZ_BB = crossprod(tmp) %+=% (tmp2 + Vsmooth0)        # E(Z(-1)'Z_(-1))
  EZZ_FB = crossprod(Zsmooth, tmp) %+=% sum3(VVsmooth)  # E(Z'Z_(-1))

  A_new = A
  A_new[sr, ] = EZZ_FB[sr, , drop = FALSE] %*% ainv(EZZ_BB)
  Q_new = Q

  if(rQi) {
    Qsr = (EZZ[sr, sr] - tcrossprod(A_new[sr, , drop = FALSE], EZZ_FB[sr, , drop = FALSE])) / TT
    Q_new[sr, sr] = if(rQi == 2L) Qsr else diag(diag(Qsr))
  } else Q_new[sr, sr] = diag(r)

  Zsmooth = Zsmooth[, sr, drop = FALSE]
  Vsmooth = Vsmooth[sr, sr,, drop = FALSE]

  # E(X'X) & E(X'Z)
  denom = numeric(n*r^2)
  nom = matrix(0, n, r)
  for (t in 1:TT) {
    tmp = Zsmooth[t, ]
    nom %+=% tcrossprod(XW0[t, ], tmp)
    tmp2 = tcrossprod(tmp) + Vsmooth[,, t]
    dim(tmp2) = NULL
    denom %+=% tcrossprod(tmp2, !W[t, ])
  }

  dim(denom) = c(r, r, n)
  dnkron[dnkron_ind] = aperm.default(denom, c(1L, 3L, 2L))
  C_new = solve.default(dnkron, unattrib(nom)) # ainv -> slower...
  dim(C_new) = c(n, r)

  if(rRi) {
    R_new = matrix(0, n, n)
    Rdg = R[dgind]
    for (t in 1:TT) {
      nanYt = W[t, ]
      tmp = C_new * !nanYt
      R[dgind] = Rdg * nanYt # If R is not diagonal
      tmp2 = tmp %*% tcrossprod(Vsmooth[,, t], tmp)
      tmp2 %+=% tcrossprod(XW0[t, ] - tmp %*% Zsmooth[t, ])
      tmp2 %+=% R # If R is not diagonal
      # tmp2[dgind] = tmp2[dgind] + (Rdg * nanYt) # If R is diagonal...
      R_new %+=% tmp2
    }
    if(rRi == 2L) { # Unrestricted
      R_new %/=% TT
      R_new[R_new < 1e-7] = 1e-7
    } else { # Diagonal
      # R_new = diag(R_new[dgind] / TT)
      RR = R_new[dgind] / TT
      RR[RR < 1e-7] = 1e-7 # RR(RR<1e-2) = 1e-2;
      R_new = diag(RR)
    }
  } else R_new = diag(n)

  C[, sr] = C_new

  return(list(A = A_new,
              C = C,
              Q = Q_new,
              R = R_new,
              F_0 = drop(Zsmooth0),
              P_0 = Vsmooth0,
              loglik = kfs_res$loglik))
}



# Estep(X, C, Q, R, A, F_0, P_0)
# F_0 = F_0
# P_0 = P_0
# EMstepBM <- function(X, A, C, Q, R, F_0, P_0) { # [C_new, R_new, A_new, Q_new, F_0, P_0, loglik]
#
#   # Kalman filer BM:
#   ###  Inputs:
#   ###    X: k-by-nobs matrix of input data
#   ###    A: m-by-m transition matrix
#   ###    C: k-by-m measurement matrix
#   ###    Q: m-by-m covariance matrix for transition equation residuals (mu_t)
#   ###    R: k-by-k covariance for measurement matrix residuals (e_t)
#   ###    F_0: 1-by-m vector, initial value of state
#   ###    P_0: m-by-m matrix, initial value of factor covariance matrix
#   ###
#   ###  Outputs:
#   ###    Fsmooth: k-by-(nobs+1) matrix, smoothed factor estimates (i.e. Fsmooth[, t + 1] = F_t|T)
#   ###    Psmooth: k-by-k-by-(nobs+1) array, smoothed factor covariance matrices (i.e. Psmooth[, , t + 1) = Cov(F_t|T))
#   ###    PPsmooth: k-by-k-by-nobs array, lag 1 factor covariance matrices (i.e. Cov(F_t, F_t-1|T))
#   ###    loglik: scalar, log-likelihood
#
#   c("T", "n") %=% dim(X)
#   c("Fsmooth", "Psmooth", "PPsmooth", "loglik") %=% SKFS(X, C, Q, R, A, F_0, P_0)
#   r <- dim(Fsmooth)[2L]
#   # TODO: Psmooth should have 1 more obs..
#   ncv <- dim(Psmooth)[3L]
#
#   # Conditional moments of the factors (F) in eq. 6 and 8
#   # TODO:: Crossprod right?
#   EFF = crossprod(Fsmooth[-1L, ]) + rowSums(Psmooth[,,-1L], dims = 2L)   # E(F'F)
#   EFF_BB = crossprod(Fsmooth[-r, ]) + rowSums(Psmooth[,,-ncv], dims = 2L)         # E(F(-1)'F_(-1))
#   EFF_FB = crossprod(Fsmooth[-1L, ], Fsmooth[-r, ]) + rowSums(PPsmooth, dims = 2L) # E(F'F_(-1))
#
#   A_new = A
#   # Equation 6: Estimate VAR(p) for factor
#   A_new[1:r, ] = EFF_FB[1:r, ] %*% ainv(EFF_BB)
#   Q_new = Q;
#   # Equation 8: Covariance matrix of residuals of VAR
#   Q_new[1:r, 1:r] = (EFF[1:r, 1:r] - tcrossprod(A_new[1:r, ], EFF_FB[1:r, ])) / T # matrix division??
#
#
#   # E(Y'Y) & E(Y'F)
#   nanY = is.na(X)
#   X[nanY] = 0
#
#   # TODO: Use fcumsum here:
#   # Eq. 11:
#   denom = matrix(0, n*r, n*r)
#   nom = matrix(0, n, r)
#   # TODO: Also add lags here where necessary...
#   for (t in 1:T) {
#     nanYt = diag(!nanY[t, ]) # What does this give exactly ??
#     denom = denom + kronecker(tcrossprod(Fsmooth[t, 1:r]) + Psmooth[1:r, 1:r, t], nanYt)
#     nom = nom + tcrossprod(X[t, ], Fsmooth[t, 1:r])
#   }
#
#   C_new = ainv(denom) %*% unattrib(nom)
#   dim(C_new) = c(n, r)
#
#   # Eq 12:
#   R_new = matrix(0, n, n)
#   for (t in 1:T) {
#     nanYt = diag(!nanY[t, ])
#     R_new = R_new + tcrossprod(X[t, ] - nanYt %*% C_new %*% Fsmooth[t, 1:r]) +
#       nanYt %*% C_new %*% Psmooth[1:r, 1:r, t] %*% crossprod(C_new, nanYt) +
#       (diag(n) - nanYt) %*% R %*% (diag(n) - nanYt)
#   }
#
#   R_new = R_new / T
#   RR = diag(R_new) # RR(RR<1e-2) = 1e-2;
#   R_new = diag(RR)
#
#   # Initial conditions
#   F_0 = Fsmooth[, ] # zeros(size(Fsmooth,1),1); %
#   P_0 = Psmooth[,, 1]
#
#   list(A = A_new,
#        C= C_new,
#        Q = Q_new,
#        R = R_new,
#        F_0 = F_0,
#        P_0 = P_0,
#        loglik = loglik)
# }
