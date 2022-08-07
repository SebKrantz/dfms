
EMstepBMOPT <- function(X, A, C, Q, R, F0, P0, XW0, W, n, r, sr, T, dnkron, dnkron_ind) {

  kfs_res = fKFS(X, A, C, Q, R, F0, P0, 2L)
  # kfs_res = tryCatch(fKFS(X, A, C, Q, R, F0, P0), error = function(e) return(list(NULL, X, A, C, Q, R, F0, P0)))
  # if(is.null(kfs_res[[1]])) return(kfs_res)
  Zsmooth = kfs_res$F_smooth
  Vsmooth = kfs_res$P_smooth
  VVsmooth = kfs_res$PPm_smooth
  Zsmooth0 = kfs_res$F_smooth_0
  Vsmooth0 = kfs_res$P_smooth_0

  tmp = rbind(Zsmooth0, Zsmooth[-T,, drop = FALSE])
  EZZ = crossprod(Zsmooth) %+=% rowSums(Vsmooth, dims = 2L)                                     # E(Z'Z)
  EZZ_BB = crossprod(tmp) %+=% (rowSums(Vsmooth[,, -T, drop = FALSE], dims = 2L) + Vsmooth0)    # E(Z(-1)'Z_(-1))
  EZZ_FB = crossprod(Zsmooth, tmp) %+=% rowSums(VVsmooth, dims = 2L)                            # E(Z'Z_(-1))

  A_new = A
  A_new[sr, ] = EZZ_FB[sr, , drop = FALSE] %*% ainv(EZZ_BB)
  Q_new = Q
  Q_new[sr, sr] = (EZZ[sr, sr] - tcrossprod(A_new[sr, , drop = FALSE], EZZ_FB[sr, , drop = FALSE])) / T

  # E(X'X) & E(X'Z)
  denom = numeric(n*r^2)
  nom = matrix(0, n, r)
  for (t in 1:T) {
    tmp = Zsmooth[t, sr]
    nom %+=% tcrossprod(XW0[t, ], tmp)
    tmp2 = tcrossprod(tmp) + Vsmooth[sr, sr, t]
    dim(tmp2) = NULL
    denom %+=% tcrossprod(tmp2, !W[t, ])
  }

  dim(denom) = c(r, r, n)
  dnkron[dnkron_ind] = aperm(denom, c(1L, 3L, 2L))
  C_new = ainv(dnkron) %*% unattrib(nom) # cinv
  dim(C_new) = c(n, r)

  R_new = matrix(0, n, n)
  for (t in 1:T) {
    nanYt = W[t, ]
    tmp = C_new * !nanYt
    R2 = R
    diag(R2) = diag(R) * nanYt
    tmp2 = tmp %*% tcrossprod(Vsmooth[sr, sr, t], tmp)
    tmp2 %+=% tcrossprod(XW0[t, ] - tmp %*% Zsmooth[t, sr])
    tmp2 %+=% R2
    R_new %+=% tmp2
  }

  R_new %/=% T
  RR = diag(R_new) # RR(RR<1e-2) = 1e-2;
  R_new = diag(RR)

  C[, sr] = C_new

  # A = A_new
  # C = C
  # Q = Q_new
  # R = R_new
  # F0 = drop(Zsmooth0)
  # P0 = Vsmooth0
  # A; C; Q; diag(R); F0; P0

  return(list(A = A_new,
              C = C,
              Q = Q_new,
              R = R_new,
              F0 = drop(Zsmooth0),
              P0 = Vsmooth0,
              loglik = kfs_res$loglik))
}



# Estep(X, C, Q, R, A, F0, P0)
# F0 = F0
# P0 = P0
# EMstepBM <- function(X, A, C, Q, R, F0, P0) { # [C_new, R_new, A_new, Q_new, F0, P0, loglik]
#
#   # Kalman filer BM:
#   ###  Inputs:
#   ###    X: k-by-nobs matrix of input data
#   ###    A: m-by-m transition matrix
#   ###    C: k-by-m measurement matrix
#   ###    Q: m-by-m covariance matrix for transition equation residuals (mu_t)
#   ###    R: k-by-k covariance for measurement matrix residuals (e_t)
#   ###    F0: 1-by-m vector, initial value of state
#   ###    P0: m-by-m matrix, initial value of factor covariance matrix
#   ###
#   ###  Outputs:
#   ###    Fsmooth: k-by-(nobs+1) matrix, smoothed factor estimates (i.e. Fsmooth[, t + 1] = F_t|T)
#   ###    Psmooth: k-by-k-by-(nobs+1) array, smoothed factor covariance matrices (i.e. Psmooth[, , t + 1) = Cov(F_t|T))
#   ###    PPsmooth: k-by-k-by-nobs array, lag 1 factor covariance matrices (i.e. Cov(F_t, F_t-1|T))
#   ###    loglik: scalar, log-likelihood
#
#   c("T", "n") %=% dim(X)
#   c("Fsmooth", "Psmooth", "PPsmooth", "loglik") %=% fKFS(X, C, Q, R, A, F0, P0)
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
#   F0 = Fsmooth[, ] # zeros(size(Fsmooth,1),1); %
#   P0 = Psmooth[,, 1]
#
#   list(A = A_new,
#        C= C_new,
#        Q = Q_new,
#        R = R_new,
#        F0 = F0,
#        P0 = P0,
#        loglik = loglik)
# }
