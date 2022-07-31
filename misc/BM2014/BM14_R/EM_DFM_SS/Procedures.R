# r = 3; p = 3; optNaN = options; x = EuStockMarkets
InitCond <- function(x, r, p, optNaN) { # [A, C, Q, R, initZ, initV]

  rp = r * p
  sr = 1:r

  c('xBal', 'indNaN') %=% remNaNs(x, optNaN)
  c('T', 'N') %=% dim(xBal)

  #  Eigenval decomp of cov(x) = VDV', only r largest evals
  eigen_decomp = eigen(cov(xBal, use = "complete.obs"), symmetric = TRUE)
  v = eigen_decomp$vectors[, sr]
  d = eigen_decomp$values[sr]

  # Static predictions
  chi = xBal %*% tcrossprod(v)

  res = xBal - chi
  res[indNaN] = NA_real_

  # Observation equation
  R = diag(fvar(res))

  # Observation equation
  C = cbind(v, matrix(0, N, rp - r))

  # Transition equation

  # Estimate A & Q from stacked F(t) = A*F(t-1) + e(t);
  F = xBal %*% v
  Z = NULL # This is capital Z
  for (i in 1:p) Z = cbind(Z, F[(p-i+1L):(T-i), ]) # stacked regressors (lagged SPC)
  z = F[(p+1L):T, ]

  # run the var chi(t) = A*chi(t-1) + e(t);
  A = matrix(0, rp, rp)
  A_temp = ainv(crossprod(Z)) %*% crossprod(Z, z) # cinv
  A[sr, 1:rp] = t(A_temp)

  A[(r+1L):rp, 1:(rp-r)] = diag(rp-r)

  Q = matrix(0, rp, rp)
  e = z - Z %*% A_temp # VAR residuals
  Q[sr, sr] = cov(e) # VAR covariance matrix

  # Initial conditions
  initZ = matrix(0, ncol(Z)) # [randn(1,r*(nlag+1))]';
  initV = matrix(ainv(diag(rp^2) - kronecker(A,A)) %*% unattrib(Q), rp, rp) # solve(, tol = .Machine$double.eps/1e10)

  return(list(A = A, C = C, Q = Q, R = R, initZ = initZ, initV = initV))
}

# sum3 <- function(a) {
#   d <- dim(a)
#   if(length(d) != 3L) stop("a needs to be a 3D array")
#   dim(a) <- c(d[1L] * d[2L], d[3L])
#   sa <- rowSums(a, dims = 2)
# }

EMstep <- function(y, A, C, Q, R, Z_0, V_0, r) { # [C_new, R_new, A_new, Q_new, Z_0, V_0, loglik]

  c("n", "T") %=% dim(y)
  c("Zsmooth", "Vsmooth", "VVsmooth", "loglik") %=% runKF(y, A, C, Q, R, Z_0, V_0)
  nc = dim(Zsmooth)[2L]
  ncv = dim(Vsmooth)[3L]
  sr = 1:r


  EZZ = tcrossprod(Zsmooth[, -1L]) + rowSums(Vsmooth[,, -1L], dims = 2L)             # E(Z'Z)
  EZZ_BB = tcrossprod(Zsmooth[, -nc]) + rowSums(Vsmooth[,, -ncv], dims = 2L)         # E(Z(-1)'Z_(-1))
  EZZ_FB = tcrossprod(Zsmooth[, -1L], Zsmooth[, -nc]) + rowSums(VVsmooth, dims = 2L) # E(Z'Z_(-1))

  A_new = A
  A_new[sr, ] = EZZ_FB[sr, ] %*% ainv(EZZ_BB) # solve(, tol = .Machine$double.eps/1e10)
  Q_new = Q;
  Q_new[sr, sr] = (EZZ[sr, sr] - tcrossprod(A_new[sr, ], EZZ_FB[sr, ])) / T


  # E(Y'Y) & E(Y'Z)
  nanY = is.na(y)
  y[nanY] = 0

  denom = matrix(0, n*r, n*r)
  nom = matrix(0, n, r)
  for (t in 1:T) {
      nanYt = diag(as.double(!nanY[, t]))
      tmp = Zsmooth[sr, t+1L]
      denom = denom + kronecker(tcrossprod(tmp) + Vsmooth[sr, sr, t+1L], nanYt)
      nom = nom + tcrossprod(y[, t], tmp)
  }

  C_new = ainv(denom) %*% unattrib(nom) # solve(, tol = .Machine$double.eps/1e10)
  dim(C_new) = c(n, r)

  R_new = matrix(0, n, n)
  In = diag(n)
  for (t in 1:T) {
      nanYt = diag(as.double(!nanY[, t]))
      tmp = nanYt %*% C_new
      tmp2 = In - nanYt
      R_new = R_new + tcrossprod(y[, t] - tmp %*% Zsmooth[sr, t+1L]) +
              tmp %*% Vsmooth[sr, sr, t+1L] %*% crossprod(C_new, nanYt) +
              tmp2 %*% R %*% tmp2
  }

  R_new = R_new / T
  RR = diag(R_new) # RR(RR<1e-2) = 1e-2;
  R_new = diag(RR)

  # Initial conditions
  Z_0 = Zsmooth[, 1L] # zeros(size(Zsmooth,1),1); %
  V_0 = Vsmooth[,, 1L]

  return(list(C_new = C_new,
              R_new = R_new,
              A_new = A_new,
              Q_new = Q_new,
              Z_0 = Z_0,
              V_0 = V_0,
              loglik = loglik))
}


# threshold = thresh
em_converged <- function(loglik, previous_loglik, threshold = 1e-4, check_increased = TRUE) { # [converged, decrease]
# EM_CONVERGED Has EM converged?
# [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
#
# We have converged if the slope of the log-likelihood function falls below 'threshold',
# i.e., |f(t) - f(t-1)| / avg < threshold,
# where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
# 'threshold' defaults to 1e-4.
#
# This stopping criterion is from Numerical Recipes in C p423
#
# If we are doing MAP estimation (using priors), the likelihood can decrase,
# even though the mode of the posterior is increasing.

  converged = FALSE
  decrease = FALSE
  eps = .Machine$double.eps # Matlab floating-point relative accuracy

  if(check_increased) {
    test <- loglik - previous_loglik < -1e-3
    if(is.finite(test) && test) { # allow for a little imprecision
       sprintf('******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik)
       decrease = TRUE
    }
  }

  delta_loglik = abs(loglik - previous_loglik)
  avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2
  test = delta_loglik / avg_loglik < threshold
  if(is.finite(test) && test) converged = TRUE

  return(list(converged = converged, decrease = decrease))
}
