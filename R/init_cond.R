
init_cond <- function(X, F_pc, v, n, r, p, BMl, rRi, rQi) {
  rp <- r * p
  sr <- seq_len(r)
  srp <- seq_len(rp)

  # Observation equation -------------------------------
  # Static predictions (all.equal(unattrib(HDB(X_imp, F_pc)), unattrib(F_pc %*% t(v))))
  C <- cbind(v, matrix(0, n, rp-r))

  if(rRi) {
    res <- X - tcrossprod(F_pc, v) # residuals from static predictions
    R <- if(rRi == 2L) cov(res, use = "pairwise.complete.obs") else diag(fvar(res, na.rm = TRUE))
  } else R <- diag(n)

  # Transition equation -------------------------------
  var <- .VAR(F_pc, p)
  A <- rbind(t(var$A), diag(1, rp-r, rp)) # var$A is rp x r matrix
  Q <- matrix(0, rp, rp)
  Q[sr, sr] <- switch(rQi + 1L, diag(r),  diag(fvar(var$res)), cov(var$res))

  # Initial state and state covariance (P) ------------
  F_0 <- if(isTRUE(BMl)) rep(0, rp) else var$X[1L, ] # BM14 uses zeros, DGR12 uses the first row of PC's. Both give more or less the same...
  # Kalman gain is normally A %*% t(A) + Q, but here A is somewhat tricky...
  P_0 <- ainv(diag(rp^2) - kronecker(A,A)) %*% unattrib(Q)
  dim(P_0) <- c(rp, rp)

  return(list(A = A, C = C, Q = Q, R = R, F_0 = F_0, P_0 = P_0))
}


init_cond_idio_ar1 <- function(X, F_pc, v, n, r, p, BMl, rRi, rQi, anymiss, tol) {
  rp <- r * p
  sr <- seq_len(r)
  srp <- seq_len(rp)
  end <- (rp+1L):(rp+n)

  # Observation equation -------------------------------
  C <- cbind(v, matrix(0, n, rp-r), diag(n))

  if(rRi) {
    res <- X - tcrossprod(F_pc, v) # residuals from static predictions
    # If AR(1) residuals, need to estimate coefficient and clean residuals from autocorrelation
    res_AC1 <- AC1(res, anymiss)
    res <- res[-1L, ] %-=% setop(res[-nrow(res), ], "*", res_AC1, rowwise = TRUE)
    R <- if(rRi == 2L) cov(res, use = "pairwise.complete.obs") else diag(fvar(res, na.rm = TRUE))
  } else R <- diag(n)

  # Transition equation -------------------------------
  var <- .VAR(F_pc, p)
  P_0 <- A <- Q <- matrix(0, rp+n, rp+n)
  A[sr, srp] <- t(var$A)
  if(p > 1L) A[(r+1L):rp, srp] <- diag(1, rp-r, rp)
  A[end, end] <- diag(res_AC1) # Estimates of residual autocorrelation
  Q[sr, sr] <- switch(rQi + 1L, diag(r),  diag(fvar(var$res)), cov(var$res))
  Q[end, end] <- R # Observation covariance is estimated in State Equation

  # Initial state and state covariance (P) ------------
  # BM14 uses zeros, DGR12 uses the first row of PC's. Both give more or less the same...
  F_0 <- c(if(isTRUE(BMl)) rep(0, rp) else var$X[1L, ], rep(0, n))
  # Kalman gain is normally A %*% t(A) + Q, but here A is somewhat tricky...
  tmp = A[srp, srp, drop = FALSE]
  P_0[srp, srp] <- ainv(diag(rp^2) - kronecker(tmp, tmp)) %*% unattrib(Q[srp, srp, drop = FALSE])
  P_0[end, end] <- diag(1/(1-res_AC1^2) * diag(R))
  if(rRi == 2L) R <- diag(n)
  diag(R) <- tol # The actual observation covariance is a very small fixed number (kappa)

  return(list(A = A, C = C, Q = Q, R = R, F_0 = F_0, P_0 = P_0))
}

