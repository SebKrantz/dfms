
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


# Global variable in the package
Rcon <- matrix(c(2, -1, 0, 0, 0,
                 3, 0, -1, 0, 0,
                 2, 0, 0, -1, 0,
                 1, 0, 0, 0, -1), ncol = 5, byrow = TRUE)

init_cond_MQ <- function(X, X_imp, F_pc, v, n, r, p, TT, nq, rRi, rQi) {

  # .c(X, X_imp, F_pc, v, n, r, p, TT, nq, rRi, rQi) %=% list(X, X_imp, F_pc, v, n, r, p, TT, nq, rRi, rQi)
  rp <- r * p
  rC <- 5L # ncol(Rcon)
  pC <- max(p, rC)
  rpC <- r * pC
  nm <- n - nq

  # Observation equation -------------------------------
  # Static predictions (all.equal(unattrib(HDB(X_imp, F_pc)), unattrib(F_pc %*% t(v))))
  C <- cbind(v, matrix(0, n, rpC-r))
  rRcon <- kronecker(Rcon, diag(r))

  # Contemporaneous factors + lags
  FF <- do.call(cbind, lapply(0:(rC-1), function(i) F_pc[(rC - i):(TT - i), ]))

  # Now looping over the quarterly variables: at the end
  for (i in (nm+1):n) {
      x_i = X[rC:TT, i]
      nna = whichNA(x_i, invert = TRUE)
      if (length(nna) < ncol(FF) + 2L) x_i = X_imp[rC:TT, i]
      ff_i = FF[nna, ]
      x_i = x_i[nna] # Quarterly observations (no interpolation)
      Iff_i = ainv(crossprod(ff_i))
      Cc = Iff_i %*% crossprod(ff_i, x_i) # Coefficients from regressing quarterly observations on lagged factors
      # This is restricted least squares with restrictions: Rcon * C_0 = (q is 0)
      # The restrictions in Rcon (with -1 in the right places) make sense!
      tmp = tcrossprod(Iff_i, rRcon)
      Cc = Cc - tmp %*% ainv(rRcon %*% tmp) %*% rRcon %*% Cc
      C[i, 1:(rC*r)] = Cc # This replaces the corresponding row.
  }

  if(rRi) {
    # This computes residuals based on the new C matrix
    res <- X[rC:TT, ] - tcrossprod(FF, C[, 1:(rC*r)]) # residuals from static predictions
    R <- if(rRi == 2L) cov(res, use = "pairwise.complete.obs") else diag(fvar(res, na.rm = TRUE))
  } else R <- diag(n)
  R[(nm+1):n, (nm+1):n] <- 0 # Note: should have tol on diagonal?
  # diag(R)[diag(R) == 0] <- 1e-6

  # Note: should allow for additional zeros?
  C = cbind(C, rbind(matrix(0, nm, rC*nq), kronecker(t(c(1, 2, 3, 2, 1)), diag(nq))))

  # Transition equation -------------------------------
  rpC5nq <- rpC + 5 * nq
  var <- .VAR(F_pc, p)
  A <- Q <- matrix(0, rpC5nq, rpC5nq)
  A[1:r, 1:rp] <- t(var$A) # var$A is rp x r matrix
  A[(r+1L):rpC, 1:(rpC-r)] <- diag(rpC-r)
  A[(rpC+nq+1L):rpC5nq, (rpC+nq+1L):rpC5nq] <- diag(4*nq)

  Q[1:r, 1:r] <- switch(rQi + 1L, diag(r), diag(fvar(var$res)), cov(var$res))
  Q[(rpC+1):(rpC+nq), (rpC+1):(rpC+nq)] <- if(rRi == 2L)
      cov(res[, -seq_len(nm), drop = FALSE], use = "pairwise.complete.obs") else if(rRi == 1L)
      diag(fvar(res[, -seq_len(nm), drop = FALSE], na.rm = TRUE), nrow = nq) else diag(nq)
  diag(Q)[diag(Q) == 0] <- 1e-6 # Prevent singularity in Kalman Filter


  # Initial state and state covariance (P) ------------
  F_0 <- rep(0, rpC5nq) # BM14 uses zeros, DGR12 uses the first row of PC's. Both give more or less the same...
  # Kalman gain is normally A %*% t(A) + Q, but here A is somewhat tricky...
  # A_sparse <- as(A, "sparseMatrix")
  # M_sparse <- diag(rpC5nq^2) - kronecker(A_sparse, A_sparse) + 1e-6
  # P_0 <- solve(M_sparse) %*% unattrib(Q)
  M <- diag(rpC5nq^2) - kronecker(A, A)
  diag(M) <- diag(M) + 1e-4 # Ensure matrix is non-singular
  P_0 <- ainv(M) %*% unattrib(Q)
  dim(P_0) <- c(rpC5nq, rpC5nq)

  return(list(A = A, C = C, Q = Q, R = R, F_0 = F_0, P_0 = P_0))
}


# Mixed frequency with AR(1) idiosyncratic errors
# Based on MATLAB: EM_DFM_SS_idio_restrMQ.m InitCond()
# State vector: [factors(rpC), monthly_errors(nm), quarterly_error_lags(5*nq)]
init_cond_MQ_idio <- function(X, X_imp, F_pc, v, n, r, p, TT, nq, rRi, rQi, anymiss, tol) {

  # Dimensions
  rp <- r * p
  rC <- 5L # ncol(Rcon)
  pC <- max(p, rC)
  rpC <- r * pC
  nm <- n - nq
  state_dim <- rpC + nm + 5L * nq

  sr <- seq_len(r)
  srp <- seq_len(rp)
  snm <- seq_len(nm)

  # Observation equation -------------------------------
  # 1. Factor loadings with Rcon constraints for quarterly variables
  C <- cbind(v, matrix(0, n, rpC - r))
  rRcon <- kronecker(Rcon, diag(r))

  # Contemporaneous factors + lags for quarterly variable regression
  FF <- do.call(cbind, lapply(0:(rC-1L), function(i) F_pc[(rC - i):(TT - i), ]))

  # Restricted least squares for quarterly loadings
  for (i in (nm+1L):n) {
    x_i <- X[rC:TT, i]
    nna <- whichNA(x_i, invert = TRUE)
    if (length(nna) < ncol(FF) + 2L) x_i <- X_imp[rC:TT, i]
    ff_i <- FF[nna, , drop = FALSE]
    x_i <- x_i[nna]
    Iff_i <- ainv(crossprod(ff_i))
    Cc <- Iff_i %*% crossprod(ff_i, x_i)
    tmp <- tcrossprod(Iff_i, rRcon)
    Cc <- Cc - tmp %*% ainv(rRcon %*% tmp) %*% rRcon %*% Cc
    C[i, 1:(rC*r)] <- Cc
  }

  # Compute residuals for AR(1) estimation
  res <- X[rC:TT, ] - tcrossprod(FF, C[, 1:(rC*r)])
  resNaN <- res
  resNaN[is.na(X[rC:TT, ])] <- NA

  # Initial R for variance estimation
  R_init <- if(rRi) {
    if(rRi == 2L) cov(resNaN, use = "pairwise.complete.obs")
    else diag(fvar(resNaN, na.rm = TRUE))
  } else diag(n)

  # 2. Monthly error identity block: eyeN in MATLAB
  # C = [C, eyeN] where eyeN is N x NM (identity for monthly rows, zeros for quarterly)
  C <- cbind(C, rbind(diag(nm), matrix(0, nq, nm)))

  # 3. Quarterly temporal aggregation block
  # C = [C, [zeros(NM,5*NQ); kron(eye(NQ),[1 2 3 2 1])]]
  # Note: MATLAB kron(eye(NQ),[1 2 3 2 1]) treats [1 2 3 2 1] as row vector -> nq x 5*nq matrix
  C <- cbind(C, rbind(matrix(0, nm, 5L*nq), kronecker(diag(nq), matrix(c(1, 2, 3, 2, 1), nrow = 1))))

  # Transition equation -------------------------------
  # Monthly AR(1) estimation: BM and SM matrices
  BM <- SM <- numeric(nm)

  for (i in snm) {
    res_i <- resNaN[, i]
    # Find valid (non-NA) observations, excluding leading/trailing NAs
    valid <- !is.na(res_i)
    if(sum(valid) > 2L) {
      res_clean <- res[valid, i]
      n_clean <- length(res_clean)
      # AR(1) coefficient: rho = cov(e_t, e_{t-1}) / var(e_{t-1})
      BM[i] <- sum(res_clean[-n_clean] * res_clean[-1L]) /
               sum(res_clean[-n_clean]^2)
      # Innovation variance: var(e_t - rho * e_{t-1})
      SM[i] <- var(res_clean[-1L] - res_clean[-n_clean] * BM[i])
    } else {
      BM[i] <- 0.1 # Default initial value
      SM[i] <- R_init[i, i]
    }
  }
  # Stationary variance for monthly errors
  initViM <- SM / (1 - BM^2)

  # Quarterly initial values (from MATLAB lines 317-328)
  # sig_e = Rdiag(NM+1:N)/19 where 19 = 1^2 + 2^2 + 3^2 + 2^2 + 1^2
  sig_e <- diag(R_init)[(nm+1L):n] / 19
  rho0 <- 0.1 # Initial AR(1) coefficient for quarterly errors

  # Build A matrix: blkdiag(A_factors, BM, BQ)
  var <- .VAR(F_pc, p)
  A <- matrix(0, state_dim, state_dim)

  # Factor VAR block
  A[sr, srp] <- t(var$A)
  if(pC > 1L) A[(r+1L):rpC, 1:(rpC-r)] <- diag(rpC - r)

  # Monthly AR1 block (diagonal)
  monthly_idx <- (rpC+1L):(rpC+nm)
  A[monthly_idx, monthly_idx] <- diag(BM)

  # Quarterly AR1 lag structure: BQ = kron(eye(NQ), [[rho0, 0,0,0,0]; [eye(4), zeros(4,1)]])
  # Each 5x5 block shifts lags and updates current with AR1
  BQ_block <- rbind(c(rho0, rep(0, 4)), cbind(diag(4), rep(0, 4)))
  for (j in seq_len(nq)) {
    idx <- rpC + nm + (j-1L)*5L + 1:5
    A[idx, idx] <- BQ_block
  }

  # Build Q matrix: blkdiag(Q_factors, SM, SQ)
  Q <- matrix(0, state_dim, state_dim)
  Q[sr, sr] <- switch(rQi + 1L, diag(r), diag(fvar(var$res)), cov(var$res))

  # Monthly innovation variance (diagonal)
  Q[monthly_idx, monthly_idx] <- diag(SM)

  # Quarterly innovation variance: SQ = kron(diag((1-rho0^2)*sig_e), temp) where temp[1,1]=1
  # Only the (1,1) position of each 5x5 block has innovation variance
  for (j in seq_len(nq)) {
    idx <- rpC + nm + (j-1L)*5L + 1L
    Q[idx, idx] <- (1 - rho0^2) * sig_e[j]
  }

  # Prevent singularity
  diag(Q)[diag(Q) < 1e-6] <- 1e-6

  # Observation covariance R: fixed small value (kappa)
  # Actual observation error variance is modeled in state equation
  R <- tol * diag(n)

  # Initial state and state covariance (P) ------------
  F_0 <- rep(0, state_dim)

  # P_0 factor block: stationary VAR covariance
  srpC <- seq_len(rpC)
  tmp_A <- A[srpC, srpC, drop = FALSE]
  P_0_factors <- ainv(diag(rpC^2) - kronecker(tmp_A, tmp_A)) %*% unattrib(Q[srpC, srpC, drop = FALSE])
  dim(P_0_factors) <- c(rpC, rpC)

  # P_0 monthly errors: stationary AR1 variance
  P_0_monthly <- diag(initViM)

  # P_0 quarterly errors: initViQ = reshape(inv(eye((5*NQ)^2)-kron(BQ,BQ))*SQ(:),5*NQ,5*NQ)
  # For each quarterly variable, compute stationary covariance
  BQ_kron <- kronecker(BQ_block, BQ_block)
  P_0_block_inv <- ainv(diag(25) - BQ_kron)

  P_0 <- matrix(0, state_dim, state_dim)
  P_0[srpC, srpC] <- P_0_factors
  P_0[monthly_idx, monthly_idx] <- P_0_monthly

  for (j in seq_len(nq)) {
    idx <- rpC + nm + (j-1L)*5L + 1:5
    SQ_j <- matrix(0, 5, 5)
    SQ_j[1, 1] <- (1 - rho0^2) * sig_e[j]
    P_0[idx, idx] <- matrix(P_0_block_inv %*% c(SQ_j), 5, 5)
  }

  return(list(A = A, C = C, Q = Q, R = R, F_0 = F_0, P_0 = P_0))
}