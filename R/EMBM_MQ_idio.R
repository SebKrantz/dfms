# EM step for Mixed-Frequency with Idiosyncratic AR(1) Errors
# Based on MATLAB: EM_DFM_SS_idio_restrMQ.m EMstep()
#
# State vector: [factors(rpC), monthly_errors(nm), quarterly_error_lags(5*nq)]
#
# Parameters:
# X: T x n matrix of observed data (standardized, with NAs)
# A: state_dim x state_dim transition matrix
# C: n x state_dim observation matrix
# Q: state_dim x state_dim state covariance matrix
# R: n x n observation covariance matrix (fixed small value kappa)
# Z_0 (F_0): state_dim x 1 initial state vector
# V_0 (P_0): state_dim x state_dim initial state covariance
# XW0: T x n data with NAs replaced by zeros
# NW: T x n logical matrix (TRUE = observed, FALSE = missing)
# dgind: diagonal indices for n x n matrix
# dnkron: pre-allocated Kronecker product matrix for C update
# dnkron_ind: indices for updating dnkron
# r: number of factors
# p: number of lags in factor VAR
# R_mat: constraint matrix for quarterly loadings (kronecker(Rcon, diag(r)))
# n: total number of series
# nq: number of quarterly series
# sr: 1:r sequence
# TT: number of time periods
# rQi: restriction on Q (0=identity, 1=diagonal, 2=unrestricted)
# rRi: restriction on R (0=identity, 1=diagonal, 2=unrestricted)

EMstepBMMQidio <- function(X, A, C, Q, R, Z_0, V_0, XW0, NW, dgind, dnkron, dnkron_ind, r, p, R_mat, n, nq, sr, TT, rQi, rRi) {

  # Dimensions
  rC <- 5L
  pC <- max(p, rC)
  rp <- r * p
  rpC <- r * pC
  nm <- n - nq
  rCr <- rC * r  # = 5 * r (columns for quarterly factor loadings)
  state_dim <- rpC + nm + 5L * nq

  # Index sequences
  snm <- seq_len(nm)
  srp <- seq_len(rp)
  srpC <- seq_len(rpC)

  # State block indices
  monthly_idx <- (rpC + 1L):(rpC + nm)  # Monthly AR1 errors
  quarterly_all_idx <- (rpC + nm + 1L):state_dim  # All quarterly error lags (5*nq)
  # Indices for first lag of each quarterly error (for AR1 coefficient updates)
  quarterly_first_idx <- rpC + nm + seq(1L, 5L * nq, by = 5L)

  # Run Kalman filter with current parameter estimates
  kfs_res <- SKFS(X, A, C, Q, R, Z_0, V_0, TRUE)

  # Extract Kalman smoother output
  Zsmooth <- kfs_res$F_smooth
  Vsmooth <- kfs_res$P_smooth
  VVsmooth <- kfs_res$PPm_smooth
  Zsmooth0 <- kfs_res$F_smooth_0
  Vsmooth0 <- kfs_res$P_smooth_0
  loglik <- kfs_res$loglik

  # ==========================================
  # Factor expectations (indices 1:rpC)
  # ==========================================
  tmp <- rbind(Zsmooth0[srpC], Zsmooth[-TT, srpC, drop = FALSE])
  tmp2 <- sum3(Vsmooth[srpC, srpC, -TT, drop = FALSE])
  EZZ <- crossprod(Zsmooth[, srpC, drop = FALSE]) %+=% (tmp2 + Vsmooth[srpC, srpC, TT])  # E(Z'Z)
  EZZ_BB <- crossprod(tmp) %+=% (tmp2 + Vsmooth0[srpC, srpC])  # E(Z(-1)'Z(-1))
  EZZ_FB <- crossprod(Zsmooth[, srpC, drop = FALSE], tmp) %+=% sum3(VVsmooth[srpC, srpC, , drop = FALSE])  # E(Z'Z(-1))

  # ==========================================
  # Error expectations (for AR1 coefficient updates)
  # MATLAB: EZZ2, EZZ_BB2, EZZ_FB2 over rpC+1:end
  # These are diagonal matrices - we compute only diagonal elements
  # ==========================================
  # Indices for all errors that have AR1 dynamics:
  # - Monthly errors: positions rpC+1 to rpC+nm
  # - Quarterly errors first lag: positions rpC+nm+1, rpC+nm+6, rpC+nm+11, ... (every 5th)
  all_error_idx <- c(monthly_idx, quarterly_first_idx)
  n_errors <- length(all_error_idx)

  tmp_e <- rbind(Zsmooth0[all_error_idx], Zsmooth[-TT, all_error_idx, drop = FALSE])
  tmp2_e <- sum3(Vsmooth[all_error_idx, all_error_idx, -TT, drop = FALSE])

  EZZ_u <- diag(crossprod(Zsmooth[, all_error_idx, drop = FALSE])) +
           diag(tmp2_e + Vsmooth[all_error_idx, all_error_idx, TT])  # E(Z'Z)
  EZZ_BB_u <- diag(crossprod(tmp_e)) +
              diag(tmp2_e + Vsmooth0[all_error_idx, all_error_idx])  # E(Z(-1)'Z(-1))
  EZZ_FB_u <- diag(diag(crossprod(Zsmooth[, all_error_idx, drop = FALSE], tmp_e))) +
              diag(diag(rowSums(VVsmooth[all_error_idx, all_error_idx, , drop = FALSE], dims = 2L)))  # E(Z'Z(-1))

  # ==========================================
  # Update A and Q matrices
  # ==========================================
  A_new <- A
  Q_new <- Q

  # Factor transition matrix update: A(1:r, 1:rp)
  A_new[sr, srp] <- EZZ_FB[sr, srp, drop = FALSE] %*% ainv(EZZ_BB[srp, srp, drop = FALSE])

  # Factor covariance update: Q(1:r, 1:r)
  if(rQi) {
    Qsr <- (EZZ[sr, sr] - tcrossprod(A_new[sr, srp, drop = FALSE], EZZ_FB[sr, srp, drop = FALSE])) / TT
    Q_new[sr, sr] <- if(rQi == 2L) Qsr else diag(diag(Qsr))
  } else {
    Q_new[sr, sr] <- diag(r)
  }

  # Error AR1 coefficient and variance updates
  # MATLAB: A_new2 = EZZ_FB2 * diag(1./diag((EZZ_BB2)))
  #         Q_new2 = (EZZ2 - A_new2*EZZ_FB2') / T
  A_new2 <- EZZ_FB_u / EZZ_BB_u
  Q_new2 <- (EZZ_u - A_new2 * EZZ_FB_u) / TT

  # Monthly error updates: A(rpC+1:rpC+nM, rpC+1:rpC+nM)
  diag(A_new[monthly_idx, monthly_idx]) <- A_new2[1:nm]
  if(rRi) {
    if(rRi == 2L) stop("Cannot estimate unrestricted observation covariance matrix together with AR(1) serial correlation")
    diag(Q_new[monthly_idx, monthly_idx]) <- Q_new2[1:nm]
  } else {
    Q_new[monthly_idx, monthly_idx] <- diag(nm)
  }

  # ==========================================
  # C matrix update for monthly variables
  # MATLAB: nom = nom + y(1:nM,t)*Zsmooth(1:r,t+1)'
  #               - nanYt*(Zsmooth(rpC+1:rpC+nM,t+1)*Zsmooth(1:r,t+1)'+Vsmooth(rpC+1:rpC+nM,1:r,t+1))
  # ==========================================
  denom <- numeric(nm * r^2)
  nom <- matrix(0, nm, r)

  for (t in seq_len(TT)) {
    nmiss <- as.double(NW[t, snm])
    tmp <- t(Zsmooth[t, sr])
    tmp2 <- crossprod(tmp) + Vsmooth[sr, sr, t]
    dim(tmp2) <- NULL
    denom %+=% tcrossprod(tmp2, nmiss)
    nom %+=% (XW0[t, snm] %*% tmp)
    # Subtract idiosyncratic error contribution for monthly variables
    nom %-=% ((Zsmooth[t, monthly_idx] %*% tmp + Vsmooth[monthly_idx, sr, t]) * nmiss)
  }

  dim(denom) <- c(r, r, nm)
  dnkron[dnkron_ind] <- aperm.default(denom, c(1L, 3L, 2L))
  C_new <- C
  C_new[snm, sr] <- solve.default(dnkron, unattrib(nom))

  # ==========================================
  # C matrix update for quarterly variables
  # And quarterly AR1 parameter updates
  # MATLAB: for i=n-nQ+1:n
  #           nom = nom - nanYt*([1 2 3 2 1]*Zsmooth(i_idio_jQ,t+1)*Zsmooth(1:rC,t+1)'+
  #                              [1 2 3 2 1]*Vsmooth(i_idio_jQ,1:rC,t+1))
  # ==========================================
  V_0_new <- V_0
  weights <- c(1, 2, 3, 2, 1)

  for (i in (nm + 1L):n) {
    denom_q <- numeric(rCr^2)
    nom_q <- numeric(rCr)

    j <- i - nm  # Quarterly variable index (1 to nq)
    i_idio_jQ <- rpC + nm + (j - 1L) * 5L + 1:5  # 5 error lag indices for this quarterly var

    for (t in seq_len(TT)) {
      if(NW[t, i]) {
        denom_q %+=% (crossprod(Zsmooth[t, 1:rCr, drop = FALSE]) + Vsmooth[1:rCr, 1:rCr, t])
        nom_q %+=% (XW0[t, i] * Zsmooth[t, 1:rCr])
        # Subtract quarterly idiosyncratic contribution with temporal aggregation weights
        nom_q %-=% (weights %*% (tcrossprod(Zsmooth[t, i_idio_jQ], Zsmooth[t, 1:rCr]) +
                                  Vsmooth[i_idio_jQ, 1:rCr, t]))
      }
    }

    dim(denom_q) <- c(rCr, rCr)
    denom_inv <- ainv(denom_q)
    C_i <- denom_inv %*% nom_q
    # Apply Rcon constraint via restricted least squares
    tmp <- tcrossprod(denom_inv, R_mat)
    C_i_constr <- C_i - tmp %*% ainv(R_mat %*% tmp) %*% R_mat %*% C_i
    C_new[i, 1:rCr] <- C_i_constr

    # Update quarterly AR1 parameters
    # MATLAB: V_0(i_idio_jQ,i_idio_jQ) = Vsmooth(i_idio_jQ,i_idio_jQ,1)
    #         A_new(i_idio_jQ(1),i_idio_jQ(1)) = A_new2(i_idio_jQ(1)-rpC,i_idio_jQ(1)-rpC)
    #         Q_new(i_idio_jQ(1),i_idio_jQ(1)) = Q_new2(i_idio_jQ(1)-rpC,i_idio_jQ(1)-rpC)
    V_0_new[i_idio_jQ, i_idio_jQ] <- Vsmooth[i_idio_jQ, i_idio_jQ, 1]
    A_new[i_idio_jQ[1], i_idio_jQ[1]] <- A_new2[nm + j]
    Q_new[i_idio_jQ[1], i_idio_jQ[1]] <- Q_new2[nm + j]
  }

  # ==========================================
  # Update initial conditions
  # MATLAB: V_0(1:rpC,1:rpC) = Vsmooth(1:rpC,1:rpC,1)
  #         V_0(rpC+1:end,rpC+1:end) = diag(diag(Vsmooth(rpC+1:end,rpC+1:end,1)))
  # ==========================================
  V_0_new[srpC, srpC] <- Vsmooth0[srpC, srpC]
  V_0_new[monthly_idx, monthly_idx] <- diag(diag(Vsmooth0[monthly_idx, monthly_idx, drop = FALSE]))
  # Note: quarterly V_0_new blocks were updated inside the loop above

  return(list(
    A = A_new,
    C = C_new,
    Q = Q_new,
    R = R,  # R stays fixed at kappa * I
    F_0 = drop(Zsmooth0),
    P_0 = V_0_new,
    loglik = loglik
  ))
}
