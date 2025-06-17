# Define the inputs of the function
# X: a matrix with n rows and TT columns representing the observed data at each time t.
# A: a transition matrix with size r*p+r*r.
# C: an observation matrix with size n*r.
# Q: a covariance matrix for the state equation disturbances with size r*p+r*r.
# R: a covariance matrix for the observation disturbances with size n*n.
# Z_0: the initial value of the state variable.
# V_0: the initial value of the variance-covariance matrix of the state variable.
# r: order of the state autoregression (AR) process.
# p: number of lags in the state AR process.
# i_idio: a vector of indicators specifying which variables are idiosyncratic noise variables.'

# TODO: FINISH and VERIFY!!
# Z_0 = F_0; V_0 = P_0
EMstepBMMQidio = function(X, A, C, Q, R, Z_0, V_0, XW0, NW, dgind, dnkron, dnkron_ind, r, p, R_mat, n, nq, sr, TT, rQi, rRi) {

  # Define dimensions of the input arguments
  rC = 5L # ncol(Rcon)
  pC = max(p, rC)
  rp = r * p
  rpC = r * pC
  nm = n - nq
  rC = rC * r # Note here replacing rC!!

  snm = seq_len(nm)
  srp = seq_len(rp)
  srpC = seq_len(rpC)
  rpC1nq = (rpC+1L):(rpC+nq)

  # Run Kalman filter with current parameter estimates
  kfs_res = SKFS(X, A, C, Q, R, Z_0, V_0, TRUE)

  # Extract values from the Kalman Filter output
  Zsmooth = kfs_res$F_smooth
  Vsmooth = kfs_res$P_smooth
  VVsmooth = kfs_res$PPm_smooth
  Zsmooth0 = kfs_res$F_smooth_0
  Vsmooth0 = kfs_res$P_smooth_0
  loglik = kfs_res$loglik

  # Normal Expectations for factors
  tmp = rbind(Zsmooth0[srpC], Zsmooth[-TT, srpC, drop = FALSE])
  tmp2 = sum3(Vsmooth[srpC, srpC, -TT, drop = FALSE])
  EZZ = crossprod(Zsmooth[, srpC, drop = FALSE]) %+=% (tmp2 + Vsmooth[srpC, srpC, TT])                   # E(Z'Z)
  EZZ_BB = crossprod(tmp) %+=% (tmp2 + Vsmooth0[srpC, srpC])                                             # E(Z(-1)'Z(-1))
  EZZ_FB = crossprod(Zsmooth[, srpC, drop = FALSE], tmp) %+=% sum3(VVsmooth[srpC, srpC,, drop = FALSE])  # E(Z'Z(-1))

  # Expectations for idiosyncratic errors
  tmp = rbind(Zsmooth0[rpC1nq], Zsmooth[-TT, rpC1nq, drop = FALSE])
  tmp2 = sum3(Vsmooth[rpC1nq, rpC1nq, -TT, drop = FALSE])
  EZZ_u = diag(crossprod(Zsmooth[, rpC1nq, drop = FALSE])) + diag(tmp2 + Vsmooth[rpC1nq, rpC1nq, TT])                     # E(Z'Z)
  EZZ_BB_u = diag(crossprod(tmp)) + diag(tmp2 + Vsmooth0[rpC1nq, rpC1nq])                                                # E(Z(-1)'Z(-1))
  EZZ_FB_u = diag(diag(crossprod(Zsmooth[, rpC1nq, drop = FALSE], tmp)) + diag(rowSums(VVsmooth[rpC1nq, rpC1nq,, drop = FALSE], dims = 2L)))  # E(Z'Z(-1))

  # Update matrices A and Q
  A_new = A
  Q_new = Q

  # System matrices for factors
  A_new[sr, srp] = EZZ_FB[sr, , drop = FALSE] %*% ainv(EZZ_BB)
  if(rQi) {
    Qsr = (EZZ[sr, sr] - tcrossprod(A_new[sr, srp, drop = FALSE], EZZ_FB[sr,, drop = FALSE])) / TT
    Q_new[sr, sr] = if(rQi == 2L) Qsr else diag(diag(Qsr))
  } else Q_new[sr, sr] = diag(r)

  # System matrices for errors
  A_new[rpC1nq, rpC1nq] = EZZ_FB_u / EZZ_BB_u
  if(rRi) {
    if(rRi == 2L) stop("Cannot estimate unrestricted observation covariance matrix together with AR(1) serial correlation")
    Q_new[rpC1nq, rpC1nq] = (diag(EZZ_u) - A_new[rpC1nq, rpC1nq] * EZZ_FB_u) / TT
  } else Q_new[rpC1nq, rpC1nq] = diag(nm)

  # Estimate matrix C using maximum likelihood approach
  denom = numeric(nm*r^2)
  nom = matrix(0, nm, r)

  for (t in seq_len(TT)) {
    nmiss = as.double(NW[t, 1:nm])
    tmp = t(Zsmooth[t, sr])
    tmp2 = crossprod(tmp) + Vsmooth[sr, sr, t]
    dim(tmp2) = NULL
    denom %+=% tcrossprod(tmp2, nmiss)
    nom %+=% (XW0[t, 1:nm] %*% tmp)
    nom %-=% ((Zsmooth[t, rpC1nq] %*% tmp + Vsmooth[rpC1nq, sr, t]) * nmiss)
  }

  dim(denom) = c(r, r, nm)
  dnkron[dnkron_ind] = aperm.default(denom, c(1L, 3L, 2L))
  vec_C = solve.default(dnkron, unattrib(nom))
  C_new = C
  C_new[1:nm, sr] = vec_C

  # Now updating the quarterly observation matrix C
  for (i in (nm+1):n) {
    denom = numeric(rC^2)
    nom = numeric(rC)
    i_idio_jQ = rpC + nm + 5*((i-nm)-1) + 1:5

    for (t in 1:TT) {
      if(NW[t, i]) {
        denom %+=% (crossprod(Zsmooth[t, 1:rC, drop = FALSE]) + Vsmooth[1:rC, 1:rC, t])
        nom %+=% (XW0[t, i] * Zsmooth[t, 1:rC])
        nom %-=% (c(1,2,3,2,1) %*% (Zsmooth[t, i_idio_jQ] %*% t(Zsmooth[t, 1:rC]) + Vsmooth[i_idio_jQ, 1:rC, t]))
      }
    }

    dim(denom) = c(rC, rC)
    denom_inv = ainv(denom)
    C_i = denom_inv %*% nom
    tmp = tcrossprod(denom_inv, R_mat)
    C_i_constr = C_i - tmp %*% ainv(R_mat %*% tmp) %*% R_mat %*% C_i
    C_new[i, 1:rC] = C_i_constr

    # Update AR parameters for quarterly idiosyncratic errors
    V_0_new = V_0
    V_0_new[i_idio_jQ, i_idio_jQ] = Vsmooth[i_idio_jQ, i_idio_jQ, 1]
    A_new[i_idio_jQ[1], i_idio_jQ[1]] = A_new[rpC1nq[i-nm], rpC1nq[i-nm]]
    Q_new[i_idio_jQ[1], i_idio_jQ[1]] = Q_new[rpC1nq[i-nm], rpC1nq[i-nm]]
  }

  # Set initial conditions
  V_0_new = V_0
  V_0_new[srpC, srpC] = Vsmooth0[srpC, srpC]
  V_0_new[rpC1nq, rpC1nq] = diag(diag(Vsmooth0[rpC1nq, rpC1nq, drop = FALSE]))

  return(list(A = A_new,
              C = C_new,
              Q = Q_new,
              R = R,  # R stays fixed
              F_0 = drop(Zsmooth0),
              P_0 = V_0_new,
              loglik = loglik))
}
