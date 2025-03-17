
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

# Z_0 = F_0; V_0 = P_0
EMstepBMMQ = function(X, A, C, Q, R, Z_0, V_0, XW0, NW, dgind, dnkron, dnkron_ind, r, p, R_mat, n, nq, sr, TT, rQi, rRi) {

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

  # Compute expected sufficient statistics for a single Kalman filter sequence
  # Run the Kalman filter with the current parameter estimates
  kfs_res = SKFS(X, A, C, Q, R, Z_0, V_0, TRUE)

  # Extract values from the Kalman Filter output
  Zsmooth = kfs_res$F_smooth
  Vsmooth = kfs_res$P_smooth
  VVsmooth = kfs_res$PPm_smooth
  Zsmooth0 = kfs_res$F_smooth_0
  Vsmooth0 = kfs_res$P_smooth_0
  loglik = kfs_res$loglik

  # Perform algebraic operations involving E(Z_t), E(Z_{t-1}), Y_t, Y_{t-1} to update A_new, Q_new, C_new, and R_new

  # Normal Expectations to get Z(t+1) and A(t+1).
  tmp = rbind(Zsmooth0[srpC], Zsmooth[-TT, srpC, drop = FALSE])
  tmp2 = sum3(Vsmooth[srpC, srpC, -TT, drop = FALSE])
  EZZ = crossprod(Zsmooth[, srpC, drop = FALSE]) %+=% (tmp2 + Vsmooth[srpC, srpC, TT])                   # E(Z'Z)
  EZZ_BB = crossprod(tmp) %+=% (tmp2 + Vsmooth0[srpC, srpC])                                             # E(Z(-1)'Z(-1))
  EZZ_FB = crossprod(Zsmooth[, srpC, drop = FALSE], tmp) %+=% sum3(VVsmooth[srpC, srpC,, drop = FALSE])  # E(Z'Z(-1))

  # Update matrices A and Q
  A_new = A
  Q_new = Q

  # System matrices
  A_new[sr, srp] = EZZ_FB[sr, srp, drop = FALSE] %*% ainv(EZZ_BB[srp, srp, drop = FALSE])
  if(rQi) {
    Qsr = (EZZ[sr, sr] - tcrossprod(A_new[sr, srp, drop = FALSE], EZZ_FB[sr, srp, drop = FALSE])) / TT
    Q_new[sr, sr] = if(rQi == 2L) Qsr else diag(diag(Qsr))
  } else Q_new[sr, sr] = diag(r)
  Q_new[rpC1nq, rpC1nq] = diag(diag(crossprod(Zsmooth[, rpC1nq, drop = FALSE]) + sum3(Vsmooth[rpC1nq, rpC1nq,, drop = FALSE])) / TT)

  # E(X'X) & E(X'Z)
  # Estimate matrix C using maximum likelihood approach
  denom = numeric(nm*r^2)
  nom = matrix(0, nm, r)

  # nanYt = diagm(nanY[1:nM,t] .== 0)
  # denom += kron(Zsmooth[1:r,t+1] * Zsmooth[1:r,t+1]' + Vsmooth[1:r,1:r,t+1], nanYt)
  # nom += y[1:nM,t]*Zsmooth[1:r,t+1]';
  for (t in 1:TT) {
    tmp = Zsmooth[t, sr]
    nom %+=% tcrossprod(XW0[t, snm], tmp)
    tmp2 = tcrossprod(tmp) + Vsmooth[sr, sr, t]
    dim(tmp2) = NULL
    denom %+=% tcrossprod(tmp2, NW[t, snm])
  }

  dim(denom) = c(r, r, nm)
  dnkron[dnkron_ind] = aperm.default(denom, c(1L, 3L, 2L))
  # Solve for vec_C and reshape into C_new
  C_new = C
  C_new[1:nm, sr] = solve.default(dnkron, unattrib(nom)) # ainv() -> slower...

  # Now updating the quarterly observation matrix C
  for (i in (nm+1):n) {
      denom = numeric(rC^2)
      nom = numeric(rC)

      for (t in 1:TT) {
          if(NW[t, i]) denom %+=% (crossprod(Zsmooth[t, 1:rC, drop = FALSE]) + Vsmooth[1:rC, 1:rC , t])
          nom %+=% (XW0[t, i] * Zsmooth[t, 1:rC])
      }
      dim(denom) = c(rC, rC)
      denom_inv = ainv(denom)
      C_i = denom_inv %*% nom # Is this a scalar? -> yes, otherwise the replacement with C_new doesn't make sense
      tmp = tcrossprod(denom_inv, R_mat)
      C_i_constr = C_i - tmp %*% ainv(R_mat %*% tmp) %*% R_mat %*% C_i
      C_new[i, 1:rC] = C_i_constr
  }

  if(rRi) {
    R_new = matrix(0, n, n)
    Rdg = R[dgind]
    for (t in 1:TT) {
      nnanYt = NW[t, ]
      tmp = C_new * nnanYt
      R[dgind] = Rdg * !nnanYt # If R is not diagonal
      tmp2 = tmp %*% tcrossprod(Vsmooth[,, t], tmp)
      tmp2 %+=% tcrossprod(XW0[t, ] - tmp %*% Zsmooth[t, ])
      tmp2 %+=% R # If R is not diagonal
      # tmp2[dgind] = tmp2[dgind] + (Rdg * !nnanYt) # If R is diagonal...
      R_new %+=% tmp2
    }
    if(rRi == 2L) { # Unrestricted
      R_new %/=% TT
      R_new[R_new < 1e-7] = 1e-7
      R_new[(nm+1):n, ] = 0
      R_new[, (nm+1):n] = 0
    } else { # Diagonal
      RR = R_new[dgind] / TT
      RR[RR < 1e-7] = 1e-7 # RR(RR<1e-2) = 1e-2;
      R_new %*=% 0
      R_new[dgind[snm]] = RR[snm] # TODO: set quarterly obs to a small number?
    }
  } else {
    R_new = diag(n)
    diag(R_new)[(nm+1):n] = 0
  }

  # Set initial conditions
  V_0_new = V_0
  V_0_new[srpC, srpC] = Vsmooth0[srpC, srpC]
  diag(V_0_new)[-srpC] = diag(Vsmooth0)[-srpC]

  return(list(A = A_new,
              C = C_new,
              Q = Q_new,
              R = R_new,
              F_0 = drop(Zsmooth0),
              P_0 = V_0_new,
              loglik = kfs_res$loglik))

}
