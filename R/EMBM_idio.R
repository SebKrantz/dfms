
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
EMstepBMidio = function(X, A, C, Q, R, Z_0, V_0, XW0, W, dgind, dnkron, dnkron_ind, r, p, n, sr, TT, rQi, rRi) {

  # Define dimensions of the input arguments
  srp = seq_len(r*p)
  rp1nr = (r*p+1L):ncol(A)

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
  tmp = rbind(Zsmooth0[srp], Zsmooth[-TT, srp, drop = FALSE])
  tmp2 = sum3(Vsmooth[srp, srp, -TT, drop = FALSE])
  EZZ = crossprod(Zsmooth[, srp, drop = FALSE]) %+=% (tmp2 + Vsmooth[srp, srp, TT])                   # E(Z'Z)
  EZZ_BB = crossprod(tmp) %+=% (tmp2 + Vsmooth0[srp, srp])                                            # E(Z(-1)'Z(-1))
  EZZ_FB = crossprod(Zsmooth[, srp, drop = FALSE], tmp) %+=% sum3(VVsmooth[srp, srp,, drop = FALSE])  # E(Z'Z(-1))

  # Expectations on the erros (u_t).
  tmp = rbind(Zsmooth0[rp1nr], Zsmooth[-TT, rp1nr, drop = FALSE])
  tmp2 = sum3(Vsmooth[rp1nr, rp1nr, -TT, drop = FALSE])
  EZZ_u = diag(crossprod(Zsmooth[, rp1nr, drop = FALSE])) + diag(tmp2 + Vsmooth[rp1nr, rp1nr, TT])                     # E(Z'Z)
  EZZ_BB_u = diag(crossprod(tmp)) + diag(tmp2 + Vsmooth0[rp1nr, rp1nr])                                                # E(Z(-1)'Z(-1))
  EZZ_FB_u = diag(diag(crossprod(Zsmooth[, rp1nr, drop = FALSE], tmp)) + diag(rowSums(VVsmooth[rp1nr, rp1nr,, drop = FALSE], dims = 2L)))  # E(Z'Z(-1))

  # Update matrices A and Q
  A_new = A
  Q_new = Q

  # System matrices
  A_new[sr, srp] = EZZ_FB[sr, , drop = FALSE] %*% ainv(EZZ_BB)
  if(rQi) {
    Qsr = (EZZ[sr, sr] - tcrossprod(A_new[sr, srp, drop = FALSE], EZZ_FB[sr,, drop = FALSE])) / TT
    Q_new[sr, sr] = if(rQi == 2L) Qsr else diag(diag(Qsr))
  } else Q_new[sr, sr] = diag(r)

  # Errors
  A_new[rp1nr, rp1nr] = EZZ_FB_u / EZZ_BB_u  # %r*% (1/EZZ_BB_u) # Same as EZZ_FB_u is diagonal...
  if(rRi) {
    if(rRi == 2L) stop("Cannot estimate unrestricted observation covariance matrix together with AR(1) serial correlation")
    Q_new[rp1nr, rp1nr] = (diag(EZZ_u) - A_new[rp1nr, rp1nr] * EZZ_FB_u) / TT # tcrossprod(A_new[rp1nr, rp1nr, drop = FALSE], EZZ_FB_u) # Same as EZZ_FB_u is diagonal...
  } else Q_new[rp1nr, rp1nr] = diag(n)

  # E(X'X) & E(X'Z)
  # Estimate matrix C using maximum likelihood approach
  denom = numeric(n*r^2)
  nom = matrix(0, n, r)

  # nanYt = diag(~nanY(:,t));
  # denom = denom + kron(Zsmooth(1:r,t+1)*Zsmooth(1:r,t+1)'+Vsmooth(1:r,1:r,t+1),nanYt);
  #   nom = nom + y(:,t)*Zsmooth(1:r,t+1)'-nanYt(:,i_idio)*(Zsmooth(r*p+1:end,t+1)*Zsmooth(1:r,t+1)'+Vsmooth(r*p+1:end,1:r,t+1));

  for (t in seq_len(TT)) {
    # select non-missing columns of Y for time t
    nmiss = as.double(!W[t, ])
    tmp = t(Zsmooth[t, sr])
    # Add components to denom (same as EMBM.R)
    tmp2 = crossprod(tmp) + Vsmooth[sr, sr, t]
    dim(tmp2) = NULL
    denom %+=% tcrossprod(tmp2, nmiss)
    # add components to nom
    nom %+=% (XW0[t, ] %*% tmp)
    nom %-=% ((Zsmooth[t, rp1nr] %*% tmp + Vsmooth[rp1nr, sr, t]) * nmiss)
  }

  dim(denom) = c(r, r, n)
  dnkron[dnkron_ind] = aperm.default(denom, c(1L, 3L, 2L))
  # Solve for vec_C and reshape into C_new
  vec_C = solve.default(dnkron, unattrib(nom)) # ainv() -> slower...
  C_new = C
  C_new[, sr] = vec_C

  # R_new = zeros(n,n);
  # for t=1:T
  #    nanYt = diag(~nanY(:,t));
  #    R_new = R_new + (y(:,t)-nanYt*C_new*Zsmooth(:,t+1))*(y(:,t)-nanYt*C_new*Zsmooth(:,t+1))'+
  #             nanYt*C_new*Vsmooth(:,:,t+1)*C_new'*nanYt+(eye(n)-nanYt)*R*(eye(n)-nanYt);
  # end
  # R_new = R_new/T;
  # R_new = diag(diag(R_new));

  # # Needed ??? -> Nope, gives same result as fixed R
  # R_new = matrix(0, n, n)
  # Rdg = R[dgind]
  #
  # for (t in seq_len(TT)) {
  #   nanYt = W[t, ]
  #   tmp = C_new * !nanYt
  #   R[dgind] = Rdg * nanYt # If R is not diagonal
  #   tmp2 = tmp %*% tcrossprod(Vsmooth[,, t], tmp)
  #   tmp2 %+=% tcrossprod(XW0[t, ] - tmp %*% Zsmooth[t, ])
  #   tmp2 %+=% R # If R is not diagonal
  #   # tmp2[dgind] = tmp2[dgind] + (Rdg * nanYt) # If R is diagonal...
  #   R_new %+=% tmp2
  # }
  #
  # RR = R_new[dgind] / TT
  # RR[RR < 1e-7] = 1e-7 # RR(RR<1e-2) = 1e-2;
  # R_new = diag(RR)

  # Set initial conditions
  V_0_new = V_0
  V_0_new[srp, srp] = Vsmooth0[srp, srp]
  V_0_new[rp1nr, rp1nr] = diag(diag(Vsmooth0[rp1nr, rp1nr, drop = FALSE]))

  return(list(A = A_new,
              C = C_new,
              Q = Q_new,
              R = R, # R_new,
              F_0 = drop(Zsmooth0),
              P_0 = V_0_new,
              loglik = kfs_res$loglik))

}
