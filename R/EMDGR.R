
EMstepDGR <- function(X, A, C, Q, R, F0, P0, cpX, n, r, sr, T, rQi, rRi) {

  ## E-step will return a list of sufficient statistics, namely second
  ## (cross)-moments for latent and observed data. This is then plugged back
  ## into M-step.
  list2env(Estep(X, C, Q, R, A, F0, P0), envir = environment())
  betasr <- beta[sr, , drop = FALSE]

  ## M-step computes model parameters as a function of the sufficient
  ## statistics that were computed with the E-step. Iterate the procedure
  ## until convergence. Due to the model specification, likelihood maximiation
  ## in the M-step is just an OLS estimation. In particular, X_t = C*F_t and
  ## F_t = A*F_(t-1).

  C[, sr] <- delta[, sr] %*% apinv(gamma[sr, sr, drop = FALSE])
  A_update <- betasr %*% ainv(gamma1)
  A[sr, ] <- A_update
  if(rQi) {
    Qsr <- (gamma2[sr, sr] - tcrossprod(A_update, betasr)) / (T-1L)
    Q[sr, sr] <- if(rQi == 2L) Qsr else diag(diag(Qsr))
  } else Q[sr, sr] <- diag(r)

  if(rRi) {
    R <- (cpX - tcrossprod(C, delta)) / T
    if(rRi == 2L) R[R < 1e-7] <- 1e-7 else {
      RR <- diag(R)
      RR[RR < 1e-7] <- 1e-7
      R <- diag(RR)
    }
  } else R <- diag(n)

  return(list(A = A, C = C, Q = Q, R = R, F0 = F0, P0 = P0, loglik = loglik))

}