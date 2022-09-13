library(DFM)
library(collapse)


dynFA_OPT <- function(X, r, p, max_iter, thresh = 1e-4) { # [F_hat,F_pc,F_kal,num_iter]
  ########################################################################################################
  # "A Quasi?Maximum Likelihood Approach for Large, Approximate Dynamic Factor Models,"
  # The Review of Economics and Statistics, MIT Press, vol. 94(4), pages 1014-1024, November 2012.
  # Catherine Doz, Universite' Cergy-Pontoise
  # Domenico Giannone, Universite' Libre de Bruxelles, ECARES and CEPR
  # Lucrezia Reichlin, London Business School and CEPR
  #
  #
  # Programs are also available at: http://homepages.ulb.ac.be/~dgiannon/
  #
  ########################################################################################################
  #
  #
  # DynFA:            Extracts the unobservable factors using three different methods
  #
  #       - QML:      Max Likelihood estimates using the Expectation Maximization (EM) algorithm
  #                    (Doz, Giannone and Reichlin, 2012)
  #
  #       - TWO STEP: Principal components + Kalman filtering
  #                   Doz, Catherine & Giannone, Domenico & Reichlin, Lucrezia, 2011.
  #                   "A two-step estimator for large approximate dynamic factor models based on Kalman filtering,"
  #                   Journal of Econometrics, Elsevier, vol. 164(1), pages 188-205, September.
  #
  #       - PC:       principal components
  #
  # INPUTS
  # X - matrix of observable variables
  # r - # of static factors
  # q - # of dynamic factors
  # p - # length of ar filter on common factors
  # max_iter - max # of iterations in estimation
  #
  # OUTPUTS
  # F_hat -   factors from QML
  # F_pc  -   factors using principal components
  # F_kal -   factors using from two steps

  X = qM(X)
  .c(T, N) %=% dim(X)

  # x is the matrix of standardized observable variables
  x = fscale(X)

  # the number of static factors cannot be great of the dynamic ones
  q = r
  # if (r < q) stop('q has to be less or equal to r')

  nlag = p-1
  rp = r * p

  #extract the first r eigenvectors and eigenvalues from cov(x)
  eigen_decomp = eigen(cov(x), symmetric = TRUE)
  v = eigen_decomp$vectors[, 1:r, drop = FALSE]
  # d = eigen_decomp$values[1:r]
  chi = x %*% tcrossprod(v)         # common component
  F = x %*% v
  F_pc = F         # factors using principal components
  var = .VAR(F, p) # run the var chi(t) = A*chi(t-1) + e(t);
  A = rbind(t(var$A), diag(1, rp-r, rp))
  Q = matrix(0, rp, rp)
  Q[1:r, 1:r] = cov(var$res) # VAR covariance matrix
  R = diag(diag(cov(x-chi))) # R diagonal
  initx = var$X[1, ]         # Initial state mean

  # initial state covariance
  initV = apinv(diag(rp^2)-kronecker(A,A)) %*% unattrib(Q)
  dim(initV) = c(rp, rp)

  C = cbind(v, matrix(0, N, r*nlag))

  # initialize the estimation and define ueful quantities
  previous_loglik = -Inf
  loglik = 0
  num_iter = 0
  y = t(x)
  LL = -Inf
  converged = FALSE

  # estimation of the factors with the Kalman filter using as initial values
  # for A, C, Q, R, initx, initV the ones computed with the principal
  # components
  kfs_res = SKFS(x, A, C, Q, R, initx, initV)
  xsmooth = kfs_res$F_smooth
  Vsmooth = kfs_res$P_smooth
  VVsmooth = kfs_res$PPm_smooth
  F_kal =  xsmooth[, 1:r]

  # factors estimation with the EM algorithm

  #repeat the algorithm until convergence
  while ((num_iter < max_iter) & !converged) {

    ### E step : compute the sufficient statistics

    # In our model the sufficient statistics are
    # delta = sum_t=1^T (x_t * f'_t)
    # gamma = sum_t=1^T (f_t * f'_t)
    # gamma1 = sum_t=2^T (f_t-1 * f'_t-1)
    # gamma2 = sum_t=2^T (f_t * f'_t)
    # beta = sum_t=1^T (f_t * f'_t-1)
    # P1sum  variance of the initial state
    # x1sum  expected value of the initial state

    # initialize to zero all the sufficient statistics
    delta = matrix(0, N, rp)
    gamma = matrix(0, rp, rp)
    gamma1 = matrix(0, rp, rp)
    gamma2 = matrix(0, rp, rp)
    beta = matrix(0, rp, rp)
    P1sum = matrix(0, rp, rp)
    x1sum = numeric(rp)
    loglik = 0

    # use the function Estep to update the expected sufficient statistics
    # note that at the first iteration  we use as initial values for A, C, Q,
    # R, initx, initV the ones computed with the principal components
    .c(beta_t, gamma_t, delta_t, gamma1_t, gamma2_t, x1, V1, loglik_t, xsmooth) %=% Estep_OPT(y, A, C, Q, R, initx, initV)

    # fix the expected sufficient statistics equal to the one computed with
    # the function Estep
    beta %+=% beta_t
    gamma %+=% gamma_t
    delta %+=% delta_t
    gamma1 %+=% gamma1_t
    gamma2 %+=% gamma2_t
    P1sum %+=% (V1 + tcrossprod(x1))
    x1sum = x1sum + x1

    # update the loglikelihood
    loglik = loglik + loglik_t

    # update the counter for the iterations
    num_iter =  num_iter + 1L

    ### M step
    # compute the parameters of the model as a function of the sufficient
    # statistics (computed with the function Estep)

    # The formulas for the parameters derive from the maximum likelihood
    # method. In the EM algorithm we substitute in the ML estimator the
    # sufficient statistics (computed in the E step) and then iterate the
    # procedure till the maximization of the likelihood

    # C = (sum_t=1^T x_t*f'_t)* (sum_t=1^T f_t*f'_t)^-1
    # substituting for the sufficient statistics

    C[, 1:r] = delta[, 1:r] %*% apinv(gamma[1:r, 1:r])

    # A = (sum_t=2^T f_t*f'_t-1)* (sum_2=1^T f_t-1*f'_t-1)^-1
    Atemp = beta[1:r, 1:rp] %*% ainv(gamma1[1:rp, 1:rp])
    A[1:r, 1:rp] = Atemp

    # Q = ( (sum_t=2^T f_t*f'_t) - A * (sum_2=1^T f_t-1*f'_t) )/(T-1)
    Q[1:r, 1:r] = (gamma2[1:r, 1:r] - Atemp %*% t(beta[1:r, 1:rp])) / (T-1)

    # R = ( sum_t=1^T (x_t*x'_t) - C * f_t*x'_t) )/T
    R = (crossprod(x) - C %*% t(delta))/T

    RR = diag(R); RR[RR < 1e-7] = 1e-7; R = diag(RR)

    R = diag(diag(R))                  # R diagonal

    loglik = loglik_t

    LL = c(LL, loglik)

    initx = x1sum
    initV = P1sum - tcrossprod(initx)

    converged = em_converged(loglik, previous_loglik, thresh)[1L]

    previous_loglik = loglik
  }

  # Final run of the Kamlman Filter and Smoother
  F_hat = SKFS(x, A, C, Q, R, initx, initV)$F_smooth[, 1:r]

  return(list(F_hat = F_hat, F_pc = F_pc, F_kal = F_kal, num_iter = num_iter))

}


#############################################################
# Estep
#############################################################


Estep_OPT <- function(y, A, C, Q, R, initx, initV) { # [beta, gamma, delta, gamma1, gamma2, x1, V1, loglik_t, xsmooth]

  # This function computes the (expected) sufficient statistics for a single Kalman filter sequence.
  #
  # y is the observable and x the hidden state

  # INPUTS
  # y(:,t) - the observation at time t
  # A - the system matrix
  # C - the observation matrix
  # Q - the system covariance
  # R - the observation covariance
  # initx - the initial state (column) vector
  # initV - the initial state covariance

  # OUTPUTS: the expected sufficient statistics, i.e.
  # beta = sum_t=1^T (x_t * x'_t-1)
  # gamma = sum_t=1^T (x_t * x'_t)
  # delta = sum_t=1^T (y_t * x'_t)
  # gamma1 = sum_t=2^T (x_t-1 * x'_t-1)
  # gamma2 = sum_t=2^T (x_t * x'_t)
  # x1  expected value of the initial state
  # V1  variance of the initial state
  # loglik value of the loglikelihood
  # xsmooth expected value of the state

  .c(N, T) %=% dim(y)
  rp = ncol(A)

  #   xsmooth = y;
  #   Vsmooth = zeros(rp, rp, T); # no uncertainty about the hidden states
  #   VVsmooth = zeros(rp, rp, T);
  #   loglik = 0;

  # use the Kalman smoother to compute
  # xsmooth = E[X(:,t) | y(:,1:T)]
  # Vsmooth = Cov[X(:,t) | y(:,1:T)]
  # VVsmooth = Cov[X(:,t), X(:,t-1) | y(:,1:T)] t >= 2
  # loglik = sum{t=1}^T log P(y(:,t))

  kfs_res = SKFS(t(y), A, C, Q, R, initx, initV, TRUE)
  loglik_t = kfs_res$loglik
  xsmooth = kfs_res$F_smooth
  Vsmooth = kfs_res$P_smooth
  VVsmooth = kfs_res$PPm_smooth

  # compute the expected sufficient statistics
  delta = matrix(0, N, rp)
  gamma = matrix(0, rp, rp)
  beta = matrix(0, rp, rp)

  for (t in 1:T) {
    delta %+=% tcrossprod(y[, t], xsmooth[t, ])
    gamma %+=% (tcrossprod(xsmooth[t, ]) + Vsmooth[,, t])
    if(t > 1) beta %+=% (tcrossprod(xsmooth[t, ], xsmooth[t-1, ]) + VVsmooth[,, t])
  }

  gamma1 = gamma - tcrossprod(xsmooth[T, ]) - Vsmooth[,, T]
  gamma2 = gamma - tcrossprod(xsmooth[1, ]) - Vsmooth[,, 1]

  x1 = xsmooth[1, ]
  V1 = Vsmooth[,, 1]

  return(list(beta = beta,
              gamma = gamma,
              delta = delta,
              gamma1 = gamma1,
              gamma2 = gamma2,
              x1 = x1, V1 = V1,
              loglik_t = loglik_t,
              xsmooth = xsmooth))
}



