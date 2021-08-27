#' Estimate a Dynamic Factor Model
#'
#' Efficient estimation of a Dynamic Factor Model via the EM Algorithm - on stationary data of a single frequency,
#' with time-invariant system matrices and classical assumptions, while permitting arbitrary patterns of missing data.
#'
#' @param X data matrix or frame.
#' @param r number of factors.
#' @param p number of lags in factor VAR.
#' @param rQ restrictions on the state (transition) covariance matrix (Q).
#' @param rR restrictions on the observation (measurement) covariance matrix (R).
#' @param max.inter maximum number of EM iterations.
#' @param tol EM convergence tolerance.
#' @param miss.threshold proportions of series missing for a case to be removed.
#'
#' @details
#' This function efficiently estimates a Dynamic Factor Model with the following classical assumptions:
#' \enumerate{
#' \item Linearity
#' \item Idiosynchratic measurement (observation) errors
#' \item No relationship between series and lagged factors (\emph{ceteris paribus})
#' \item No relationship between lagged error terms in the either measurement or transition equation (no serial correlation)
#' }
#' Factors are allowed to evolve in a \eqn{VAR(p)} process, and data is standardized (scaled and centered) before estimation (removing the need of intercept terms).
#' By assumptions 1-4, this translates into the following dynamic form:
#'
#' \deqn{\textbf{x}_t = \textbf{C}_0 \textbf{f}_t + \textbf{e}_t \tilde N(\textbf{0}, \textbf{R})}{xt = C0 ft + et ~ N(0, R)}
#' \deqn{\textbf{f}_t = \sum_{i=1}^p \textbf{A}_p \textbf{f}_{t-p} + \textbf{u}_t \tilde N(0, \textbf{Q}_0)}{ft = A1 ft-1 + \dots + Ap ft-p + ut ~ N(0, Q0)}
#'
#' where the first equation is called the measurement or observation equation and the second equation is called transition, state or process equation, and
#'
#' \tabular{llll}{
#'  \eqn{n} \tab\tab number of series in \eqn{\textbf{x}_t}{xt} (\eqn{r} and \eqn{p} as the arguments to \code{DFM}).\cr\cr
#'  \eqn{\textbf{x}_t}{xt} \tab\tab \eqn{n \times 1}{n x 1} vector of observed series at time \eqn{t}{t}: \eqn{(x_{1t}, \dots, x_{nt})'}{(x1t, \dots, xnt)'}. Some observations can be missing.  \cr\cr
#'  \eqn{\textbf{f}_t}{ft} \tab\tab \eqn{r \times 1}{r x 1} vector of factors at time \eqn{t}{t}: \eqn{(f_{1t}, \dots, f_{rt})'}{(f1t, \dots, frt)'}.\cr\cr
#'  \eqn{\textbf{C}_0}{C0} \tab\tab \eqn{n \times r}{n x r} measurement (observation) matrix.\cr\cr
#'  \eqn{\textbf{A}_j}{Aj} \tab\tab \eqn{r \times r}{r x r} state transition matrix at lag \eqn{j}{j}. \cr\cr
#'  \eqn{\textbf{Q}_0}{Q0} \tab\tab \eqn{r \times r}{r x r} state covariance matrix.\cr\cr
#'  \eqn{\textbf{R}}{R} \tab\tab \eqn{n \times n}{n x n} measurement (observation) covariance matrix. It is diagonal by assumption 2 that \eqn{E[\textbf{x}_{it}|\textbf{x}_{-i,t},\textbf{x}_{i,t-1}, \dots, \textbf{f}_t, \textbf{f}_{t-1}, \dots] = \textbf{Cf}_t \forall i}{E[xit|x(-i)t, xt-1, \dots, ft, ft-1, \dots] = C ft}. This assumption is also referred to as 'Exact DFM' by Stock & Watson (2016), where all correlation between the series is accounted for by the latent factors.\cr\cr
#' }
#'
#'
#' This model can be estimated using a classical form of the Kalman Filter and the Expectation Maximization (EM) algorithm, after transforming it to State-Space (stacked, VAR(1)) form:
#'
#' \deqn{\textbf{x}_t = \textbf{C} \textbf{F}_t + \textbf{e}_t \tilde N(\textbf{0}, \textbf{R})}{xt = C Ft + et ~ N(0, R)}
#' \deqn{\textbf{F}_t = \textbf{A F}_{t-1} + \textbf{u}_t \tilde N(0, \textbf{Q})}{Ft = A Ft-1 + ut ~ N(0, Q)}
#'
#' where
#' \tabular{llll}{
#'  \eqn{n} \tab\tab number of series in \eqn{\textbf{x}_t}{xt} (\eqn{r} and \eqn{p} as the arguments to \code{DFM}).\cr\cr
#'  \eqn{\textbf{x}_t}{xt} \tab\tab \eqn{n \times 1}{n x 1} vector of observed series at time \eqn{t}{t}: \eqn{(x_{1t}, \dots, x_{nt})'}{(x1t, \dots, xnt)'}. Some observations can be missing.  \cr\cr
#'  \eqn{\textbf{F}_t}{Ft} \tab\tab \eqn{rp \times 1}{rp x 1} vector of stacked factors at time \eqn{t}{t}: \eqn{(f_{1t}, \dots, f_{rt}, f_{1,t-1}, \dots, f_{r,t-1}, \dots, f_{1,t-p}, \dots, f_{r,t-p})'}{(f1t, \dots, frt, f1t-1, \dots, frt-1, \dots, f1t-p, \dots, frt-p)'}.\cr\cr
#'  \eqn{\textbf{C}}{C} \tab\tab \eqn{n \times rp}{n x rp} observation matrix. Only the first \eqn{n \times r}{n x r} terms are non-zero, by assumption 3 that \eqn{E[\textbf{x}_t|\textbf{F}_t] = E[\textbf{x}_t|\textbf{f}_t]}{E[Xt|Ft] = E[Xt|ft]} (no relationship of observed series with lagged factors given contemporaneous factors).\cr\cr
#'  \eqn{\textbf{A}}{A} \tab\tab stacked \eqn{rp \times rp}{rp x rp} state transition matrix consisting of 3 parts: the top \eqn{r \times rp}{r x rp} part provides the dynamic relationships captured by \eqn{(\textbf{A}_1, \dots, \textbf{A}_p)}{(A1, \dots, Ap)} in the dynamic form, the terms \code{A[(r+1):rp, 1:(rp-r)]} constitute an \eqn{(rp-r) \times (rp-r)}{(rp-r) x (rp-r)} identity matrix mapping all lagged factors to their known values at times t. The remining part \code{A[(rp-r+1):rp, (rp-r+1):rp]} is an \eqn{r \times r}{r x r} matrix of zeros. \cr\cr
#'  \eqn{\textbf{Q}}{Q} \tab\tab \eqn{rp \times rp}{rp x rp} state covariance matrix. The top \eqn{r \times r}{r x r} part gives the contemporaneous relationships, the rest are zeros by assumption 4.\cr\cr % that \eqn{E[\textbf{f}_t|\textbf{F}_{t-1}] = E[\textbf{f}_t|\textbf{f}_{t-1}] = \textbf{A}_1 \textbf{f}_{t-1}}{E[ft|Ft-1] = E[ft|ft-1] = A1 ft-1} (all relationships between lagged factors are captured in \eqn{\textbf{A}_1}{A1}).\cr\cr
#'  \eqn{\textbf{R}}{R} \tab\tab \eqn{n \times n}{n x n} observation covariance matrix. It is diagonal by assumption 2 and identical to \eqn{\textbf{R}}{R} as stated in the dynamic form.\cr\cr
#' }
#' @useDynLib DFM, .registration = TRUE
#' @export

DFM <- function(X, r, p = 1L,
                rQ = c("none", "diagonal", "identity"),
                rR = c("diagonal", "identity"),
                max.iter = 100L, tol = 1e-4,
                miss.threshold = 0.8) {

  rp = r * p
  sr <- 1:r
  srp <- 1:rp
  X = fscale(qM(X))
  T = dim(X)[1L]
  n = dim(X)[2L]

  # Missing values
  X_imp = X
  anymiss = anyNA(X)
  if(anymiss) { # Simple median imputation
    W = is.na(X)
    X_imp[W] = fmedian(X, TRA = 1L)[W]
  }

  # Run PCA to get initial factor estimates:
  v = svd(X_imp, nu = 0L, nv = min(as.integer(r), n, T))$v
  F_pc = X_imp %*% v

  # Observation equation -------------------------------
  # Static predictions (all.equal(unattrib(HDB(X_imp, F_pc)), unattrib(F_pc %*% t(v))))
  C <- cbind(v, matrix(0, n, rp-r))
  res = X_imp - F_pc %*% t(v) # residuals from static predictions
  if(anymiss) res[W] = NA
  R = diag(fvar(res)) # Covariance (assumed idiosynchratic)
  # R = cov(res) # unrestricted covariance estimate

  # Transition equation -------------------------------
  var = VAR(F_pc, p)
  A = rbind(t(var$A), diag(1, rp-r, rp)) # var$A is rp x r matrix
  Q = matrix(0, rp, rp)
  Q[sr, sr] = cov(var$res) # unrestricted covariance estimate

  # Initial state and state covariance (P) ------------
  x0 = var$X[1L, ] # rep(0, rp) # This should better be called f0, the factors are the state
  # Kalman gain is normally A %*% t(A) + Q, but here A is somewhat tricky...
  P0 = matrix(ginv(kronecker(A, A)) %*% unattrib(Q), rp, rp)
  # BM2014:
  # P0 = matrix(solve(diag(rp^2) - kronecker(A, A)) %*% unattrib(Q), rp, rp)

  # FKF::fkf(x0, P0, x0, rep(0, n), A, C, Q, R, X)

  ## Run standartized data through Kalman filter and smoother once
  kf_res <- KalmanFilter(X, C, Q, R, A, x0, P0)
  ks_res <- with(kf_res, KalmanSmoother(A, C, R, xF, xP, Pf, Pp))

  ## Two-step solution is state mean from the Kalman smoother
  F_kal <- ks_res$xS
  # return(F_kal[, sr])

  previous_loglik <- -.Machine$double.xmax
  num_iter <- 0L
  converged <- FALSE

  xx <- if(anymiss) na_omit(X) else X
  while ((num_iter < max.iter) & !converged) {

    ## E-step will return a list of sufficient statistics, namely second
    ## (cross)-moments for latent and observed data. This is then plugged back
    ## into M-step.
    em_res <- Estep(X, C, Q, R, A, x0, P0)
    beta <- em_res$beta_t
    gamma <- em_res$gamma_t
    delta <- em_res$delta_t
    gamma1 <- em_res$gamma1_t
    gamma2 <- em_res$gamma2_t
    P1sum <- em_res$V1 + tcrossprod(em_res$x1)
    x1sum <- em_res$x1
    loglik <- em_res$loglik_t

    num_iter <- num_iter + 1L

    ## M-step computes model parameters as a function of the sufficient
    ## statistics that were computed with the E-step. Iterate the procedure
    ## until convergence. Due to the model specification, likelihood maximiation
    ## in the M-step is just an OLS estimation. In particular, X_t = C*F_t and
    ## F_t = A*F_(t-1).

    C[, sr] <- delta[, sr] %*% ginv(gamma[sr, sr])

    A_update <- beta[sr, srp, drop = FALSE] %*% solve(gamma1[srp, srp])
    A[sr, srp] <- A_update
    Q[sr, sr] <- (gamma2[sr, sr] - tcrossprod(A_update, beta[sr, srp, drop = FALSE])) / (T-1)

    R <- (crossprod(xx) - tcrossprod(C, delta)) / T
    RR <- diag(R); RR[RR < 1e-7] <- 1e-7; R <- diag(RR)
    R <- diag(diag(R))

    ## Assign new initial values for next EM-algorithm step
    x0 <- x1sum
    P0 <- P1sum - tcrossprod(x0)

    converged <- em_converged(loglik, previous_loglik, threshold = tol)
    previous_loglik <- loglik

    ## Iterate at least 25 times
    if(num_iter < 25L) converged <- FALSE
  }

  if(converged) message("Converged after ", num_iter, " iterations.")
  else warning("Maximum number of iterations reached.")

  ## Run the Kalman filtering and smoothing step for the last time
  ## with optimal estimates
  kf <- KalmanFilter(X, C, Q, R, A, x0, P0)
  F_hat <- KalmanSmoother(A, C, R, kf$xF, kf$xP, kf$Pf, kf$Pp)$xS
  final_object <- list(pca = F_pc,
                       twostep = F_kal[, sr],
                       qml = F_hat[, sr],
                       A = A[sr, ],
                       C = C[, sr],
                       Q = Q[sr, sr],
                       R = R)

  class(final_object) <- "dfm"
  final_object
}
