# Quoting some functions that need to be evaluated iteratively
.EM_DGR <- quote(EMstepDGR(X, A, C, Q, R, F0, P0, cpX, n, r, sr, T, rQi, rRi))
.EM_BM <- quote(EMstepBM(X, A, C, Q, R, F0, P0))
.KFS <- quote(KalmanFilterSmoother(X, C, Q, R, A, F0, P0))


#' Estimate a Dynamic Factor Model
#'
#' Efficient estimation of a Dynamic Factor Model via the EM Algorithm - on stationary data of a single frequency,
#' with time-invariant system matrices and classical assumptions, while permitting arbitrary patterns of missing data.
#'
#' @param X data matrix or frame.
#' @param r number of factors.
#' @param p number of lags in factor VAR.
#' @param \dots further arguments to be added here in the future, such as further estimation methods or block-structures.
#' @param rQ restrictions on the state (transition) covariance matrix (Q).
#' @param rR restrictions on the observation (measurement) covariance matrix (R).
#' @param min.inter minimum number of EM iterations (to ensure a convergence path).
#' @param max.inter maximum number of EM iterations.
#' @param tol EM convergence tolerance.
#' @param max.missing proportion of series missing for a case to be considered missing.
#' @param na.rm.method method to apply concerning missing cases selected through \code{max.missing}: \code{"LE"} only removes cases at the beginning or end of the sample, whereas \code{"all"} always removes missing cases.
#' @param na.impute method to impute missing values for the PCA estimates used to initialize the EM algorithm. Note that data are standardized (scaled and centered) beforehand. Available options are:
#' \tabular{llll}{
#' \code{"median"} \tab\tab simple series-wise median imputation. \cr\cr
#' \code{"rnrom"} \tab\tab imputation with random numbers drawn from a standard normal distribution. \cr\cr
#' \code{"med_MA"} \tab\tab values are initially imputed with the median, but then a moving average is applied to smooth the estimates. \cr\cr
#' \code{"med_MA_spline"} \tab\tab "internal" missing values (not at the beginning or end of the sample) are imputed using a cubic spline, whereas missing values at the beginning and end are imputed with the median of the series and smoothed with a moving average.\cr\cr
#' }
#' @param na.impute.MA the order of the (2-sided) moving average applied in \code{na.impute} methods \code{"med_MA"} and \code{"med_MA_spline"}.
#'
#' @details
#' This function efficiently estimates a Dynamic Factor Model with the following classical assumptions:
#' \enumerate{
#' \item Linearity
#' \item Idiosynchratic measurement (observation) errors (no cross-sectional correlation)
#' \item No direct relationship between series and lagged factors (\emph{ceteris paribus})
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
#' @importFrom collapse fscale qsu fvar fmedian qM unattrib na_omit
#' @export

DFM <- function(X, r, p = 1L, ...,
                rQ = c("none", "diagonal", "identity"),
                rR = c("diagonal", "identity", "none"),
                EM.method = c("DGR", "BM", "none"),
                min.iter = 25L, max.iter = 100L, tol = 1e-4,
                max.missing = 0.5,
                na.rm.method = c("LE", "all"),
                na.impute = c("median", "rnrom", "med_MA", "med_MA_spline"),
                na.impute.MA = 3L) {

  rRi <- switch(rR[1L], identity = 0L, diagonal = 1L, none = 2L, stop("Unknown rR option:", rR[1L]))
  rQi <- switch(rQ[1L], identity = 0L, diagonal = 1L, none = 2L, stop("Unknown rQ option:", rQ[1L]))
  BMl <- switch(EM.method[1L], DGR = FALSE, BM = TRUE, none = NA, stop("Unknown EM option:", EM.method[1L]))

  rp <- r * p
  sr <- 1:r
  # srp <- 1:rp
  ax <- attributes(X)
  ilX <- is.list(X)
  Xstat <- qsu(X)
  X <- fscale(qM(X))
  T <- dim(X)[1L]
  n <- dim(X)[2L]

  # Missing values
  X_imp <- X
  na.rm <- NULL
  anymiss <- anyNA(X)
  if(anymiss) { # Missing value removal / imputation
    W <- NULL
    list2env(tsremimpNA(X, max.missing, na.rm.method, na.impute, na.impute.MA),
             envir = environment())
    if(length(na.rm)) X <- X[-na.rm, ]
  }

  # Run PCA to get initial factor estimates:
  v <- svd(X_imp, nu = 0L, nv = min(as.integer(r), n, T))$v
  F_pc <- X_imp %*% v

  # Observation equation -------------------------------
  # Static predictions (all.equal(unattrib(HDB(X_imp, F_pc)), unattrib(F_pc %*% t(v))))
  C <- cbind(v, matrix(0, n, rp-r))
  if(rRi) {
    res <- X_imp - F_pc %*% t(v) # residuals from static predictions
    if(anymiss) res[W] <- NA # Good??? -> Yes, BM do the same...
    R <- if(rRi == 2L) cov(res, use = "pairwise.complete.obs") else diag(fvar(res))
  } else R <- diag(n)

  # Transition equation -------------------------------
  var <- fVAR(F_pc, p)
  A <- rbind(t(var$A), diag(1, rp-r, rp)) # var$A is rp x r matrix
  Q <- matrix(0, rp, rp)
  Q[sr, sr] <- switch(rQi + 1L, diag(r),  diag(fvar(var$res)), cov(var$res))

  # Initial state and state covariance (P) ------------
  F0 <- var$X[1L, ] # rep(0, rp)
  # Kalman gain is normally A %*% t(A) + Q, but here A is somewhat tricky...
  P0 <- matrix(apinv(kronecker(A, A)) %*% unattrib(Q), rp, rp)
  # BM2014: P0 <- matrix(solve(diag(rp^2) - kronecker(A, A)) %*% unattrib(Q), rp, rp)

  ## Run standartized data through Kalman filter and smoother once
  ks_res <- KalmanFilterSmoother(X, C, Q, R, A, F0, P0)

  ## Two-step solution is state mean from the Kalman smoother
  F_kal <- ks_res$Fs

  # Results object for the two-step case
  object_init <- list(X_imp = structure(X_imp,
                                        stats = Xstat,
                                        missing = if(anymiss) W else NULL,
                                        attributes = ax,
                                        is.list = ilX),
                       pca = F_pc,
                       twostep = F_kal[, sr, drop = FALSE],
                       anyNA = anymiss,
                       na.rm = na.rm,
                       EM.method = EM.method[1L],
                       call = match.call())

  # We only report two-step solution
  if(is.na(BMl)) {
  # TODO: Better solution for system matrix estimation after Kalman Filtering and Smoothing?
    final_object <- c(object_init[1:3],
                      list(A = A[sr, , drop = FALSE],
                           C = C[, sr, drop = FALSE],
                           Q = Q[sr, sr, drop = FALSE],
                           R = R),
                      object_init[-(1:3)])
    class(final_object) <- "dfm"
    return(final_object)
  }

  previous_loglik <- -.Machine$double.xmax
  loglik_all <- NULL
  num_iter <- 0L
  converged <- FALSE

  # TODO: What is the good solution with missing values here???
  cpX <- crossprod(X_imp) # <- crossprod(if(anymiss) na_omit(X) else X)
  em_res <- list()
  expr <- if(BMl) .EM_BM else .EM_DGR
  encl <- environment()
  while(num_iter < max.iter && !converged) {

    em_res <- eval(expr, em_res, encl)
    loglik <- em_res$loglik

    ## Iterate at least min.iter times
    converged <- if(num_iter < min.iter) FALSE else
        em_converged(loglik, previous_loglik, tol)
    previous_loglik <- loglik
    loglik_all <- c(loglik_all, loglik)
    num_iter <- num_iter + 1L
  }

  if(converged) message("Converged after ", num_iter, " iterations.")
  else warning("Maximum number of iterations reached.")

  ## Run the Kalman filtering and smoothing step for the last time
  ## with optimal estimates
  F_hat <- eval(.KFS, em_res, encl)$Fs
  final_object <- c(object_init[1:3],
               list(qml = F_hat[, sr, drop = FALSE],
                    A = em_res$A[sr, , drop = FALSE],
                    C = em_res$C[, sr, drop = FALSE],
                    Q = em_res$Q[sr, sr, drop = FALSE],
                    R = em_res$R,
                    loglik = loglik_all,
                    tol = tol,
                    converged = converged),
                    object_init[-(1:3)])

  class(final_object) <- "dfm"
  return(final_object)
}

