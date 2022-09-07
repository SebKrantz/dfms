# Quoting some functions that need to be evaluated iteratively
.EM_DGR <- quote(EMstepDGR(X, A, C, Q, R, F_0, P_0, cpX, n, r, sr, T, rQi, rRi))
.EM_BM <- quote(EMstepBMOPT(X, A, C, Q, R, F_0, P_0, XW0, W, n, r, sr, T, dgind, dnkron, dnkron_ind, rQi, rRi))
.KFS <- quote(SKFS(X, A, C, Q, R, F_0, P_0))


#' Estimate a Dynamic Factor Model
#'
#' Efficient estimation of a Dynamic Factor Model via the EM Algorithm - on stationary data
#' with time-invariant system matrices and classical assumptions, while permitting missing data.
#'
#' @param X a \code{T x n} data matrix or frame.
#' @param r number of factors.
#' @param p number of lags in factor VAR.
#' @param \dots further arguments to be added here in the future, such as further estimation methods or block-structures.
#' @param rQ restrictions on the state (transition) covariance matrix (Q).
#' @param rR restrictions on the observation (measurement) covariance matrix (R).
#' @param em.method character. The implementation of the Expectation Maximization Algorithm used. The options are:
#'    \tabular{llll}{
#' \code{"DGR"} \tab\tab The classical EM implementation of Doz, Giannone and Reichlin (2012). This implementation is efficient and quite robust, missing values are removed on a casewise basis in the Kalman Filter and Smoother, but not explicitly accounted for in EM iterations. \cr\cr
#' \code{"BM"} \tab\tab The modified EM algorithm of Banbura and Modugno (2014) which also accounts for missing data in the EM iterations. Optimal for datasets with arbitrary patterns of missing data e.g. datasets with series at different frequencies.  \cr\cr
#' \code{"none"} \tab\tab Performs no EM iterations and just returns the Two-Step estimates from running the data through the Kalman Filter and Smoother once as in
#' Doz, Giannone and Reichlin (2011) (the Kalman Filter is Initialized with system matrices obtained from a regression and VAR on PCA factor estimates).
#' This yields significant performance gains over the iterative methods. Final system matrices are estimated by running a regression and a VAR on the smoothed factors.  \cr\cr
#' }
#' @param min.iter integer. Minimum number of EM iterations (to ensure a convergence path).
#' @param max.iter integer. Maximum number of EM iterations.
#' @param tol numeric. EM convergence tolerance.
#' @param check.increased logical. Check if likelihood has increased. Passed to \code{\link{em_converged}}.
#' @param max.missing numeric. Proportion of series missing for a case to be considered missing.
#' @param na.rm.method character. Method to apply concerning missing cases selected through \code{max.missing}: \code{"LE"} only removes cases at the beginning or end of the sample, whereas \code{"all"} always removes missing cases.
#' @param na.impute character. Method to impute missing values for the PCA estimates used to initialize the EM algorithm. Note that data are standardized (scaled and centered) beforehand. Available options are:
#'    \tabular{llll}{
#' \code{"median"} \tab\tab simple series-wise median imputation. \cr\cr
#' \code{"rnorm"} \tab\tab imputation with random numbers drawn from a standard normal distribution. \cr\cr
#' \code{"median.ma"} \tab\tab values are initially imputed with the median, but then a moving average is applied to smooth the estimates. \cr\cr
#' \code{"median.ma.spline"} \tab\tab "internal" missing values (not at the beginning or end of the sample) are imputed using a cubic spline, whereas missing values at the beginning and end are imputed with the median of the series and smoothed with a moving average.\cr\cr
#' }
#' @param ma.terms the order of the (2-sided) moving average applied in \code{na.impute} methods \code{"median.ma"} and \code{"median.ma.spline"}.
#'
#' @details
#' This function efficiently estimates a Dynamic Factor Model with the following classical assumptions:
#' \enumerate{
#' \item Linearity
#' \item Idiosynchratic measurement (observation) errors (\emph{R} is diagonal)
#' \item No direct relationship between series and lagged factors (\emph{ceteris paribus} contemporaneous factors)
#' \item No relationship between lagged error terms in the either measurement or transition equation (no serial correlation)
#' }
#' Factors are allowed to evolve in a \eqn{VAR(p)} process, and data is standardized (scaled and centered) before estimation (removing the need of intercept terms).
#' By assumptions 1-4, this translates into the following dynamic form:
#'
#' \deqn{\textbf{x}_t = \textbf{C}_0 \textbf{f}_t + \textbf{e}_t \sim N(\textbf{0}, \textbf{R})}{x(t) = C0 f(t) + e(t) ~ N(0, R)}
#' \deqn{\textbf{f}_t = \sum_{i=1}^p \textbf{A}_p \textbf{f}_{t-p} + \textbf{u}_t \sim N(0, \textbf{Q}_0)}{f(t) = A1 f(t-1) + \dots + Ap f(t-p) + u(t) ~ N(0, Q0)}
#'
#' where the first equation is called the measurement or observation equation and the second equation is called transition, state or process equation, and
#'
#' \tabular{llll}{
#'  \eqn{n} \tab\tab number of series in \eqn{\textbf{x}_t}{x(t)} (\eqn{r} and \eqn{p} as the arguments to \code{DFM}).\cr\cr
#'  \eqn{\textbf{x}_t}{x(t)} \tab\tab \eqn{n \times 1}{n x 1} vector of observed series at time \eqn{t}{t}: \eqn{(x_{1t}, \dots, x_{nt})'}{(x1(t), \dots, xn(t))'}. Some observations can be missing.  \cr\cr
#'  \eqn{\textbf{f}_t}{f(t)} \tab\tab \eqn{r \times 1}{r x 1} vector of factors at time \eqn{t}{t}: \eqn{(f_{1t}, \dots, f_{rt})'}{(f1(t), \dots, fr(t))'}.\cr\cr
#'  \eqn{\textbf{C}_0}{C0} \tab\tab \eqn{n \times r}{n x r} measurement (observation) matrix.\cr\cr
#'  \eqn{\textbf{A}_j}{Aj} \tab\tab \eqn{r \times r}{r x r} state transition matrix at lag \eqn{j}{j}. \cr\cr
#'  \eqn{\textbf{Q}_0}{Q0} \tab\tab \eqn{r \times r}{r x r} state covariance matrix.\cr\cr
#'  \eqn{\textbf{R}}{R} \tab\tab \eqn{n \times n}{n x n} measurement (observation) covariance matrix. It is diagonal by assumption 2 that \eqn{E[\textbf{x}_{it}|\textbf{x}_{-i,t},\textbf{x}_{i,t-1}, \dots, \textbf{f}_t, \textbf{f}_{t-1}, \dots] = \textbf{Cf}_t \forall i}{E[xi(t)|x(-i)(t), x(t-1), \dots, f(t), f(t-1), \dots] = C f(t)}.\cr\cr  % This assumption is also referred to as 'Exact DFM' by Stock & Watson (2016), where all correlation between the series is accounted for by the latent factors.
#' }
#'
#'
#' This model can be estimated using a classical form of the Kalman Filter and the Expectation Maximization (EM) algorithm, after transforming it to State-Space (stacked, VAR(1)) form:
#'
#' \deqn{\textbf{x}_t = \textbf{C} \textbf{F}_t + \textbf{e}_t \sim N(\textbf{0}, \textbf{R})}{x(t) = C F(t) + e(t) ~ N(0, R)}
#' \deqn{\textbf{F}_t = \textbf{A F}_{t-1} + \textbf{u}_t \sim N(0, \textbf{Q})}{F(t) = A F(t-1) + u(t) ~ N(0, Q)}
#'
#' where
#' \tabular{llll}{
#'  \eqn{n} \tab\tab number of series in \eqn{\textbf{x}_t}{x(t)} (\eqn{r} and \eqn{p} as the arguments to \code{DFM}).\cr\cr
#'  \eqn{\textbf{x}_t}{x(t)} \tab\tab \eqn{n \times 1}{n x 1} vector of observed series at time \eqn{t}{t}: \eqn{(x_{1t}, \dots, x_{nt})'}{(x1(t), \dots, xn(t))'}. Some observations can be missing.  \cr\cr
#'  \eqn{\textbf{F}_t}{F(t)} \tab\tab \eqn{rp \times 1}{rp x 1} vector of stacked factors at time \eqn{t}{t}: \eqn{(f_{1t}, \dots, f_{rt}, f_{1,t-1}, \dots, f_{r,t-1}, \dots, f_{1,t-p}, \dots, f_{r,t-p})'}{(f1(t), \dots, fr(t), f1(t-1), \dots, fr(t-1), \dots, f1(t-p), \dots, fr(t-p))'}.\cr\cr
#'  \eqn{\textbf{C}}{C} \tab\tab \eqn{n \times rp}{n x rp} observation matrix. Only the first \eqn{n \times r}{n x r} terms are non-zero, by assumption 3 that \eqn{E[\textbf{x}_t|\textbf{F}_t] = E[\textbf{x}_t|\textbf{f}_t]}{E[X(t)|F(t)] = E[X(t)|f(t)]} (no relationship of observed series with lagged factors given contemporaneous factors).\cr\cr
#'  \eqn{\textbf{A}}{A} \tab\tab stacked \eqn{rp \times rp}{rp x rp} state transition matrix consisting of 3 parts: the top \eqn{r \times rp}{r x rp} part provides the dynamic relationships captured by \eqn{(\textbf{A}_1, \dots, \textbf{A}_p)}{(A1, \dots, Ap)} in the dynamic form, the terms \code{A[(r+1):rp, 1:(rp-r)]} constitute an \eqn{(rp-r) \times (rp-r)}{(rp-r) x (rp-r)} identity matrix mapping all lagged factors to their known values at times t. The remaining part \code{A[(rp-r+1):rp, (rp-r+1):rp]} is an \eqn{r \times r}{r x r} matrix of zeros. \cr\cr
#'  \eqn{\textbf{Q}}{Q} \tab\tab \eqn{rp \times rp}{rp x rp} state covariance matrix. The top \eqn{r \times r}{r x r} part gives the contemporaneous relationships, the rest are zeros by assumption 4.\cr\cr % that \eqn{E[\textbf{f}_t|\textbf{F}_{t-1}] = E[\textbf{f}_t|\textbf{f}_{t-1}] = \textbf{A}_1 \textbf{f}_{t-1}}{E[f(t)|F(t-1)] = E[f(t)|f(t-1)] = A1 f(t-1)} (all relationships between lagged factors are captured in \eqn{\textbf{A}_1}{A1}).\cr\cr
#'  \eqn{\textbf{R}}{R} \tab\tab \eqn{n \times n}{n x n} observation covariance matrix. It is diagonal by assumption 2 and identical to \eqn{\textbf{R}}{R} as stated in the dynamic form.\cr\cr
#' }
#'
#' @returns A list-like object of class 'dfm' with the following elements:
#' \tabular{llll}{
#'  \code{X_imp} \tab\tab \eqn{T \times n}{T x n} matrix with the imputed and standardized (scaled and centered) data - with attributes attached allowing reconstruction of the original data:
#'    \itemize{
#'       \item \code{"stats"} is a \eqn{n \times 5}{n x 5} matrix of summary statistics of class \code{"qsu"} (see \code{\link[collapse]{qsu}}).\cr
#'       \item \code{"missing"} is a \eqn{T \times n}{T x n} logical matrix indicating missing or infinite values in the original data (which are imputed in \code{X_imp}).\cr
#'       \item \code{"attributes"} contains the \code{\link{attributes}} or the original data input.\cr
#'       \item \code{"is.list"} is a logical value indicating whether the original data input was a list / data frame. \cr
#'    } \cr\cr
#'  \code{eigen} \tab\tab \code{eigen(cov(X_imp))}. \cr\cr
#'  \code{F_pca} \tab\tab \eqn{T \times r}{T x r} matrix of principal component factor estimates - \code{X_imp \%*\% eigen$vectors}. \cr\cr
#'  \code{P_0} \tab\tab \eqn{r \times r}{r x r} initial factor covariance matrix estimate based on PCA results. \cr\cr
#'  \code{F_twostep} \tab\tab \eqn{T \times r}{T x r} matrix two-step factor estimates as in Doz, Giannone and Reichlin (2011) - obtained from running the data through the Kalman Filter and Smoother once, where the Filter is initialized with results from PCA. \cr\cr
#'  \code{P_twostep} \tab\tab \eqn{r \times r \times T}{r x r x T} covariance matrices of two-step factor estimates. \cr\cr
#'  \code{F_qml} \tab\tab \eqn{T \times r}{T x r} matrix of quasi-maximum likelihood factor estimates - obtained by iteratively Kalman Filtering and Smoothing the factor estimates until EM convergence. \cr\cr
#'  \code{P_qml} \tab\tab \eqn{r \times r \times T}{r x r x T} covariance matrices of QML factor estimates. \cr\cr
#'  \code{A} \tab\tab \eqn{r \times rp}{r x rp} factor transition matrix.\cr\cr
#'  \code{C} \tab\tab \eqn{n \times r}{n x r} observation matrix.\cr\cr
#'  \code{Q} \tab\tab \eqn{r \times r}{r x r} state (error) covariance matrix.\cr\cr
#'  \code{R} \tab\tab \eqn{n \times n}{n x n} observation (error) covariance matrix.\cr\cr
#'  \code{loglik} \tab\tab vector of log-likelihoods - one for each EM iteration. The final value corresponds to the log-likelihood of the reported model.\cr\cr
#'  \code{tol} \tab\tab The numeric convergence tolerance used.\cr\cr
#'  \code{converged} \tab\tab single logical valued indicating whether the EM algorithm converged (within \code{max.iter} iterations subject to \code{tol}).\cr\cr
#'  \code{anyNA} \tab\tab single logical valued indicating whether there were any missing values in the data. If \code{FALSE}, \code{X_imp} is simply the original data in matrix form, and does not have the \code{"missing"} attribute attached.\cr\cr
#'  \code{na.rm} \tab\tab vector of any cases that were removed beforehand (subject to \code{max.missing} and \code{na.rm.method}). If no cases were removed the slot is \code{NULL}. \cr\cr
#'  \code{em.method} \tab\tab The EM method used.\cr\cr
#'  \code{call} \tab\tab call object obtained from \code{match.call()}.\cr\cr
#' }
#'
#' @references
#' Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large approximate dynamic factor models based on Kalman filtering. \emph{Journal of Econometrics, 164}(1), 188-205.
#'
#' Doz, C., Giannone, D., & Reichlin, L. (2012). A quasi-maximum likelihood approach for large, approximate dynamic factor models. \emph{Review of economics and statistics, 94}(4), 1014-1024.
#'
#' Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor models on datasets with arbitrary pattern of missing data. \emph{Journal of Applied Econometrics, 29}(1), 133-160.
#'
#' Stock, J. H., & Watson, M. W. (2016). Dynamic Factor Models, Factor-Augmented Vector Autoregressions, and Structural Vector Autoregressions in Macroeconomics. \emph{Handbook of Macroeconomics, 2}, 415–525. https://doi.org/10.1016/bs.hesmac.2016.04.002
#'
#' @useDynLib DFM, .registration = TRUE
#' @importFrom collapse fscale qsu fvar fmedian fmedian.default qM unattrib na_omit %=% %+=% %/=% whichv setColnames setDimnames
#' @importFrom grDevices rainbow
#' @importFrom graphics abline legend lines par
#' @importFrom stats filter residuals rnorm spline ts.plot
#'
#' @examples
#' library(magrittr)
#' library(xts)
#' library(vars)
#'
#' # BM14 Replication Data. Constructing the database:
#' BM14 = merge(BM14_M, BM14_Q)
#' BM14[, BM14_Models$log_trans] %<>% log()
#' BM14[, BM14_Models$freq == "M"] %<>% diff()
#' BM14[, BM14_Models$freq == "Q"] %<>% diff(3)
#'
#'
#' ### Small Model ---------------------------------------
#'
#' # IC for number of factors
#' IC_small = ICr(BM14[, BM14_Models$small], max.r = 5)
#' plot(IC_small)
#' screeplot(IC_small)
#'
#' # I take 2 factors. Now number of lags
#' VARselect(IC_small$F_pca[, 1:2])
#'
#' # Estimating the model with 2 factors and 3 lags
#' dfm_small = DFM(BM14[, BM14_Models$small], 2, 3, em.method = "BM")
#'
#' # Inspecting the model
#' summary(dfm_small)
#' plot(dfm_small)  # Factors and data
#' plot(dfm_small, method = "all", type = "individual") # Factor estimates
#' plot(dfm_small, type = "residual") # Residuals from factor predictions
#'
#' # 10 periods ahead forecast
#' plot(predict(dfm_small), xlim = c(300, 370))
#'
#'
#' ### Medium-Sized Model ---------------------------------
#'
#' # IC for number of factors
#' IC_medium = ICr(BM14[, BM14_Models$medium])
#' plot(IC_medium)
#' screeplot(IC_medium)
#'
#' # I take 3 factors. Now number of lags
#' VARselect(IC_medium$F_pca[, 1:3])
#'
#' # Estimating the model with 3 factors and 3 lags
#' dfm_medium = DFM(BM14[, BM14_Models$medium], 3, 3, em.method = "BM")
#'
#' # Inspecting the model
#' summary(dfm_medium)
#' plot(dfm_medium)  # Factors and data
#' plot(dfm_medium, method = "all", type = "individual") # Factor estimates
#' plot(dfm_medium, type = "residual") # Residuals from factor predictions
#'
#' # 10 periods ahead forecast
#' plot(predict(dfm_medium), xlim = c(300, 370))
#'
#'
#' ### Large Model ---------------------------------
#' \donttest{
#' # IC for number of factors
#' IC_large = ICr(BM14)
#' plot(IC_large)
#' screeplot(IC_large)
#'
#' # I take 6 factors. Now number of lags
#' VARselect(IC_large$F_pca[, 1:6])
#'
#' # Estimating the model with 6 factors and 3 lags
#' dfm_large = DFM(BM14, 6, 3, em.method = "BM")
#'
#' # Inspecting the model
#' summary(dfm_large)
#' plot(dfm_large)  # Factors and data
#' # plot(dfm_large, method = "all", type = "individual") # Factor estimates
#' plot(dfm_large, type = "residual") # Residuals from factor predictions
#'
#' # 10 periods ahead forecast
#' plot(predict(dfm_large), xlim = c(300, 370))
#' }
#' @export

DFM <- function(X, r, p = 1L, ...,
                rQ = c("none", "diagonal", "identity"),
                rR = c("diagonal", "identity", "none"),
                em.method = c("DGR", "BM", "none"),
                min.iter = 25L,
                max.iter = 100L,
                tol = 1e-4,
                check.increased = FALSE,
                max.missing = 0.8,
                na.rm.method = c("LE", "all"),
                na.impute = c("median.ma.spline", "median.ma", "median", "rnorm"),
                ma.terms = 3L) {

  rRi <- switch(rR[1L], identity = 0L, diagonal = 1L, none = 2L, stop("Unknown rR option:", rR[1L]))
  rQi <- switch(rQ[1L], identity = 0L, diagonal = 1L, none = 2L, stop("Unknown rQ option:", rQ[1L]))
  BMl <- switch(em.method[1L], DGR = FALSE, BM = TRUE, none = NA, stop("Unknown EM option:", em.method[1L]))

  rp <- r * p
  sr <- 1:r
  fnam <- paste0("f", sr)
  unam <- paste0("u", sr)
  # srp <- 1:rp
  ax <- attributes(X)
  ilX <- is.list(X)
  Xstat <- qsu(X)
  X <- fscale(qM(X))
  Xnam <- dimnames(X)[[2L]]
  dimnames(X) <- NULL
  n <- ncol(X)

  # Missing values
  X_imp <- X
  na.rm <- NULL
  anymiss <- anyNA(X)
  if(anymiss) { # Missing value removal / imputation
    W <- NULL
    list2env(tsremimpNA(X, max.missing, na.rm.method, na.impute, ma.terms),
             envir = environment())
    if(length(na.rm)) X <- X[-na.rm, , drop = FALSE]
  }
  T <- nrow(X)

  # Run PCA to get initial factor estimates:
  # v <- svd(X_imp, nu = 0L, nv = min(as.integer(r), n, T))$v # Not faster than eigen...
  eigen_decomp = eigen(cov(X_imp), symmetric = TRUE)
  v = eigen_decomp$vectors[, sr, drop = FALSE]
  # d = eigen_decomp$values[sr]
  F_pc <- X_imp %*% v

  # Observation equation -------------------------------
  # Static predictions (all.equal(unattrib(HDB(X_imp, F_pc)), unattrib(F_pc %*% t(v))))
  C <- cbind(v, matrix(0, n, rp-r))
  if(rRi) {
    res <- X_imp - tcrossprod(F_pc, v) # residuals from static predictions
    if(anymiss) res[W] <- NA
    R <- if(rRi == 2L) cov(res, use = "pairwise.complete.obs") else diag(fvar(res))
  } else R <- diag(n)

  # Transition equation -------------------------------
  var <- .VAR(F_pc, p)
  A <- rbind(t(var$A), diag(1, rp-r, rp)) # var$A is rp x r matrix
  Q <- matrix(0, rp, rp)
  Q[sr, sr] <- switch(rQi + 1L, diag(r),  diag(fvar(var$res)), cov(var$res))

  # Initial state and state covariance (P) ------------
  F_0 <- rep(0, rp) # ar$X[1L, ] #
  # TODO: Kalman gain is normally A %*% t(A) + Q, but here A is somewhat tricky...
  # -> cannot be that this influences the twostep estimates. Need to chosse one way to initialize. Better BM14..
  P_0 <- matrix(ainv(diag(rp^2) - kronecker(A,A)) %*% unattrib(Q), rp, rp)
  # if(!is.na(BMl) && BMl)  else matrix(apinv(kronecker(A,A)) %*% unattrib(Q), rp, rp)

  ## Run standartized data through Kalman filter and smoother once
  kfs_res <- SKFS(X, A, C, Q, R, F_0, P_0, FALSE)

  ## Two-step solution is state mean from the Kalman smoother
  F_kal <- kfs_res$F_smooth[, sr, drop = FALSE]

  # Results object for the two-step case
  object_init <- list(X_imp = structure(X_imp,
                                        dimnames = list(NULL, Xnam),
                                        stats = Xstat,
                                        missing = if(anymiss) W else NULL,
                                        attributes = ax,
                                        is.list = ilX),
                       eigen = eigen_decomp,
                       F_pca = setColnames(F_pc, paste0("PC", sr)),
                       P_0 = setDimnames(P_0[sr, sr, drop = FALSE], list(fnam, fnam)),
                       F_twostep = setColnames(F_kal, fnam),
                       P_twostep = setDimnames(kfs_res$P_smooth[sr, sr,, drop = FALSE], list(fnam, fnam, NULL)),
                       anyNA = anymiss,
                       na.rm = na.rm,
                       em.method = em.method[1L],
                       call = match.call())

  # We only report two-step solution
  if(is.na(BMl)) {
  # TODO: Better solution for system matrix estimation after Kalman Filtering and Smoothing? (could take matrices from Kalman Filter, but that would be before smoothing)
    var <- .VAR(F_kal, p)
    beta <- ainv(crossprod(F_kal)) %*% crossprod(F_kal, if(anymiss) replace(X_imp, W, 0) else X_imp) # good??
    Q <- switch(rQi + 1L, diag(r),  diag(fvar(var$res)), cov(var$res))
    if(rRi) {
      res <- X_imp - F_kal %*% beta
      if(anymiss) res[W] <- NA
      R <- if(rRi == 2L) cov(res, use = "pairwise.complete.obs") else diag(fvar(res))
    } else R <- diag(n)
    final_object <- c(object_init[1:6],
                      list(A = setDimnames(t(var$A), lagnam(fnam, p)), # A[sr, , drop = FALSE],
                           C = setDimnames(t(beta), list(Xnam, fnam)), # C[, sr, drop = FALSE],
                           Q = setDimnames(Q, list(unam, unam)),       # Q[sr, sr, drop = FALSE],
                           R = setDimnames(R, list(Xnam, Xnam))),
                      object_init[-(1:6)])
    class(final_object) <- "dfm"
    return(final_object)
  }

  previous_loglik <- -Inf # .Machine$double.xmax
  loglik_all <- NULL
  num_iter <- 0L
  converged <- FALSE

  if(BMl) {
    expr <- .EM_BM
    dnkron <- matrix(1, r, r) %x% diag(n) # Used to be inside EMstep
    dnkron_ind <- whichv(dnkron, 1)
    XW0 <- X_imp
    dgind <- 0:(n-1) * n + 1:n
    if(anymiss) XW0[W] <- 0 else W <- is.na(X) # TODO: think about this...
  } else {
    expr <- .EM_DGR
    # TODO: What is the good solution with missing values here?? -> Zeros are ignored in crossprod, so it's like skipping those obs
    cpX <- crossprod(if(anymiss) replace(X_imp, W, 0) else X_imp) # <- crossprod(if(anymiss) na_omit(X) else X)
  }
  em_res <- list()
  encl <- environment()
  while(num_iter < max.iter && !converged) {
    # print(num_iter)
    em_res <- eval(expr, em_res, encl)
    # if(is.null(em_res[[1]])) return(em_res)
    loglik <- em_res$loglik

    ## Iterate at least min.iter times
    converged <- if(num_iter < min.iter) FALSE else
        em_converged(loglik, previous_loglik, tol, check.increased)[1L]
    previous_loglik <- loglik
    loglik_all <- c(loglik_all, loglik)
    num_iter <- num_iter + 1L
  }

  if(converged) message("Converged after ", num_iter, " iterations.")
  else warning("Maximum number of iterations reached.")

  ## Run the Kalman filtering and smoothing step for the last time
  ## with optimal estimates
  kfs_res <- eval(.KFS, em_res, encl)

  final_object <- c(object_init[1:6],
               list(F_qml = setColnames(kfs_res$F_smooth[, sr, drop = FALSE], fnam),
                    P_qml = setDimnames(kfs_res$P_smooth[sr, sr,, drop = FALSE], list(fnam, fnam, NULL)),
                    A = setDimnames(em_res$A[sr, , drop = FALSE], lagnam(fnam, p)),
                    C = setDimnames(em_res$C[, sr, drop = FALSE], list(Xnam, fnam)),
                    Q = setDimnames(em_res$Q[sr, sr, drop = FALSE], list(unam, unam)),
                    R = setDimnames(em_res$R, list(Xnam, Xnam)),
                    loglik = loglik_all,
                    tol = tol,
                    converged = converged),
                    object_init[-(1:6)])

  class(final_object) <- "dfm"
  return(final_object)
}

