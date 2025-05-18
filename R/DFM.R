# Quoting some functions that need to be evaluated iteratively
.EM_DGR <- quote(EMstepDGR(X, A, C, Q, R, F_0, P_0, cpX, n, r, sr, TT, rQi, rRi))
.EM_BM <- quote(EMstepBMOPT(X, A, C, Q, R, F_0, P_0, XW0, W, n, r, sr, TT, dgind, dnkron, dnkron_ind, rQi, rRi))
.EM_BM_idio <- quote(EMstepBMidio(X, A, C, Q, R, F_0, P_0, XW0, W, dgind, dnkron, dnkron_ind, r, p, n, sr, TT, rQi, rRi))
.EM_BM_MQ <- quote(EMstepBMMQ(X, A, C, Q, R, F_0, P_0, XW0, NW, dgind, dnkron, dnkron_ind, r, p, R_mat, n, nq, sr, TT, rQi, rRi))
.EM_BM_MQ_idio <- quote(EMstepBMMQidio(X, A, C, Q, R, F_0, P_0, XW0, NW, dgind, dnkron, dnkron_ind, r, p, R_mat, n, nq, sr, TT, rQi, rRi))
.KFS <- quote(SKFS(X, A, C, Q, R, F_0, P_0))

#' Estimate a Dynamic Factor Model
#'
#' Efficient estimation of a Dynamic Factor Model via the EM Algorithm - on stationary data
#' with time-invariant system matrices and classical assumptions, while permitting missing data.
#'
#' @param X a \code{T x n} numeric data matrix or frame of stationary time series. May contain missing values. \emph{Note} that data is internally standardized (scaled and centered) before estimation.
#' @param r integer. Number of factors.
#' @param p integer. Number of lags in factor VAR.
#' @param \dots (optional) arguments to \code{\link{tsnarmimp}}. The default settings impute internal missing values with a cubic spline and the edges with the median and a 3-period moving average.
#' @param idio.ar1 logical. Model observation errors as AR(1) processes: \eqn{e_t = \rho e_{t-1} + v_t}{e(t) = rho e(t-1) + v(t)}. \emph{Note} that this substantially increases computation time, and is generally not needed if \code{n} is large (>30). See theoretical vignette for details.
#' @param quarterly.vars character. Names of quarterly variables in \code{X} (if any). Monthly variables should be to the left of the quarterly variables in the data matrix and quarterly observations should be provided every 3rd period.
#' @param rQ character. Restrictions on the state (transition) covariance matrix (Q).
#' @param rR character. Restrictions on the observation (measurement) covariance matrix (R).
#' @param em.method character. The implementation of the Expectation Maximization Algorithm used. The options are:
#'    \tabular{llll}{
#' \code{"auto"} \tab\tab Automatic selection: \code{"BM"} if \code{anyNA(X)}, else \code{"DGR"}. \cr\cr
#' \code{"DGR"} \tab\tab The classical EM implementation of Doz, Giannone and Reichlin (2012). This implementation is efficient and quite robust, missing values are removed on a casewise basis in the Kalman Filter and Smoother, but not explicitly accounted for in EM iterations. \cr\cr
#' \code{"BM"} \tab\tab The modified EM algorithm of Banbura and Modugno (2014) which also accounts for missing data in the EM iterations. Optimal for datasets with systematically missing data e.g. datasets with ragged edges or series at different frequencies.  \cr\cr
#' \code{"none"} \tab\tab Performs no EM iterations and just returns the Two-Step estimates from running the data through the Kalman Filter and Smoother once as in
#' Doz, Giannone and Reichlin (2011) (the Kalman Filter is Initialized with system matrices obtained from a regression and VAR on PCA factor estimates).
#' This yields significant performance gains over the iterative methods. Final system matrices are estimated by running a regression and a VAR on the smoothed factors.  \cr\cr
#' }
#' @param min.iter integer. Minimum number of EM iterations (to ensure a convergence path).
#' @param max.iter integer. Maximum number of EM iterations.
#' @param tol numeric. EM convergence tolerance.
#' @param pos.corr logical. Increase the likelihood that factors correlate positively with the data, by scaling the eigenvectors such that the principal components (used to initialize the Kalman Filter) co-vary positively with the row-means of the standardized data.
#' @param check.increased logical. Check if likelihood has increased. Passed to \code{\link{em_converged}}. If \code{TRUE}, the algorithm only terminates if convergence was reached with decreasing likelihood.
#'
#' @details
#' This function efficiently estimates a Dynamic Factor Model with the following classical assumptions:
#' \enumerate{
#' \item Linearity
#' \item Idiosynchratic measurement (observation) errors (\emph{R} is diagonal)
#' \item No direct relationship between series and lagged factors (\emph{ceteris paribus} contemporaneous factors)
#' \item No relationship between lagged error terms in the either measurement or transition equation (no serial correlation), unless explicitly modeled as AR(1) processes using \code{idio.ar1 = TRUE}.
#' }
#' Factors are allowed to evolve in a \eqn{VAR(p)} process, and data is internally standardized (scaled and centered) before estimation (removing the need of intercept terms).
#' By assumptions 1-4, this translates into the following dynamic form:
#'
#' \deqn{\textbf{x}_t = \textbf{C}_0 \textbf{f}_t + \textbf{e}_t \ \sim\  N(\textbf{0}, \textbf{R})}{x(t) = C0 f(t) + e(t) ~ N(0, R)}
#' \deqn{\textbf{f}_t = \sum_{j=1}^p \textbf{A}_j \textbf{f}_{t-j} + \textbf{u}_t \ \sim\  N(\textbf{0}, \textbf{Q}_0)}{f(t) = A1 f(t-1) + \dots + Ap f(t-p) + u(t) ~ N(0, Q0)}
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
#'  \eqn{\textbf{R}}{R} \tab\tab \eqn{n \times n}{n x n} measurement (observation) covariance matrix. It is diagonal by assumption 2 that \eqn{E[\textbf{x}_{it}|\textbf{x}_{-i,t},\textbf{x}_{i,t-1}, \dots, \textbf{f}_t, \textbf{f}_{t-1}, \dots] = \textbf{Cf}_t \forall i}{E[xi(t)|x(-i)(t), x(t-1), \dots, f(t), f(t-1), \dots] = C f(t)}.\cr\cr
#' }
# This assumption is also referred to as 'Exact DFM' by Stock & Watson (2016), where all correlation between the series is accounted for by the latent factors.
#'
#' This model can be estimated using a classical form of the Kalman Filter and the Expectation Maximization (EM) algorithm, after transforming it to State-Space (stacked, VAR(1)) form:
#'
#' \deqn{\textbf{x}_t = \textbf{C} \textbf{F}_t + \textbf{e}_t \ \sim\  N(\textbf{0}, \textbf{R})}{x(t) = C F(t) + e(t) ~ N(0, R)}
#' \deqn{\textbf{F}_t = \textbf{A F}_{t-1} + \textbf{u}_t \ \sim\  N(\textbf{0}, \textbf{Q})}{F(t) = A F(t-1) + u(t) ~ N(0, Q)}
#'
#' where
#' \tabular{llll}{
#'  \eqn{n} \tab\tab number of series in \eqn{\textbf{x}_t}{x(t)} (\eqn{r} and \eqn{p} as the arguments to \code{DFM}).\cr\cr
#'  \eqn{\textbf{x}_t}{x(t)} \tab\tab \eqn{n \times 1}{n x 1} vector of observed series at time \eqn{t}{t}: \eqn{(x_{1t}, \dots, x_{nt})'}{(x1(t), \dots, xn(t))'}. Some observations can be missing.  \cr\cr
#'  \eqn{\textbf{F}_t}{F(t)} \tab\tab \eqn{rp \times 1}{rp x 1} vector of stacked factors at time \eqn{t}{t}: \eqn{(f_{1t}, \dots, f_{rt}, f_{1,t-1}, \dots, f_{r,t-1}, \dots, f_{1,t-p}, \dots, f_{r,t-p})'}{(f1(t), \dots, fr(t), f1(t-1), \dots, fr(t-1), \dots, f1(t-p), \dots, fr(t-p))'}.\cr\cr
#'  \eqn{\textbf{C}}{C} \tab\tab \eqn{n \times rp}{n x rp} observation matrix. Only the first \eqn{n \times r}{n x r} terms are non-zero, by assumption 3 that \eqn{E[\textbf{x}_t|\textbf{F}_t] = E[\textbf{x}_t|\textbf{f}_t]}{E[X(t)|F(t)] = E[X(t)|f(t)]} (no relationship of observed series with lagged factors given contemporaneous factors).\cr\cr
#'  \eqn{\textbf{A}}{A} \tab\tab stacked \eqn{rp \times rp}{rp x rp} state transition matrix consisting of 3 parts: the top \eqn{r \times rp}{r x rp} part provides the dynamic relationships captured by \eqn{(\textbf{A}_1, \dots, \textbf{A}_p)}{(A1, \dots, Ap)} in the dynamic form, the terms \code{A[(r+1):rp, 1:(rp-r)]} constitute an \eqn{(rp-r) \times (rp-r)}{(rp-r) x (rp-r)} identity matrix mapping all lagged factors to their known values at times t. The remaining part \code{A[(rp-r+1):rp, (rp-r+1):rp]} is an \eqn{r \times r}{r x r} matrix of zeros. \cr\cr
#'  \eqn{\textbf{Q}}{Q} \tab\tab \eqn{rp \times rp}{rp x rp} state covariance matrix. The top \eqn{r \times r}{r x r} part gives the contemporaneous relationships, the rest are zeros by assumption 4.\cr\cr
#'  \eqn{\textbf{R}}{R} \tab\tab \eqn{n \times n}{n x n} observation covariance matrix. It is diagonal by assumption 2 and identical to \eqn{\textbf{R}}{R} as stated in the dynamic form.\cr\cr
#' }
#  that \eqn{E[\textbf{f}_t|\textbf{F}_{t-1}] = E[\textbf{f}_t|\textbf{f}_{t-1}] = \textbf{A}_1 \textbf{f}_{t-1}}{E[f(t)|F(t-1)] = E[f(t)|f(t-1)] = A1 f(t-1)} (all relationships between lagged factors are captured in \eqn{\textbf{A}_1}{A1}).\cr\cr
#' The filter is initialized with PCA estimates on the imputed dataset---see \code{\link{SKFS}} for a complete code example.
#'
#'
#' @returns A list-like object of class 'dfm' with the following elements:
#'  \item{\code{X_imp}}{\eqn{T \times n}{T x n} matrix with the imputed and standardized (scaled and centered) data---after applying \code{\link{tsnarmimp}}. It has attributes attached allowing for reconstruction of the original data:
#'  \tabular{llll}{
#'      \code{"stats"} \tab\tab is a \eqn{n \times 5}{n x 5} matrix of summary statistics of class \code{"qsu"} (see \code{\link[collapse]{qsu}}).\cr\cr
#'      \code{"missing"} \tab\tab is a \eqn{T \times n}{T x n} logical matrix indicating missing or infinite values in the original data (which are imputed in \code{X_imp}).\cr\cr
#'      \code{"attributes"} \tab\tab contains the \code{\link{attributes}} of the original data input.\cr\cr
#'      \code{"is.list"} \tab\tab is a logical value indicating whether the original data input was a list / data frame. \cr\cr }
#'      }
#'  \item{\code{eigen}}{\code{eigen(cov(X_imp))}. }
#'  \item{\code{F_pca}}{\eqn{T \times r}{T x r} matrix of principal component factor estimates - \code{X_imp \%*\% eigen$vectors}. }
#'  \item{\code{P_0}}{\eqn{r \times r}{r x r} initial factor covariance matrix estimate based on PCA results. }
#'  \item{\code{F_2s}}{\eqn{T \times r}{T x r} matrix two-step factor estimates as in Doz, Giannone and Reichlin (2011) - obtained from running the data through the Kalman Filter and Smoother once, where the Filter is initialized with results from PCA. }
#'  \item{\code{P_2s}}{\eqn{r \times r \times T}{r x r x T} covariance matrices of two-step factor estimates. }
#'  \item{\code{F_qml}}{\eqn{T \times r}{T x r} matrix of quasi-maximum likelihood factor estimates - obtained by iteratively Kalman Filtering and Smoothing the factor estimates until EM convergence. }
#'  \item{\code{P_qml}}{\eqn{r \times r \times T}{r x r x T} covariance matrices of QML factor estimates. }
#'  \item{\code{A}}{\eqn{r \times rp}{r x rp} factor transition matrix.}
#'  \item{\code{C}}{\eqn{n \times r}{n x r} observation matrix.}
#'  \item{\code{Q}}{\eqn{r \times r}{r x r} state (error) covariance matrix.}
#'  \item{\code{R}}{\eqn{n \times n}{n x n} observation (error) covariance matrix.}
#'  \item{\code{e}}{\eqn{T \times n}{T x n} estimates of observation errors \eqn{\textbf{e}_t}{e(t)}. Only available if \code{idio.ar1 = TRUE}.}
#'  \item{\code{rho}}{\eqn{n \times 1}{n x 1} estimates of AR(1) coefficients (\eqn{\rho}{rho}) in observation errors: \eqn{e_t = \rho e_{t-1} + v_t}{e(t) = rho e(t-1) + v(t)}. Only available if \code{idio.ar1 = TRUE}.}
#'  \item{\code{loglik}}{vector of log-likelihoods - one for each EM iteration. The final value corresponds to the log-likelihood of the reported model.}
#'  \item{\code{tol}}{The numeric convergence tolerance used.}
#'  \item{\code{converged}}{single logical valued indicating whether the EM algorithm converged (within \code{max.iter} iterations subject to \code{tol}).}
#'  \item{\code{anyNA}}{single logical valued indicating whether there were any (internal) missing values in the data (determined after removal of rows with too many missing values). If \code{FALSE}, \code{X_imp} is simply the original data in matrix form, and does not have the \code{"missing"} attribute attached.}
#'  \item{\code{rm.rows}}{vector of any cases (rows) that were removed beforehand (subject to \code{max.missing} and \code{na.rm.method}). If no cases were removed the slot is \code{NULL}. }
#'  \item{\code{quarterly.vars}}{names of the quarterly variables (if any).}
#'  \item{\code{em.method}}{The EM method used.}
#'  \item{\code{call}}{call object obtained from \code{match.call()}.}
#'
#' @references
#' Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large approximate dynamic factor models based on Kalman filtering. \emph{Journal of Econometrics, 164}(1), 188-205.
#'
#' Doz, C., Giannone, D., & Reichlin, L. (2012). A quasi-maximum likelihood approach for large, approximate dynamic factor models. \emph{Review of Economics and Statistics, 94}(4), 1014-1024.
#'
#' Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor models on datasets with arbitrary pattern of missing data. \emph{Journal of Applied Econometrics, 29}(1), 133-160.
#'
#' Stock, J. H., & Watson, M. W. (2016). Dynamic Factor Models, Factor-Augmented Vector Autoregressions, and Structural Vector Autoregressions in Macroeconomics. \emph{Handbook of Macroeconomics, 2}, 415–525. https://doi.org/10.1016/bs.hesmac.2016.04.002
#'
#' @seealso \link{dfms-package}
#'
#' @useDynLib dfms, .registration = TRUE
#' @importFrom collapse fscale qsu fvar flag fmedian fmedian.default qM unattrib na_omit %=% %-=% %+=% %/=% %*=% %r*% whichv whichNA vec setColnames setDimnames
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
#' BM14 <- merge(BM14_M, BM14_Q)
#' BM14[, BM14_Models$log_trans] %<>% log()
#' BM14[, BM14_Models$freq == "M"] %<>% diff()
#' BM14[, BM14_Models$freq == "Q"] %<>% diff(3)
#'
#' \donttest{
#' ### Small Model ---------------------------------------
#'
#' # IC for number of factors
#' IC_small <- ICr(BM14[, BM14_Models$small], max.r = 5)
#' plot(IC_small)
#' screeplot(IC_small)
#'
#' # I take 2 factors. Now number of lags
#' VARselect(IC_small$F_pca[, 1:2])
#'
#' # Estimating the model with 2 factors and 3 lags
#' dfm_small <- DFM(BM14[, BM14_Models$small], r = 2, p = 3,
#'     quarterly.vars = BM14_Models %$% series[freq == "Q" & small])
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
#' IC_medium <- ICr(BM14[, BM14_Models$medium])
#' plot(IC_medium)
#' screeplot(IC_medium)
#'
#' # I take 3 factors. Now number of lags
#' VARselect(IC_medium$F_pca[, 1:3])
#'
#' # Estimating the model with 3 factors and 3 lags
#' dfm_medium <- DFM(BM14[, BM14_Models$medium], r = 3, p = 3,
#'     quarterly.vars = BM14_Models %$% series[freq == "Q" & medium])
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
#'
#' # IC for number of factors
#' IC_large <- ICr(BM14)
#' plot(IC_large)
#' screeplot(IC_large)
#'
#' # I take 6 factors. Now number of lags
#' VARselect(IC_large$F_pca[, 1:6])
#'
#' # Estimating the model with 6 factors and 3 lags
#' dfm_large <- DFM(BM14, r = 6, p = 3,
#'     quarterly.vars = BM14_Models %$% series[freq == "Q"])
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
                idio.ar1 = FALSE,
                quarterly.vars = NULL, # blocks = c(),
                rQ = c("none", "diagonal", "identity"),
                rR = c("diagonal", "identity", "none"),
                em.method = c("auto", "DGR", "BM", "none"),
                min.iter = 25L,
                max.iter = 100L,
                tol = 1e-4,
                pos.corr = TRUE,
                check.increased = FALSE) {

  # Strict checking of inputs: as demanded by rOpenSci
  #' @srrstats {G2.1} *Implement assertions on types of inputs (see the initial point on nomenclature above).*
  #' @srrstats {G2.1a} *Provide explicit secondary documentation of expectations on data types of all vector inputs.*
  #' @srrstats {G2.2} *Appropriately prohibit or restrict submission of multivariate input to parameters expected to be univariate.*
  #' @srrstats {G2.3} *For univariate character input:*
  #' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
  #' @srrstats {G2.3b} *Either: use `tolower()` or equivalent to ensure input of character parameters is not case dependent; or explicitly document that parameters are strictly case-sensitive.*
  #' @srrstats {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
  #' @srrstats {G2.4a} *explicit conversion to `integer` via `as.integer()`*
  #' @srrstats {G2.4b} *explicit conversion to continuous via `as.numeric()`*
  #' @srrstats {G2.4c} *explicit conversion to character via `as.character()` (and not `paste` or `paste0`)*
  #' @srrstats {G2.4d} *explicit conversion to factor via `as.factor()`*
  #' @srrstats {G2.4e} *explicit conversion from factor via `as...()` functions*
  #' @srrstats {G2.5} *Where inputs are expected to be of `factor` type, secondary documentation should explicitly state whether these should be `ordered` or not, and those inputs should provide appropriate error or other routines to ensure inputs follow these expectations.*
  #' @srrstats {G2.6} *Software which accepts one-dimensional input should ensure values are appropriately pre-processed regardless of class structures.*
  #' @srrstats {G2.7} *Software should accept as input as many of the above standard tabular forms as possible, including extension to domain-specific forms.*
  #' @srrstats {G2.8} *Software should provide appropriate conversion or dispatch routines as part of initial pre-processing to ensure that all other sub-functions of a package receive inputs of a single defined class or type.*
  #' @srrstats {G2.9} *Software should issue diagnostic messages for type conversion in which information is lost (such as conversion of variables from factor to character; standardisation of variable names; or removal of meta-data such as those associated with [`sf`-format](https://r-spatial.github.io/sf/) data) or added (such as insertion of variable or column names where none were provided).*
  #' @srrstats {G2.10} *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.*
  #' @srrstats {G2.11} *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*
  #' @srrstats {G2.12} *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.*
  #' -> In general, all input to DFM() is converted to matrix by collapse::qM(), which also removes classes from objects that are already matrix based. All internal code is based on plain matrices.
  rRi <- switch(tolower(rR[1L]), identity = 0L, diagonal = 1L, none = 2L, stop("Unknown rR option:", rR[1L]))
  rQi <- switch(tolower(rQ[1L]), identity = 0L, diagonal = 1L, none = 2L, stop("Unknown rQ option:", rQ[1L]))
  #' @srrstats {G2.0} *Implement assertions on lengths of inputs, particularly through asserting that inputs expected to be single- or multi-valued are indeed so.*
  #' @srrstats {G2.0a} Provide explicit secondary documentation of any expectations on lengths of inputs
  if(sum(length(r), length(p), length(idio.ar1), length(min.iter), length(max.iter), length(tol), length(pos.corr), length(check.increased)) != 8L)
    stop("Parameters r, p, idio.ar1, min.iter, max.iter, tol, pos.corr and check.increased need to be length 1")
  if(!is.numeric(r) || r <= 0L) stop("r needs to be integer > 0")
  if(!is.integer(r)) r <- as.integer(r)
  if(!is.numeric(p) || p <= 0L) stop("p needs to be integer > 0")
  if(!is.integer(p)) p <- as.integer(p)
  if(!is.numeric(min.iter) || min.iter < 0L) stop("min.iter needs to be integer >= 0")
  if(!is.integer(min.iter)) min.iter <- as.integer(min.iter)
  if(!is.numeric(max.iter) || max.iter < min.iter) stop("max.iter needs to be integer >= min.iter")
  if(!is.integer(max.iter)) max.iter <- as.integer(max.iter)
  if(!is.numeric(tol) || tol <= 0) stop("tol needs to be numeric > 0")
  if(!is.logical(idio.ar1) || is.na(idio.ar1)) stop("idio.ar1 needs to be logical")
  if(!is.logical(pos.corr) || is.na(pos.corr)) stop("pos.corr needs to be logical")
  if(!is.logical(check.increased) || is.na(check.increased)) stop("check.increased needs to be logical")
  if(!is.null(quarterly.vars) && !is.character(quarterly.vars)) stop("quarterly.vars needs to be a character vector with the names of quarterly variables in X")

  #' @srrstats {G5.4a} *For new methods, it can be difficult to separate out correctness of the method from the correctness of the implementation, as there may not be reference for comparison. In this case, testing may be implemented against simple, trivial cases or against multiple implementations such as an initial R implementation compared with results from a C/C++ implementation.*
  #' @srrstats {G5.4b} *For new implementations of existing methods, correctness tests should include tests against previous implementations. Such testing may explicitly call those implementations in testing, preferably from fixed-versions of other software, or use stored outputs from those where that is not possible.*
  #' @srrstats {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original implementations is not available*
  #' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
  #' @srrstats {G5.6} **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
  #' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
  #' @srrstats {G5.6b} *Parameter recovery tests should be run with multiple random seeds when either data simulation or the algorithm contains a random component. (When long-running, such tests may be part of an extended, rather than regular, test suite; see G4.10-4.12, below).*
  #' @srrstats {G5.7} **Algorithm performance tests** *to test that implementation performs as expected as properties of data change. For instance, a test may show that parameters approach correct estimates within tolerance as data size increases, or that convergence times decrease for higher convergence thresholds.*
  #' @srrstats {G5.8} **Edge condition tests** *to test that these conditions produce expected behaviour such as clear warnings or errors when confronted with data with extreme properties including but not limited to:*
  #' @srrstats {G5.8a} *Zero-length data*
  #' @srrstats {G5.8b} *Data of unsupported types (e.g., character or complex numbers in for functions designed only for numeric data)*
  #' @srrstats {G5.8c} *Data with all-`NA` fields or columns or all identical fields or columns*
  #' @srrstats {G5.8d} *Data outside the scope of the algorithm (for example, data with more fields (columns) than observations (rows) for some regression algorithms)*
  #' @srrstats {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
  #' @srrstats {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
  #' @srrstats {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results*
  #' @srrstats {G5.10} *Extended tests should included and run under a common framework with other tests but be switched on by flags such as as a `<MYPKG>_EXTENDED_TESTS="true"` environment variable.* - The extended tests can be then run automatically by GitHub Actions for example by adding the following to the `env` section of the workflow:
  #' @srrstats {G5.11} *Where extended tests require large data sets or other assets, these should be provided for downloading and fetched as part of the testing workflow.*
  #' @srrstats {G5.11a} *When any downloads of additional data necessary for extended tests fail, the tests themselves should not fail, rather be skipped and implicitly succeed with an appropriate diagnostic message.*
  #' @srrstats {G5.12} *Any conditions necessary to run extended tests such as platform requirements, memory, expected runtime, and artefacts produced that may need manual inspection, should be described in developer documentation such as a `CONTRIBUTING.md` or `tests/README.md` file.*
  #' @srrstats {UL1.0} *Unsupervised Learning Software should explicitly document expected format (types or classes) for input data, including descriptions of types or classes which are not accepted; for example, specification that software accepts only numeric inputs in `vector` or `matrix` form, or that all inputs must be in `data.frame` form with both column and row names.*
  #' @srrstats {UL1.1} *Unsupervised Learning Software should provide distinct sub-routines to assert that all input data is of the expected form, and issue informative error messages when incompatible data are submitted.*
  #' @srrstats {UL1.2} *Unsupervised learning which uses row or column names to label output objects should assert that input data have non-default row or column names, and issue an informative message when these are not provided.*
  #' @srrstats {UL1.3} *Unsupervised Learning Software should transfer all relevant aspects of input data, notably including row and column names, and potentially information from other `attributes()`, to corresponding aspects of return objects.*
  #' @srrstats {UL1.3a} *Where otherwise relevant information is not transferred, this should be explicitly documented.*
  #' @srrstats {UL1.4} *Unsupervised Learning Software should document any assumptions made with regard to input data; for example assumptions about distributional forms or locations (such as that data are centred or on approximately equivalent distributional scales). Implications of violations of these assumptions should be both documented and tested, in particular:*
  #' @srrstats {UL1.4a} *Software which responds qualitatively differently to input data which has components on markedly different scales should explicitly document such differences, and implications of submitting such data.*
  #' -> Overall, I think the software is well behaved. It is not possible to use it with unsupported types (e.g. character or complex) or zero length data, it reliably produces the same output from the same inputs, it is meant to deal with noise (that's what factor models are for after all)
  #' @srrstats {UL2.0} *Routines likely to give unreliable or irreproducible results in response to violations of assumptions regarding input data (see UL1.6) should implement pre-processing steps to diagnose potential violations, and issue appropriately informative messages, and/or include parameters to enable suitable transformations to be applied.*
  #' @srrstats {UL2.1} *Unsupervised Learning Software should document any transformations applied to input data, for example conversion of label-values to `factor`, and should provide ways to explicitly avoid any default transformations (with error or warning conditions where appropriate).*
  #' @srrstats {UL2.2} *Unsupervised Learning Software which accepts missing values in input data should implement explicit parameters controlling the processing of missing values, ideally distinguishing `NA` or `NaN` values from `Inf` values.*
  #' @srrstats {UL2.3} *Unsupervised Learning Software should implement pre-processing routines to identify whether aspects of input data are perfectly collinear.*
  #' -> It does not matter for DFMs if series are collinear. Users may deliberately choose to duplicate quarterly series to increase their weight in mixed-frequency estimations.
  #' @srrstats {UL3.0} *Algorithms which apply sequential labels to input data (such as clustering or partitioning algorithms) should ensure that the sequence follows decreasing group sizes (so labels of "1", "a", or "A" describe the largest group, "2", "b", or "B" the second largest, and so on.)*
  #' @srrstats {UL3.1} *Dimensionality reduction or equivalent algorithms which label dimensions should ensure that that sequences of labels follows decreasing "importance" (for example, eigenvalues or variance contributions).*
  #' @srrstats {UL3.2} *Unsupervised Learning Software for which input data does not generally include labels (such as `array`-like data with no row names) should provide an additional parameter to enable cases to be labelled.*
  #' @srrstats {UL4.0} *Unsupervised Learning Software should return some form of "model" object, generally through using or modifying existing class structures for model objects, or creating a new class of model objects.*
  #' @srrstats {UL4.2} *The return object from Unsupervised Learning Software should include, or otherwise enable immediate extraction of, all parameters used to control the algorithm used.*
  #' @srrstats {UL4.3} *Model objects returned by Unsupervised Learning Software should implement or appropriately extend a default `print` method which provides an on-screen summary of model (input) parameters and methods used to generate results. The `print` method may also summarise statistical aspects of the output data or results.*
  #' @srrstats {UL4.3a} *The default `print` method should always ensure only a restricted number of rows of any result matrices or equivalent are printed to the screen.*
  #' @srrstats {UL4.4} *Unsupervised Learning Software should also implement `summary` methods for model objects which should summarise the primary statistics used in generating the model (such as numbers of observations, parameters of methods applied). The `summary` method may also provide summary statistics from the resultant model.*
  #' @srrstats {UL6.0} *Objects returned by Unsupervised Learning Software should have default `plot` methods, either through explicit implementation, extension of methods for existing model objects, through ensuring default methods work appropriately, or through explicit reference to helper packages such as [`factoextra`](https://github.com/kassambara/factoextra) and associated functions.*
  #' @srrstats {UL6.2} *Where default plot methods include labelling components of return objects (such as cluster labels), routines should ensure that labels are automatically placed to ensure readability, and/or that appropriate diagnostic messages are issued where readability is likely to be compromised (for example, through attempting to place too many labels).*
  #' @srrstats {UL7.0} *Inappropriate types of input data are rejected with expected error messages.*
  #' @srrstats {UL7.1} *Tests should demonstrate that violations of assumed input properties yield unreliable or invalid outputs, and should clarify how such unreliability or invalidity is manifest through the properties of returned objects.*
  #' @srrstats {UL7.2} *Demonstrate that labels placed on output data follow decreasing group sizes (**UL3.0**)*
  #' @srrstats {UL7.3} *Demonstrate that labels on input data are propagated to, or may be recovered from, output data.
  #' @srrstats {UL7.4} *Demonstrate that submission of new data to a previously fitted model can generate results more efficiently than initial model fitting.*
  #' @srrstats {TS1.0} *Time Series Software should use and rely on explicit class systems developed for representing time series data, and should not permit generic, non-time-series input*
  #' @srrstats {TS1.1} *Time Series Software should explicitly document the types and classes of input data able to be passed to each function.*
  #' @srrstats {TS1.2} *Time Series Software should implement validation routines to confirm that inputs are of acceptable classes (or represented in otherwise appropriate ways for software which does not use class systems).*
  #' @srrstats {TS1.3} *Time Series Software should implement a single pre-processing routine to validate input data, and to appropriately transform it to a single uniform type to be passed to all subsequent data-processing functions (the [`tsbox` package](https://www.tsbox.help/) provides one convenient approach for this).*
  #' @srrstats {TS1.4} *The pre-processing function described above should maintain all time- or date-based components or attributes of input data.*
  #' @srrstats {TS1.5} *The software should ensure strict ordering of the time, frequency, or equivalent ordering index variable.*
  #' @srrstats {TS1.6} *Any violations of ordering should be caught in the pre-processing stages of all functions.*
  #' -> This software is at the intersection of dimensionality reduction and time series, and does not require input objects to have a certain class.
  #' For all practical purposes I believe it is convenient to not require certain object types, although certain input object attributes such as the
  #' time scale of 'ts' or 'xts' objects could be used in model plots. This is an area where I still want to improve the software e.g. keeping it class
  #' agnostic but using some of the information contained in time series objects passed as inputs.
  #' @srrstats {TS1.7} *Accept inputs defined via the [`units` package](https://github.com/r-quantities/units/) for attributing SI units to R vectors.*
  #' @srrstats {TS1.8} *Where time intervals or periods may be days or months, be explicit about the system used to represent such, particularly regarding whether a calendar system is used, or whether a year is presumed to have 365 days, 365.2422 days, or some other value.*
  #' @srrstats {TS2.0} *Time Series Software which presumes or requires regular data should only allow **explicit** missing values, and should issue appropriate diagnostic messages, potentially including errors, in response to any **implicit** missing values.*
  #' @srrstats {TS2.1} *Where possible, all functions should provide options for users to specify how to handle missing data, with options minimally including:*
  #' @srrstats {TS2.1a} *error on missing data; or.
  #' @srrstats {TS2.1b} *warn or ignore missing data, and proceed to analyse irregular data, ensuring that results from function calls with regular yet missing data return identical values to submitting equivalent irregular data with no missing values; or*
  #' @srrstats {TS2.1c} *replace missing data with appropriately imputed values.*
  #' @srrstats {TS2.2} *Consider stationarity of all relevant moments - typically first (mean) and second (variance) order, or otherwise document why such consideration may be restricted to lower orders only.*
  #' @srrstats {TS2.3} *Explicitly document all assumptions and/or requirements of stationarity*
  #' @srrstats {TS2.4} *Implement appropriate checks for all relevant forms of stationarity, and either:*
  #' @srrstats {TS2.4a} *issue diagnostic messages or warnings; or*
  #' @srrstats {TS2.4b} *enable or advise on appropriate transformations to ensure stationarity.*
  #' -> it is possible to run DFM() on non-stationary data and obtain convergent results, although this is strongly discouraged.
  rp <- r * p
  sr <- seq_len(r)
  fnam <- paste0("f", sr)
  unam <- paste0("u", sr)
  ax <- attributes(X)
  ilX <- is.list(X)
  Xstat <- qsu(X)
  #' @srrstats {UL1.4b} *Examples or other documentation should not use `scale()` or equivalent transformations without explaining why scale is applied, and explicitly illustrating and contrasting the consequences of not applying such transformations.*
  #' -> Data is always scaled and centered before estimation to remove the need for intercept terms in the Kalman Filter and Smoother. This is stated in the documentatio e.g. ?DFM
  X <- fscale(qM(X), na.rm = TRUE)
  Xnam <- dimnames(X)[[2L]]
  dimnames(X) <- NULL
  n <- ncol(X)
  if(MFl <- !is.null(quarterly.vars)) {
    qind <- ckmatch(quarterly.vars, Xnam)
    nq <- length(qind)
    if(!identical(sort(qind), (n-nq+1):n)) stop("Pleae order your data such that quarterly variables are to the right (after) the monthly variables")
  }

  # Missing values
  #' @srrstats {G2.13} *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
  #' @srrstats {G2.14} *Where possible, all functions should provide options for users to specify how to handle missing (`NA`) data, with options minimally including:*
  #' @srrstats {G2.14a} *error on missing data*
  #' @srrstats {G2.14b} *ignore missing data with default warnings or messages issued*
  #' @srrstats {G2.14c} *replace missing data with appropriately imputed values*
  #' @srrstats {G2.15} *Functions should never assume non-missingness, and should never pass data with potential missing values to any base routines with default `na.rm = FALSE`-type parameters (such as [`mean()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/mean.html), [`sd()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/sd.html) or [`cor()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html)).*
  #' -> DFM() estimation assumes missing data, all routines are designed for missing data.
  #' @srrstats {G2.16} *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.*
  #' -> NaN will be treated like NA, and Inf and -Inf will be handled as any other numeric values (in line with nearly all other software, including base R). I initially call anyNA() to efficiently check for missingness, if FALSE, the software will assume complete cases.
  anymiss <- anyNA(X)
  if(anymiss) { # Missing value removal / imputation
    # TODO: Should the data be scaled and centered again after missing value imputation?
    # I think no because we just use it to initialize the filter, and the filter runs on the standardized data with missing values
    X_imp <- tsnarmimp(X, ...)
    W <- attr(X_imp, "missing")
    rm.rows <- attr(X_imp, "rm.rows")
    attributes(X_imp) <- list(dim = dim(X_imp))
    if(length(rm.rows)) X <- X[-rm.rows, , drop = FALSE]
  } else {
    X_imp <- X
    rm.rows <- NULL
  }
  TT <- nrow(X)

  # This is because after removing missing rows, the data could be complete, e.g. differencing data with diff.xts() of collapse::fdiff() just gives a NA row
  if(anymiss && length(rm.rows)) anymiss <- any(W)
  BMl <- switch(tolower(em.method[1L]), auto = anymiss || idio.ar1 || MFl, dgr = FALSE || idio.ar1 || MFl,
                bm = TRUE, none = NA, stop("Unknown EM option:", em.method[1L]))

  # Run PCA to get initial factor estimates:
  #' @srrstats {G3.0} *Statistical software should never compare floating point numbers for equality. All numeric equality comparisons should either ensure that they are made between integers, or use appropriate tolerances for approximate equality.*
  #' @srrstats {G3.1} *Statistical software which relies on covariance calculations should enable users to choose between different algorithms for calculating covariances, and should not rely solely on covariances from the `stats::cov` function.*
  #' @srrstats {G3.1a} *The ability to use arbitrarily specified covariance methods should be documented (typically in examples or vignettes).*
  #' -> All PCA routines I know of, including the authors original Matlab code, use Pearsons Covariance for the eigen decomposition, and I currently see no reason to provide alternatives in *dfms*.
  # v <- svd(X_imp, nu = 0L, nv = min(as.integer(r), n, TT))$v # Not faster than eigen...
  eigen_decomp <- eigen(cov(X_imp), symmetric = TRUE)
  # TODO: better way to ensure factors correlate positively with data?
  if(pos.corr) {
    PCS <- X_imp %*% eigen_decomp$vectors
    setop(eigen_decomp$vectors, "*", c(-1,1)[(colSums(PCS %*=% rowMeans(X_imp)) > 0) + 1L], rowwise = TRUE)
  }
  v <- eigen_decomp$vectors[, sr, drop = FALSE]
  F_pc <- X_imp %*% v

  ## Initial System Matrices
  if(MFl) {
    init <- if(idio.ar1) stop("Mixed frequency with autocorrelated errors is not yet implemented") else
      init_cond_MQ(X, X_imp, F_pc, v, n, r, p, TT, nq, rRi, rQi)
  } else {
    init <- if(idio.ar1) init_cond_idio_ar1(X, F_pc, v, n, r, p, BMl, rRi, rQi, anymiss, tol) else
      init_cond(X, F_pc, v, n, r, p, BMl, rRi, rQi)
  }
  A <- init$A; C <- init$C; Q <- init$Q; R <- init$R; F_0 <- init$F_0; P_0 <- init$P_0

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
                       F_2s = setColnames(F_kal, fnam),
                       P_2s = setDimnames(kfs_res$P_smooth[sr, sr,, drop = FALSE], list(fnam, fnam, NULL)),
                       anyNA = anymiss, # || length(rm.rows), # This is for internal missing values only
                       rm.rows = rm.rows,
                       quarterly.vars = quarterly.vars,
                       em.method = if(is.na(BMl)) "none" else if(BMl) "BM" else "DGR",
                       call = match.call())

  # em.method = "none": only report two-step solution
  if(is.na(BMl)) {
    # Better solution for system matrix estimation after Kalman Filtering and Smoothing?
    var <- .VAR(F_kal, p)
    beta <- ainv(crossprod(F_kal)) %*% crossprod(F_kal, if(anymiss) replace(X_imp, W, 0) else X_imp) # good??
    Q <- switch(rQi + 1L, diag(r),  diag(fvar(var$res)), cov(var$res))
    if(rRi) {
      if(idio.ar1) { # If AR(1) residuals, need to estimate coefficient and clean residuals from autocorrelation
        e <- kfs_res$F_smooth[, -seq_len(rp), drop = FALSE]
        res_AC1 <- AC1(e, FALSE)
        res <- e[-1L, ] %-=% setop(e[-nrow(e), ], "*", res_AC1, rowwise = TRUE)
      } else res <- X - F_kal %*% beta
      R <- if(rRi == 2L) cov(res, use = "pairwise.complete.obs") else diag(fvar(res, na.rm = TRUE))
    } else R <- diag(n)
    final_object <- c(object_init[1:6],
                      list(A = setDimnames(t(var$A), lagnam(fnam, p)), # A[sr, , drop = FALSE],
                           C = setDimnames(t(beta), list(Xnam, fnam)), # C[, sr, drop = FALSE],
                           Q = setDimnames(Q, list(unam, unam)),       # Q[sr, sr, drop = FALSE],
                           R = setDimnames(R, list(Xnam, Xnam)),
                      e = if(idio.ar1) setColnames(e, Xnam) else NULL,
                      rho = if(idio.ar1) setNames(res_AC1, Xnam) else NULL),
                      object_init[-(1:6)])
    class(final_object) <- "dfm"
    return(final_object)
  }

  previous_loglik <- -Inf # .Machine$double.xmax
  loglik_all <- integer(max.iter)
  num_iter <- 0L
  converged <- FALSE

  if(BMl) {
    if(MFl) {
      nm <- n - nq
      rpC <- r * max(p, 5L)
      rpC1nq <- (rpC+1L):(rpC+nq)
      NW <- !W
      R_mat <- kronecker(Rcon, diag(r)) # TODO: autsource this
      expr <- if(idio.ar1) .EM_BM_MQ_idio else .EM_BM_MQ
      dnkron <- matrix(1, r, r) %x% diag(nm) # Used to be inside EMstep, taken out to speed up the algorithm
      dnkron_ind <- whichv(dnkron, 1)
    } else {
      expr <- if(idio.ar1) .EM_BM_idio else .EM_BM
      dnkron <- matrix(1, r, r) %x% diag(n) # Used to be inside EMstep, taken out to speed up the algorithm
      dnkron_ind <- whichv(dnkron, 1)
    }
    dgind <- 0:(n-1) * n + 1:n
    XW0 <- X_imp
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
    if(num_iter < min.iter) converged <- FALSE else {
      converged <- em_converged(loglik, previous_loglik, tol, check.increased)
      if(check.increased) converged <- converged[1L] && !converged[2L]
    }

    previous_loglik <- loglik
    num_iter <- num_iter + 1L
    loglik_all[num_iter] <- loglik

  }

  if(converged) message("Converged after ", num_iter, " iterations.")
  else warning("Maximum number of iterations reached.")

  ## Run the Kalman filtering and smoothing step for the last time
  ## with optimal estimates
  kfs_res <- eval(.KFS, em_res, encl)

  final_object <- c(object_init[1:6],
               list(F_qml = setColnames(kfs_res$F_smooth[, sr, drop = FALSE], fnam),
                    P_qml = setDimnames(kfs_res$P_smooth[sr, sr,, drop = FALSE], list(fnam, fnam, NULL)),
                    A = setDimnames(em_res$A[sr, seq_len(rp), drop = FALSE], lagnam(fnam, p)),
                    C = setDimnames(em_res$C[, sr, drop = FALSE], list(Xnam, fnam)),
                    Q = setDimnames(em_res$Q[sr, sr, drop = FALSE], list(unam, unam)),
                    R = setDimnames(if(MFl && idio.ar1) stop("not implemented yet") else if(MFl)
                                    `[<-`(em_res$R, (nm+1):n, (nm+1):n, value = em_res$Q[rpC1nq, rpC1nq]) else if(idio.ar1)
                                    em_res$Q[-seq_len(rp), -seq_len(rp), drop = FALSE] else em_res$R, list(Xnam, Xnam)),
                    e = if(idio.ar1) setColnames(kfs_res$F_smooth[, -seq_len(rp), drop = FALSE], Xnam) else NULL,
                    rho = if(idio.ar1) setNames(diag(em_res$A)[-seq_len(rp)], Xnam) else NULL,
                    loglik = if(num_iter == max.iter) loglik_all else loglik_all[seq_len(num_iter)],
                    tol = tol,
                    converged = converged),
                    object_init[-(1:6)])

  class(final_object) <- "dfm"
  return(final_object)
}

