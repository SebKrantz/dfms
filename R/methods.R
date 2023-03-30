#' @name summary.dfm
#' @aliases print.dfm
#' @aliases summary.dfm
#' @aliases print.dfm_summary
#'
#' @title DFM Summary Methods
#'
#' @description Summary and print methods for class 'dfm'. \code{print.dfm} just prints basic model information and the factor transition matrix \eqn{\textbf{A}}{A},
#' \code{summary.dfm} returns all system matrices and additional residual and goodness of fit statistics - with a print method allowing full or compact printout.
#'
#' @param x,object an object class 'dfm'.
#' @param digits integer. The number of digits to print out.
#' @param \dots not used.
#' @importFrom collapse qsu frange
#' @export
print.dfm <- function(x, digits = 4L, ...) {

  X <- x$X_imp
  A <- x$A
  r <- dim(A)[1L]
  p <- dim(A)[2L]/r
  cat("Dynamic Factor Model: n = ", dim(X)[2L], ", T = ", dim(X)[1L], ", r = ", r, ", p = ", p, ", %NA = ",
      if(x$anyNA) round(sum(attr(X, "missing"))/prod(dim(X))*100, digits) else 0,"\n", sep = "")
  if(length(x$rho)) cat("   with AR(1) errors: mean(abs(rho)) =", round(mean(abs(x$rho)), 3), "\n")
  fnam <- paste0("f", seq_len(r))
  cat("\nFactor Transition Matrix [A]\n")
  print(round(A, digits))
}

#' @rdname summary.dfm
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"2s"} or \code{"pca"}.
#' @param \dots not used.
#' @return Summary information following a dynamic factor model estimation.
#' @importFrom stats cov
#' @importFrom collapse pwcov
#' @export
summary.dfm <- function(object, method = switch(object$em.method, none = "2s", "qml"), ...) {

  X <- object$X_imp
  Fa <- switch(tolower(method), pca = object$F_pca, `2s` = object$F_2s, qml = object$F_qml, stop("Unkown method", method))
  A <- object$A
  r <- dim(A)[1L]
  p <- dim(A)[2L] / r
  C <- object$C
  idio_ar1 <- !is.null(object$rho)
  res <- if(idio_ar1) object[["e"]] else  X - tcrossprod(Fa, C)
  anymissing <- object$anyNA
  if(!idio_ar1 && anymissing) res[attr(X, "missing")] <- NA
  rescov <- pwcov(res, use = if(!idio_ar1 && anymissing) "pairwise.complete.obs" else "everything", P = TRUE)
  ACF <- if(idio_ar1) object$rho else AC1(res, anymissing)
  R2 <- 1 - diag(rescov[,, 1L])
  summ <- list(info = c(n = dim(X)[2L], T = dim(X)[1L], r = r, p = p,
                        `%NA` = if(anymissing) sum(attr(X, "missing")) / prod(dim(X)) * 100 else 0),
               call = object$call,
               idio_ar1 = idio_ar1,
               F_stats = msum(Fa),
               A = A,
               F_cov = pwcov(Fa, P = TRUE),
               Q = object$Q,
               C = C,
               R_diag = diag(object$R),
               res_cov = rescov,
               res_ACF = ACF,
               res_ACF_stats = msum(ACF),
               R2 = R2,
               R2_stats = msum(R2))
  class(summ) <- "dfm_summary"
  return(summ)
}

#' @rdname summary.dfm
#' @param compact integer. Display a more compact printout: \code{0} prints everything, \code{1} omits the observation matrix \eqn{\textbf{C}}{C} and residual covariance matrix \code{cov(resid(model))}, and \code{2} omits all disaggregated information on the input data. Sensible default are chosen for different sizes of the input dataset so as to limit large printouts.
#' @param \dots not used.
#'
#' @examples
#' mod = DFM(diff(BM14_Q), 2, 3)
#' print(mod)
#' summary(mod)
#'
#' @export
print.dfm_summary <- function(x,
                              digits = 4L,
                              compact = sum(x$info["n"] > 15, x$info["n"] > 40), ...) {

  inf <- as.integer(x$info[1:4])
  cat("Dynamic Factor Model: n = ", inf[1L], ", T = ", inf[2L], ", r = ", inf[3L], ", p = ", inf[4L],
      ", %NA = ", round(x$info[5L], digits), "\n", sep = "")
  if(length(x$idio_ar1)) cat("   with AR(1) errors: mean(abs(rho)) =", round(mean(abs(x$res_ACF)), 3), "\n")
  cat("\nCall: ", deparse(x$call))
  # cat("\nModel: ", ))
  cat("\n\nSummary Statistics of Factors [F]\n")
  print(x$F_stats, digits)
  cat("\nFactor Transition Matrix [A]\n")
  print(x$A, digits = digits)
  cat("\nFactor Covariance Matrix [cov(F)]\n")
  print(x$F_cov, digits)
  cat("\nFactor Transition Error Covariance Matrix [Q]\n")
  print(round(x$Q, digits))
  if(compact == 0L) {
  cat("\nObservation Matrix [C]\n")
  print(round(x$C, digits))
  }
  if(compact < 2L) {
  cat("\nObservation Error Covariance Matrix [diag(R) - Restricted]\n")
  # cat("\n Estimated Diagonal (DFM Assumes R is Diagonal)\n")
  print(round(x$R_diag, digits))
  }
  if(compact == 0L) {
  cat("\nObservation Residual Covariance Matrix [cov(resid(DFM))]\n")
  print(x$res_cov, digits)
  }
  if(compact < 2L) {
  cat("\nResidual AR(1) Serial Correlation\n")
  print(x$res_ACF, digits) # TODO: Add P-Value
  }
  cat("\nSummary of Residual AR(1) Serial Correlations\n")
  print(x$res_ACF_stats, digits)
  if(compact < 2L) {
  cat("\nGoodness of Fit: R-Squared\n")
  print(x$R2, digits)
  }
  cat("\nSummary of Individual R-Squared's\n")
  print(x$R2_stats, digits)
}


#' Plot DFM
#' @param x an object class 'dfm'.
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"2s"}, \code{"pca"} or \code{"all"} to plot all estimates.
#' @param type character. The type of plot: \code{"joint"}, \code{"individual"} or \code{"residual"}.
#' @param scale.factors logical. Standardize factor estimates, this usually improves the plot since the factor estimates corresponding to the greatest PCA eigenvalues tend to have a greater variance than the data.
#' @param \dots for \code{plot.dfm}: further arguments to \code{\link{plot}}, \code{\link{ts.plot}}, or \code{\link{boxplot}}, depending on the \code{type} of plot. For \code{screeplot.dfm}: further arguments to \code{\link{screeplot.ICr}}.
#' @returns Nothing.
#' @examples \donttest{
#' # Fit DFM with 3 factors and 3 lags in the transition equation
#' mod = DFM(diff(BM14_M), r = 3, p = 3)
#' plot(mod)
#' plot(mod, type = "individual", method = "all")
#' plot(mod, type = "residual")
#' }
#'
#' @importFrom graphics boxplot axis box mtext plot.default
#' @importFrom collapse unlist2d ckmatch na_rm seq_row
#' @export
plot.dfm <- function(x,
                     method = switch(x$em.method, none = "2s", "qml"),
                     type = c("joint", "individual", "residual"),
                     scale.factors = TRUE, ...) {
  Fa <- switch(tolower(method[1L]),
              all = cbind(x$F_pca, setColnames(x$F_2s, paste("TwoStep", colnames(x$F_2s))),
                          if(length(x$F_qml)) setColnames(x$F_qml, paste("QML", colnames(x$F_qml))) else NULL),
              pca = x$F_pca, `2s` = x$F_2s, qml = x$F_qml, stop("Unknown method:", method[1L]))

  nf <- dim(Fa)[2L]
  allests <- tolower(method[1L]) == "all"
  dots <- list(...)

  switch(tolower(type[1L]),
    joint = {
      Xr <- frange(x$X_imp)
      if(scale.factors) Fa <- fscale(Fa)
      Fr <- frange(Fa)
      ts.plot(x$X_imp, col = "grey85", ylim = c(min(Xr[1L], Fr[1L]), max(Xr[2L], Fr[2L])),
              ylab = if(is.null(dots$ylab)) "Value" else dots$ylab,
              main = if(is.null(dots$main)) "Standardized Series and Factor Estimates" else dots$main, ...)
      cols <- rainbow(nf)
      for (i in seq_len(nf)) lines(Fa[, i], col = cols[i])
      legend("topleft", colnames(Fa), col = cols, lty = 1, bty = "n", ncol = if(allests) 3L else 1L)
    },
    individual = {
      # if(allests) {
      if(scale.factors) Fa <- fscale(Fa)
      qml <- !is.null(x$F_qml)
      if(allests) nf <- nf / (2L + qml)

      # Extracted from plot.ts()...
      cex.lab = par("cex.lab")
      col.lab = par("col.lab")
      font.lab = par("font.lab")
      oldpar <- par(mar = c(0, 5.1, 0, 2.1), oma = c(6, 0, 5, 0), mfrow = c(nf, 1L))
      on.exit(par(oldpar))

      for(i in seq_len(nf)) {
        plot.default(Fa[, i], axes = FALSE, xlab = "", ylab = "", type = "n")
        lines(Fa[, i], type = 'l', col = if(allests) "red" else "black", ...)
        if(allests) {
          lines(Fa[, i + nf], type = 'l', col = "orange", ...)
          if(qml) lines(Fa[, i + 2L * nf], type = 'l', col = "blue", ...)
          if(i == 1L) legend("topleft", c("PCA", "TwoStep", if(qml) "QML"),
                             col = c("red", "orange", "blue"), lty = 1, bty = "n")
        }
        box(...)
        axis(2, xpd = NA, ...)
        if(i == nf) axis(1, xpd = NA, ...)
        mtext(paste("f", i), 2, line = 3, cex = cex.lab, col = col.lab, font = font.lab, ...)
        if(i == nf) mtext(if(is.null(dots$xlab)) "Time" else dots$xlab, side = 1, line = 3, cex = cex.lab, col = col.lab, font = font.lab, ...)
      }
      par(mfrow = c(1, 1)) # on.exit above takes care of changes to parameters.
      mtext(if(is.null(dots$main)) paste(if(scale.factors) "Standardized", "Factor Estimates") else dots$main,
            side = 3, line = 3, cex = par("cex.main"), font = par("font.main"), col = par("col.main"), ...)
      # } else {
      #   oldpar <- par(mfrow = c(nf, 1L))
      #   on.exit(par(oldpar))
      #   cnF <- colnames(Fa)
      #   for (i in seq_len(nf)) plot(Fa[, i], type = 'l', main = cnF[i], ylab = "Value",
      #                               xlab = if(i == nf) "Time" else "" , ...)
      # }
    },
    residual = {
      if(allests) stop("Need to choose a specific method for residual plots")
      oldpar <- par(mar = c(11.5, 4.1, 4.1, 2.1))
      on.exit(par(oldpar))
      boxplot(x$X_imp - tcrossprod(Fa, x$C), main = if(is.null(dots$main)) "Residuals by input variable" else dots$main, las = 2, ...)
    },
    stop("Unknown plot type: ", type[1L])
  )
}


#' Extract Factor Estimates in a Data Frame
#' @param x an object class 'dfm'.
#' @param method character. The factor estimates to use: any of \code{"qml"}, \code{"2s"}, \code{"pca"} (multiple can be supplied) or \code{"all"} for all estimates.
#' @param pivot character. The orientation of the frame: \code{"long"}, \code{"wide.factor"} or \code{"wide.method"}, \code{"wide"} or \code{"t.wide"}.
#' @param time a vector identifying the time dimension, or \code{NULL} to omit a time variable.
#' @param stringsAsFactors make factors from method and factor identifiers. Same as option to \code{\link{as.data.frame.table}}.
#' @param \dots not used.
#'
#' @return A data frame of factor estimates.
#'
#' @examples \donttest{
#' library(xts)
#' # Fit DFM with 3 factors and 3 lags in the transition equation
#' mod = DFM(diff(BM14_M), r = 3, p = 3)
#'
#' # Taking a single estimate:
#' print(head(as.data.frame(mod, method = "qml")))
#' print(head(as.data.frame(mod, method = "qml", pivot = "wide")))
#'
#' # Adding a proper time variable
#' time = index(BM14_M)[-1L]
#' print(head(as.data.frame(mod, method = "qml", time = time)))
#'
#' # All estimates: different pivoting methods
#' for (pv in c("long", "wide.factor", "wide.method", "wide", "t.wide")) {
#'    cat("\npivot = ", pv, "\n")
#'    print(head(as.data.frame(mod, pivot = pv, time = time), 3))
#' }
#' }
#'
#' @importFrom collapse ckmatch na_rm seq_row t_list unattrib
#' @importFrom stats setNames
#' @export
as.data.frame.dfm <- function(x, ...,
                              method = "all",
                              pivot = c("long", "wide.factor", "wide.method", "wide", "t.wide"),
                              time = seq_row(x$F_pca),
                              stringsAsFactors = TRUE) {

  estm <- c(PCA = "pca", TwoStep = "2s", QML = "qml")
  if(length(method) > 1L || method != "all")
     estm <- estm[ckmatch(method, estm, e = "Unknown method:")]
  estlist <- x[paste0("F_", estm)]
  names(estlist) <- names(estm)
  estlist <- na_rm(estlist) # Also removes NULL elements

  nam <- names(estlist)
  m <- length(estlist)
  TT <- nrow(estlist[[1L]])
  r <- ncol(estlist[[1L]])

  if(!is.null(time) && length(time) != TT) {
    if(length(x$rm.rows)) time <- time[-x$rm.rows]
    if(length(time) != TT) stop(sprintf("time must be a length %s vector or NULL", TT))
  }

  res <- switch(tolower(pivot[1L]),
    long = list(Method = if(stringsAsFactors) setAttrib(rep(1:m, each = TT*r), list(levels = nam, class = "factor")) else rep(nam, each = TT*r),
                Factor = if(stringsAsFactors) setAttrib(rep(1:r, times = m, each = TT), list(levels = paste0("f", 1:r), class = "factor")) else rep(paste0("f", 1:r), times = m, each = TT),
                Time = if(length(time)) rep(time, times = m*r) else NULL,
                Value = unlist(estlist, use.names = FALSE)),
    wide.factor = c(list(Method = if(stringsAsFactors) setAttrib(rep(1:m, each = TT), list(levels = nam, class = "factor")) else rep(nam, each = TT),
                         Time = if(length(time)) rep(time, times = m) else NULL),
                    setNames(lapply(t_list(unattrib(lapply(estlist, mctl))), unlist, FALSE, FALSE), paste0("f", 1:r))),
    wide.method = c(list(Factor = if(stringsAsFactors) setAttrib(rep(1:r, each = TT), list(levels = paste0("f", 1:r), class = "factor")) else rep(paste0("f", 1:r), each = TT),
                         Time = if(length(time)) rep(time, times = r) else NULL),
                    lapply(estlist, unattrib)),
    # If only one method, do not do combine names e.g. "QML_f1"? -> most of the time people just want a simple frame like this...
    wide = c(list(Time = time), setNames(unlist(lapply(estlist, mctl), FALSE, FALSE), if(length(nam) == 1L) paste0("f", 1:r) else outer(paste0("f", 1:r), nam, paste, sep = "_"))),
    t.wide = c(list(Time = time), setNames(unlist(t_list(lapply(estlist, mctl)), FALSE, FALSE), if(length(nam) == 1L) paste0("f", 1:r) else t(outer(paste0("f", 1:r), nam, paste, sep = "_")))),
    stop("Unknown pivot option:", pivot[1L])
  )

  if(is.null(time)) res <- na_rm(res)
  attr(res, "methods") <- estm
  attr(res, "row.names") <- .set_row_names(length(res[[1L]]))
  class(res) <- "data.frame"
  return(res)
}

#' @name residuals.dfm
#' @aliases residuals.dfm
#' @aliases resid.dfm
#' @aliases fitted.dfm
#'
#' @title DFM Residuals and Fitted Values
#' @description The residuals \eqn{\textbf{e}_t = \textbf{x}_t - \textbf{C} \textbf{F}_t}{e(t) = x(t) - C F(t)} or fitted values \eqn{\textbf{C} \textbf{F}_t}{C F(t)} of the DFM observation equation.
#'
#' @param object an object of class 'dfm'.
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"2s"} or \code{"pca"}.
#' @param orig.format logical. \code{TRUE} returns residuals/fitted values in a data format similar to \code{X}.
#' @param standardized logical. \code{FALSE} will put residuals/fitted values on the original data scale.
#' @param na.keep logical. \code{TRUE} inserts missing values where \code{X} is missing (default \code{TRUE} as residuals/fitted values are only defined for observed data). \code{FALSE} returns the raw prediction, which can be used to interpolate data based on the DFM. For residuals, \code{FALSE} returns the difference between the prediction and the initial imputed version of \code{X} use for PCA to initialize the Kalman Filter.
#' @param \dots not used.
#'
#' @return A matrix of DFM residuals or fitted values. If \code{orig.format = TRUE} the format may be different, e.g. a data frame.
#'
#' @examples \donttest{
#' library(xts)
#' # Fit DFM with 3 factors and 3 lags in the transition equation
#' mod = DFM(diff(BM14_M), r = 3, p = 3)
#'
#' # Residuals
#' head(resid(mod))
#' plot(resid(mod, orig.format = TRUE)) # this is an xts object
#'
#' # Fitted values
#' head(fitted(mod))
#' head(fitted(mod, orig.format = TRUE)) # this is an xts object
#' }
#'
#' @importFrom collapse TRA.matrix mctl setAttrib pad
#' @export
residuals.dfm <- function(object,
                          method = switch(object$em.method, none = "2s", "qml"),
                          orig.format = FALSE,
                          standardized = FALSE,
                          na.keep = TRUE, ...) {
  X <- object$X_imp
  Fa <- switch(tolower(method),
              pca = object$F_pca, `2s` = object$F_2s, qml = object$F_qml,
              stop("Unkown method", method))
  if(!(standardized && length(object[["e"]]))) {
    X_pred <- tcrossprod(Fa, object$C)
    if(!standardized) {  # TODO: What if AR(1) resid available?
      stats <- attr(X, "stats")
      X_pred <- unscale(X_pred, stats)
      res <- unscale(X, stats) - X_pred
    } else res <- X - X_pred
  } else res <- object[["e"]]
  if(na.keep && object$anyNA) res[attr(X, "missing")] <- NA
  if(orig.format) {
    if(length(object$rm.rows)) res <- pad(res, object$rm.rows, method = "vpos")
    if(attr(X, "is.list")) res <- mctl(res)
    return(setAttrib(res, attr(X, "attributes")))
  }
  return(qM(res))
}

#' @rdname residuals.dfm
#' @export
fitted.dfm <- function(object,
                       method = switch(object$em.method, none = "2s", "qml"),
                       orig.format = FALSE,
                       standardized = FALSE,
                       na.keep = TRUE, ...) {
  X <- object$X_imp
  Fa <- switch(tolower(method),
              pca = object$F_pca, `2s` = object$F_2s, qml = object$F_qml,
              stop("Unkown method", method))
  res <- tcrossprod(Fa, object$C)
  if(!standardized) res <- unscale(res, attr(X, "stats"))
  if(na.keep && object$anyNA) res[attr(X, "missing")] <- NA
  if(orig.format) {
    if(length(object$rm.rows)) res <- pad(res, object$rm.rows, method = "vpos")
    if(attr(X, "is.list")) res <- mctl(res)
    return(setAttrib(res, attr(X, "attributes")))
  }
  return(qM(res))
}

#% @aliases forecast.dfm
#' @name predict.dfm
#' @aliases print.dfm_forecast
#' @aliases plot.dfm_forecast
#'
#' @title DFM Forecasts
#'
#' @description This function produces h-step ahead forecasts of both the factors and the data,
#' with an option to also forecast autocorrelated residuals with a univariate method and produce a combined forecast.
#'
#' @param object an object of class 'dfm'.
#' @param h integer. The forecast horizon.
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"2s"} or \code{"pca"}.
#' @param standardized logical. \code{FALSE} will return data forecasts on the original scale.
#' @param resFUN an (optional) function to compute a univariate forecast of the residuals.
#' The function needs to have a second argument providing the forecast horizon (\code{h}) and return a vector or forecasts. See Examples.
#' @param resAC numeric. Threshold for residual autocorrelation to apply \code{resFUN}: only residual series where AC1 > resAC will be forecasted.
#' @param \dots further arguments to \code{resFUN}.
#'
#' @returns A list-like object of class 'dfm_forecast' with the following elements:
#'  \item{\code{X_fcst}}{\eqn{h \times n}{h x n} matrix with the forecasts of the variables. }
#'  \item{\code{F_fcst}}{\eqn{h \times r}{h x r} matrix with the factor forecasts. }
#'  \item{\code{X}}{\eqn{T \times n}{T x n} matrix with the standardized (scaled and centered) data - with attributes attached allowing reconstruction of the original data:
#'       \tabular{llll}{
#'       \code{"stats"} \tab\tab is a \eqn{n \times 5}{n x 5} matrix of summary statistics of class \code{"qsu"} (see \code{\link[collapse]{qsu}}). Only attached if \code{standardized = TRUE}. \cr\cr
#'       \code{"attributes"} \tab\tab contains the \code{\link{attributes}} of the original data input.\cr\cr
#'       \code{"is.list"} \tab\tab is a logical value indicating whether the original data input was a list / data frame. \cr\cr
#'       }
#'      }
#'  \item{\code{F}}{\eqn{T \times r}{T x r} matrix of factor estimates. }
#'  \item{\code{method}}{the factor estimation method used.}
#'  \item{\code{anyNA}}{logical indicating whether \code{X} contains any missing values.}
#'  \item{\code{h}}{the forecast horizon.}
#'  \item{\code{resid.fc}}{logical indicating whether a univariate forecasting function was applied to the residuals.}
#'  \item{\code{resid.fc.ind}}{indices indicating for which variables (columns of \code{X}) the residuals were forecasted using the univariate function.}
#'  \item{\code{call}}{call object obtained from \code{match.call()}.}
#'
#' @examples \donttest{
#' library(xts)
#' library(collapse)
#'
#' # Fit DFM with 3 factors and 3 lags in the transition equation
#' mod = DFM(diff(BM14_M), r = 3, p = 3)
#'
#' # 15 period ahead forecast
#' fc = predict(mod, h = 15)
#' print(fc)
#' plot(fc, xlim = c(300, 370))
#'
#' # Also forecasting autocorrelated residuals with an AR(1)
#' fcfun <- function(x, h) predict(ar(na_rm(x)), n.ahead = h)$pred
#' fcar = predict(mod, resFUN = fcfun, h = 15)
#' plot(fcar, xlim = c(300, 370))
#'
#' # Retrieving a data frame of the forecasts
#' head(as.data.frame(fcar, pivot = "wide")) # Factors
#' head(as.data.frame(fcar, use = "data"))   # Data
#' head(as.data.frame(fcar, use = "both"))   # Both
#' }
#' @export
# TODO: Option for prediction in original format??
predict.dfm <- function(object,
                        h = 10L,
                        method = switch(object$em.method, none = "2s", "qml"),
                        standardized = TRUE,
                        resFUN = NULL,
                        resAC = 0.1, ...) {

  Fa <- switch(tolower(method),
              pca = object$F_pca, `2s` = object$F_2s, qml = object$F_qml,
              stop("Unkown method", method))
  nf <- dim(Fa)[2L]
  C <- object$C
  ny <- dim(C)[1L]
  A <- object$A
  r <- dim(A)[1L]
  p <- dim(A)[2L] / r
  X <- object$X_imp

  F_fc <- matrix(NA_real_, nrow = h, ncol = nf)
  X_fc <- matrix(NA_real_, nrow = h, ncol = ny)
  F_last <- ftail(Fa, p)   # dimnames(F_last) <- list(c("L2", "L1"), c("f1", "f2"))
  spi <- p:1

  # DFM forecasting loop
  for (i in seq_len(h)) {
    F_reg <- ftail(F_last, p)
    F_fc[i, ] <- tmp <- A %*% `dim<-`(t(F_reg)[, spi, drop = FALSE], NULL)
    dim(tmp) <- NULL
    X_fc[i, ] <- C %*% tmp
    F_last <- rbind(F_last, tmp)
  }

  # Additional univariate forecasting of the residuals?
  fcr <- NULL
  if(!is.null(resFUN)) {
    if(!is.function(resFUN)) stop("resFUN needs to be a forecasting function with second argument h that produces a numeric h-step ahead forecast of a univariate time series")
    # If X is a multivariate time series object for which the univariate forecasting function could have methods.
    ofl <- !attr(X, "is.list") && length(attr(X, "attributes")[["class"]])
    rsid <- residuals(object, method, orig.format = ofl, standardized = TRUE, na.keep = FALSE) # TODO: What about missing values??
    if(ofl && length(object$rm.rows)) rsid <- rsid[-object$rm.rows, , drop = FALSE]
    ACF <- AC1(rsid, object$anyNA)
    fcr <- which(abs(ACF) >= abs(resAC)) # TODO: Check length of forecast??
    for (i in fcr) X_fc[, i] <- X_fc[, i] + as.numeric(resFUN(rsid[, i], h, ...))
  } else if(!is.null(res <- object[["e"]])) {
    rho <- object$rho
    last_res <- res[nrow(res), ]
    for (i in seq_len(h)) {
      last_res <- last_res * rho
      X_fc[i, ] <- X_fc[i, ] + last_res
    }
  }

  if(!standardized) { # Unstandardize factors with the average mean and SD??
    stats <- attr(X, "stats")
    attr(X, "stats") <- NULL
    X_fc <- unscale(X_fc, stats)
    X <- unscale(X, stats)
  }

  dimnames(X_fc) <- dimnames(X)
  dimnames(F_fc) <- dimnames(Fa)

  if(object$anyNA) {
    X[attr(X, "missing")] <- NA
    attr(X, "missing") <- NULL
  }

  # model = object, # Better only save essential objects...
  res <- list(X_fcst = X_fc,
              F_fcst = F_fc,
              X = X,
              F = Fa,
              method = method,
              anyNA = object$anyNA,
              h = h,
              resid.fc = !is.null(resFUN), # TODO: Rename list elements??
              resid.fc.ind = fcr,
              call = match.call())

  class(res) <- "dfm_forecast"
  return(res)
}
# forecast.dfm <- predict.dfm

#' @rdname predict.dfm
#' @param x object of type 'dfm_forecast', returned from \code{predict.dfm}.
#' @param digits integer. The number of digits to print out.
#' @param \dots not used.
#' @export
print.dfm_forecast <- function(x,
                               digits = 4L, ...) {
  h <- x$h
  cat(h, "Step Ahead Forecast from Dynamic Factor Model\n\n")
  cat("Factor Forecasts\n")
  F_fcst <- x$F_fcst
  dimnames(F_fcst)[[1L]] <- seq_len(h)
  print(round(F_fcst, digits))
  cat("\nSeries Forecasts\n")
  X_fcst <- x$X_fcst
  dimnames(X_fcst)[[1L]] <- seq_len(h)
  print(round(X_fcst, digits))
}

#' @rdname predict.dfm
#' @param main,xlab,ylab character. Graphical parameters passed to \code{\link{ts.plot}}.
#' @param factors integers indicating which factors to display. Setting this to \code{NA}, \code{NULL} or \code{0} will omit factor plots.
#' @param scale.factors logical. Standardize factor estimates, this usually improves the plot since the factor estimates corresponding to the greatest PCA eigenvalues tend to have a greater variance than the data.
#' @param factor.col,factor.lwd graphical parameters affecting the colour and line width of factor estimates plots. See \code{\link{par}}.
#' @param fcst.lty integer or character giving the line type of the forecasts of factors and data. See \code{\link{par}}.
#' @param data.col character vector of length 2 indicating the colours of historical data and forecasts of that data. Setting this to \code{NA}, \code{NULL} or \code{""} will not plot data and data forecasts.
#' @param legend logical. \code{TRUE} draws a legend in the top-left of the chart.
#' @param legend.items character names of factors for the legend.
#' @param grid logical. \code{TRUE} draws a grid on the background of the plot.
#' @param vline logical. \code{TRUE} draws a vertical line deliminating historical data and forecasts.
#' @param vline.lty,vline.col graphical parameters affecting the appearance of the vertical line. See \code{\link{par}}.
#' @param \dots further arguments passed to \code{\link{ts.plot}}. Sensible choices are \code{xlim} and \code{ylim} to restrict the plot range.
#' @importFrom collapse setop
#' @export
# TODO: multiple plot types...# , type = c("joint", "individual")
# also arguments show = c("both", "factors", "data"), and
# Also put plot on original timescale if ts object
# TODO: Option to unstandardize factors.
plot.dfm_forecast <- function(x,
                              main = paste(x$h, "Period Ahead DFM Forecast"),
                              xlab = "Time", ylab = "Standardized Data",
                              factors = seq_len(ncol(x$F)),
                              scale.factors = TRUE,
                              factor.col = rainbow(length(factors)),
                              factor.lwd = 1.5,
                              fcst.lty = "dashed",
                              data.col = c("grey85", "grey65"),
                              legend = TRUE,
                              legend.items = paste0("f", factors),
                              grid = FALSE, vline = TRUE,
                              vline.lty = "dotted", vline.col = "black", ...) {

  dcl <- is.character(data.col[1L]) && nzchar(data.col[1L])
  ffl <- length(factors) && !is.na(factors[1L]) && factors[1L] > 0L
  nyliml <- !(...length() && any(names(list(...)) == "ylim")) # ...names() -> Added after R 3.3.0
  if(!ffl) factors <- 1L
  Fa <- x$F[, factors, drop = FALSE]
  r <- ncol(Fa)
  TT <- nrow(Fa)
  if(ffl) {
    F_fcst <- x$F_fcst[, factors, drop = FALSE]
    if(scale.factors) {
      fcstat <- qsu(Fa)
      F_fcst <- setop(TRA.matrix(F_fcst, fcstat[, "Mean"], "-"), "/", fcstat[, "SD"], rowwise = TRUE) # Unscale ??
      Fa <- fscale(Fa)
    }
    if(nyliml) Fr <- frange(Fa)
    Fa <- rbind(Fa, matrix(NA_real_, x$h, r))
  } else Fr <- NULL
  if(dcl) {
    X <- x$X
    n <- ncol(X)
    if(nyliml) {
      Xr <- frange(X, na.rm = TRUE)
      Pr <- frange(if(ffl) c(F_fcst, x$X_fcst) else x$X_fcst)
    }
    X_fcst <- rbind(matrix(NA_real_, TT-1L, n), X[TT, , drop = FALSE], x$X_fcst)
    X <- rbind(X, matrix(NA_real_, x$h, n))
  } else {
    data.col <- Xr <- NULL
    X <- Fa[, 1L]
    if(nyliml) Pr <- frange(F_fcst)
  }
  if(ffl) F_fcst <- rbind(matrix(NA_real_, TT-1L, r), Fa[TT, , drop = FALSE], F_fcst)
  if(nyliml) {
    ts.plot(X, col = data.col[1L],
            ylim = c(min(Xr[1L], Fr[1L], Pr[1L]), max(Xr[2L], Fr[2L], Pr[1L])),
            main = main, xlab = xlab, ylab = ylab, ...)
  } else ts.plot(X, col = data.col[1L], main = main, xlab = xlab, ylab = ylab, ...)
  if(grid) grid()
  if(dcl) for (i in seq_len(n)) lines(X_fcst[, i], col = data.col[2L], lty = fcst.lty)
  if(ffl) for (i in seq_len(r)) {
    lines(Fa[, i], col = factor.col[i], lwd = factor.lwd)
    lines(F_fcst[, i], col = factor.col[i], lwd = factor.lwd, lty = fcst.lty)
  }
  if(ffl && legend) legend("topleft", legend.items, col = factor.col,
                           lwd = factor.lwd, lty = 1L, bty = "n")
  if(vline) abline(v = TT, col = vline.col, lwd = 1L, lty = vline.lty)
}


#' @rdname predict.dfm
#' @param x an object class 'dfm_forecast'.
#' @param use character. Which forecasts to use \code{"factors"}, \code{"data"} or \code{"both"}.
#' @param pivot character. The orientation of the frame: \code{"long"} or \code{"wide"}.
#' @param time a vector identifying the time dimension, must be of length T + h, or \code{NULL} to omit a time variable.
#' @param stringsAsFactors logical. If \code{TRUE} and \code{pivot = "long"} the 'Variable' column is created as a factor. Same as option to \code{\link{as.data.frame.table}}.
#' @param \dots not used.
#'
#' @export
as.data.frame.dfm_forecast <- function(x, ...,
                              use = c("factors", "data", "both"),
                              pivot = c("long", "wide"),
                              time = seq_len(nrow(x$F) + x$h),
                              stringsAsFactors = TRUE) {

  mat <- switch(tolower(use[1L]),
                factors = rbind(x$F, x$F_fcst),
                data = rbind(x$X, x$X_fcst),
                both = cbind(rbind(x$F, x$F_fcst), rbind(x$X, x$X_fcst)),
                stop("Unknown use option:", use[1L]))

  fcvec <- c(rep(FALSE, nrow(x$F)), rep(TRUE, x$h))
  TT <- nrow(mat)
  r <- ncol(mat)
  if(!is.null(time) && length(time) != TT) stop(sprintf("time must be a length %s vector or NULL", TT))

  res <- switch(tolower(pivot[1L]),
      long = list(Variable = if(stringsAsFactors) setAttrib(rep(1:r, each = TT), list(levels = dimnames(mat)[[2L]], class = "factor")) else
                                                rep(dimnames(mat)[[2L]], each = TT),
                  Time = if(length(time)) rep(time, r) else NULL,
                  Forecast = rep(fcvec, r),
                  Value = unattrib(mat)),
      wide = c(list(Time = time, Forecast = fcvec), mctl(mat, TRUE)),
      stop("Unknown pivot option:", pivot[1L])
  )

  if(is.null(time)) res <- na_rm(res)
  attr(res, "row.names") <- .set_row_names(length(res[[1L]]))
  class(res) <- "data.frame"
  return(res)
}

# interpolate.dfm <- function(x, method = "qml", interpolate = TRUE) {
#   W <- is.na(data)
#   stats <- qsu(data)
#   STDdata <- fscale(data)
#   if(nrow(x$C) != ncol(data)) stop("dimension mismatch")
#   Fcst <- tcrossprod(x[[method]], x$C)
#   # TODO: Make this work for data.table...
#   STDdata[W] <- Fcst[W]
#   STDdata <- ((STDdata %r*% stats[, "SD"]) %r+% stats[, "Mean"])
#   data[W] <- STDdata[W]
#   data
# }
#
# nowcast.dfm <- function(x, method = "qml", ...) {
# }
#
# backcast.dfm <- function(x, method = "qml", ...) {
# }

# Adapted from: https://github.com/nmecsys/nowcasting/blob/master/R/ICfactors.R
#' @title Information Criteria to Determine the Number of Factors (r)
#' @description Minimizes 3 information criteria proposed by Bai and Ng (2002) to determine the optimal number of factors r* to be used in an approximate factor model.
#' A Screeplot can also be computed to eyeball the number of factors in the spirit of Onatski (2010).
#' @param X a \code{T x n} numeric data matrix or frame of stationary time series.
#' @param max.r integer. The maximum number of factors for which IC should be computed (or eigenvalues to be displayed in the screeplot).
#' @param x an object of type 'ICr'.
#' @param \dots further arguments to \code{\link{ts.plot}} or \code{\link{plot}}.
#'
#' @return A list of 4 elements:
#' \item{F_pca}{\code{T x n} matrix of principle component factor estimates.}
#' \item{eigenvalues}{the eigenvalues of the covariance matrix of \code{X}.}
#' \item{IC}{\code{r.max x 3} 'table' containing the 3 information criteria of Bai and Ng (2002), computed for all values of \code{r} from \code{1:r.max}.}
#' \item{r.star}{vector of length 3 containing the number of factors (\code{r}) minimizing each information criterion.}
#'
#' @details Following Bai and Ng (2002) and De Valk et al. (2019), let \eqn{NSSR(r)}{NSSR(r)} be the normalized sum of squared residuals \eqn{SSR(r) / (n \times T)}{SSR(r) / (n x T)} when r factors are estimated using principal components.
#' Then the information criteria can be written as follows:
#'
#' \deqn{IC_{r1} = \ln(NSSR(r)) + r\left(\frac{n + T}{nT}\right) + \ln\left(\frac{nT}{n + T}\right)}{ICr1 = ln(NSSR(r)) + r * (n + T)/(n * T) + ln((n * T)/(n + T))}
#' \deqn{IC_{r2} = \ln(NSSR(r)) + r\left(\frac{n + T}{nT}\right) + \ln(\min(n, T))}{ICr2 = ln(NSSR(r)) + r * (n + T)/(n * T) + ln(min(n, T))}
#' \deqn{IC_{r3} = \ln(NSSR(r)) + r\left(\frac{\ln(\min(n, T))}{\min(n, T)}\right)}{ICr3 = ln(NSSR(r)) + r * ln(min(n, T))/min(n, T)}
#'
#' The optimal number of factors r* corresponds to the minimum IC. The three criteria are are asymptotically equivalent, but may give significantly
#' different results for finite samples. The penalty in \eqn{IC_{r2}}{ICr2} is highest in finite samples.
#'
#' In the Screeplot a horizontal dashed line is shown signifying an eigenvalue of 1, or a share of variance corresponding to 1 divided by the number of eigenvalues.
#'
#' @note To determine the number of lags (\code{p}) in the factor transition equation, use the function \code{vars::VARselect} with r* principle components (also returned by \code{ICr}).
#'
#' @examples
#' library(xts)
#' library(vars)
#'
#' ics = ICr(diff(BM14_M))
#' print(ics)
#' plot(ics)
#' screeplot(ics)
#'
#' # Optimal lag-order with 6 factors chosen
#' VARselect(ics$F_pca[, 1:6])
#'
#' @references
#' Bai, J., Ng, S. (2002). Determining the Number of Factors in Approximate Factor Models. \emph{Econometrica, 70}(1), 191-221. \doi{10.1111/1468-0262.00273}
#'
#' Onatski, A. (2010). Determining the number of factors from empirical distribution of eigenvalues. \emph{The Review of Economics and Statistics, 92}(4), 1004-1016.
#'
#' De Valk, S., de Mattos, D., & Ferreira, P. (2019). Nowcasting: An R package for predicting economic variables using dynamic factor models. \emph{The R Journal, 11}(1), 230-244.
#' @export
ICr <- function(X, max.r = min(20, ncol(X)-1)) {

  # Converting to matrix and standardizing
  X <- fscale(qM(X), na.rm = TRUE)
  dimnames(X) <- NULL

  if(anyNA(X)) {
    message("Missing values detected: imputing data with tsnarmimp() with default settings")
    X <- tsnarmimp(X)
    attributes(X) <- list(dim = dim(X))
  }

  n <- ncol(X)
  TT <- nrow(X)

  # defining rmax and checking if it is a positive integer
  if(!is.numeric(max.r) || max.r < 1) stop("max.r needs to be a positive integer")
  max.r <- if(max.r > n) n else as.integer(max.r)

  # Eigen decomposition
  eigen_decomp <- eigen(cov(X), symmetric = TRUE)
  evs <- eigen_decomp$vectors
  F_pca <- X %*% evs

  # Various constant terms, according to the 3 criteria of Bai and Ng (2002)
  Tn <- TT * n
  npTdTn <- (n + TT) / Tn
  minnT <- min(n, TT)
  c1 <- npTdTn * log(1/npTdTn)
  c2 <- npTdTn * log(minnT)
  c3 <- log(minnT) / minnT
  cvec <- c(c1, c2, c3)
  result <- matrix(0, max.r, 3)

  # Calculating the IC
  for (r in 1:max.r) {
    # Residuals from r PC's
    res <- X - tcrossprod(F_pca[, 1:r], evs[, 1:r])
    # Log normalized sum of squared errors
    logV <- log(sum(colSums(res^2)/Tn))
    # Computing criteria
    result[r, ] <- logV + r * cvec
  }

  dimnames(result) <- list(r = 1:max.r, IC = paste0("IC", 1:3))
  class(result) <- "table"
  colnames(F_pca) <- paste0("PC", 1:n)

  res_obj <- list(F_pca = F_pca, eigenvalues = eigen_decomp$values, IC = result, r.star = apply(result, 2, which.min))
  class(res_obj) <- "ICr"
  return(res_obj)
}

#' @rdname ICr
#' @export
print.ICr <- function(x, ...) {
  cat("Optimal Number of Factors (r) from Bai and Ng (2002) Criteria\n\n")
  print(x$r.star)
}

#' @rdname ICr
#' @importFrom collapse fmin.matrix
#' @importFrom graphics grid points
#' @export
plot.ICr <- function(x, ...) {

  ts.plot(x$IC, gpars = list(xlab = "Number of Factors (r)", ylab = "IC Value", lty = c(2L, 1L, 3L),
                             main = "Optimal Number of Factors (r) from Bai and Ng (2002) Criteria"), ...)
  # grid()
  legend("topleft", paste0(names(x$r.star), ", r* = ", x$r.star), lty = c(2L, 1L, 3L)) # , bty = "n"
  points(x = x$r.star, y = fmin.matrix(x$IC), pch = 19, col ="red")

}

#' @rdname ICr
#' @param type character. Either \code{"ev"} (eigenvalues), \code{"pve"} (percent variance explained), or \code{"cum.pve"} (cumulative PVE). Multiple plots can be requested.
#' @param show.grid logical. \code{TRUE} shows gridlines in each plot.
#' @importFrom stats screeplot
#' @export
screeplot.ICr <- function(x, type = "pve", show.grid = TRUE, max.r = 30, ...) {
  ev <- x$eigenvalues
  n <- length(ev)
  pve <- (ev / sum(ev)) * 100
  cs_pve <- cumsum(pve)

  if(length(ev) > max.r) {
    ev <- ev[1:max.r]
    pve <- pve[1:max.r]
    cs_pve <- cs_pve[1:max.r]
  }

  ## This is smarter, but less flexible...
  # if(length(ev) > 20) {
  #   if(cs_pve[21] > 90) {
  #     pve = pve[1:20]
  #     cs_pve = cs_pve[1:20]
  #   } else {
  #     pve = pve[cs_pve < 90]
  #     cs_pve = cs_pve[cs_pve < 90]
  #   }
  # }

  if(length(type) > 1L) {
    oldpar <- par(mfrow = c(1, length(type)))
    on.exit(par(oldpar))
  }
  if(any(type == "ev")) {
    plot(ev, type = "o", ylab = "Eigenvalue", xlab = "Principal Component", col = "dodgerblue4", ...)
    if(show.grid) grid()
    abline(h = 1, lty = 2)
  }
  if(any(type == "pve")) {
    plot(pve, type = "o", ylab = "% Variance Explained", xlab = "Principal Component", col = "dodgerblue4", ...)
    if(show.grid) grid()
    abline(h = 100 / n, lty = 2)
  }
  if(any(type == "cum.pve")) {
    plot(cs_pve, type = "o", ylab = "Cumulative % Variance Explained", xlab = "Number of Principal Components", col = "brown3", ...)
    if(show.grid) grid()
  }
}

#' @rdname plot.dfm
#' @export
screeplot.dfm <- function(x, ...) {
  xl <- list(eigenvalues = x$eigen$values)
  screeplot.ICr(xl, ...)
}

