#' @name summary.dfm
#' @aliases print.dfm
#' @aliases summary.dfm
#' @aliases print.dfm_summary
#'
#' @title DFM Summary Methods
#'
#' @description Summary and print methods for class 'dfm'. \code{print.dfm} just prints basic model information and the factor transition matrix [A],
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
  fnam <- paste0("f", seq_len(r))
  cat("\nFactor Transition Matrix [A]\n")
  print(round(A, digits))
}

#' @rdname summary.dfm
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"twostep"} or \code{"pca"}.
#' @param \dots not used.
#' @return Summary information following a dynamic factor model estimation.
#' @importFrom stats cov
#' @importFrom collapse pwcov
#' @export
summary.dfm <- function(object, method = switch(object$em.method, none = "twostep", "qml"), ...) {

  X <- object$X_imp
  F <- switch(method, pca = object$F_pca, twostep = object$F_twostep, qml = object$F_qml, stop("Unkown method", method))
  A <- object$A
  r <- dim(A)[1L]
  p <- dim(A)[2L] / r
  C <- object$C
  res <- X - tcrossprod(F, C)
  anymissing <- object$anyNA
  if(anymissing) res[attr(X, "missing")] <- NA
  rescov <- pwcov(res, use = if(anymissing) "pairwise.complete.obs" else "everything", P = TRUE)
  ACF <- AC1(res, anymissing)
  R2 <- 1 - diag(rescov[,, 1L])
  summ <- list(info = c(n = dim(X)[2L], T = dim(X)[1L], r = r, p = p,
                        `%NA` = if(anymissing) sum(attr(X, "missing")) / prod(dim(X)) * 100 else 0),
               call = object$call,
               F_stats = msum(F),
               A = A,
               F_cov = pwcov(F, P = TRUE),
               Q = object$Q,
               C = C,
               R_diag = diag(object$R),
               res_cov = rescov,
               res_ACF = ACF,
               R2 = R2,
               R2_stats = msum(R2))
  class(summ) <- "dfm_summary"
  return(summ)
}

#' @rdname summary.dfm
#' @param compact integer. Display a more compact printout: \code{0} prints everything, \code{1} omits the observation matrix [C] and covariance matrix [R], and \code{2} omits all disaggregated information - yielding a summary of only the factor estimates.
#' @param \dots not used.
#' @export
print.dfm_summary <- function(x,
                              digits = 4L,
                              compact = sum(x$info["n"] > 15, x$info["n"] > 40), ...) {

  inf <- as.integer(x$info[1:4])
  cat("Dynamic Factor Model: n = ", inf[1L], ", T = ", inf[2L], ", r = ", inf[3L], ", p = ", inf[4L],
      ", %NA = ", round(x$info[5L], digits), "\n", sep = "")
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
  cat("\nGoodness of Fit: R-Squared\n")
  print(x$R2, digits)
  }
  cat("\nSummary of Individual R-Squared's\n")
  print(x$R2_stats, digits)
}


#' Plot DFM
#' @param x an object class 'dfm'.
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"twostep"}, \code{"pca"} or \code{"all"} to plot all estimates.
#' @param type character. The type of plot: \code{"joint"}, \code{"individual"} or \code{"residual"}.
#' @param scale.factors logical. Standardize factor estimates, this usually improves the plot since the factor estimates corresponding to the greatest PCA eigenvalues tend to have a greater variance than the data.
#' @param \dots for \code{plot.dfm}: further arguments to \code{\link{plot}}, \code{\link{ts.plot}}, or \code{\link{boxplot}}, depending on the \code{type} of plot. For \code{screeplot.dfm}: further arguments to \code{\link{screeplot.ICr}}.
#' @importFrom graphics boxplot
#' @importFrom collapse unlist2d ckmatch na_rm seq_row
#' @export
plot.dfm <- function(x,
                     method = switch(x$em.method, none = "twostep", "qml"),
                     type = c("joint", "individual", "residual"),
                     scale.factors = TRUE, ...) {
  F <- switch(method[1L],
              all = cbind(x$F_pca, setColnames(x$F_twostep, paste("TwoStep", colnames(x$F_twostep))),
                          if(length(x$F_qml)) setColnames(x$F_qml, paste("QML", colnames(x$F_qml))) else NULL),
              pca = x$F_pca, twostep = x$F_twostep, qml = x$F_qml, stop("Unknown method:", method[1L]))
  nf <- dim(F)[2L]
  switch(type[1L],
    joint = {
      Xr <- frange(x$X_imp)
      if(scale.factors) F <- fscale(F)
      Fr <- frange(F)
      ts.plot(x$X_imp, col = "grey85", ylim = c(min(Xr[1L], Fr[1L]), max(Xr[2L], Fr[2L])),
              ylab = "Value", main = "Standardized Series and Factor Estimates", ...)
      cols <- rainbow(nf)
      for (i in seq_len(nf)) lines(F[, i], col = cols[i])
      legend("topleft", colnames(F), col = cols, lty = 1, bty = "n", ncol = if(method[1L] == "all") 3L else 1L)
    },
    individual = { # TODO: Reduce plot margins
      if(method[1L] == "all") {
        qml <- !is.null(x$F_qml)
        nf <- nf / (2L + qml)
        oldpar <- par(mfrow = c(nf, 1L))
        on.exit(par(oldpar))
        for (i in seq_len(nf)) {
          plot(F[, i], type = 'l', main = paste("Factor", i), col = "red", ylab = "Value",
               xlab = if(i == nf) "Time" else "", ...)
          lines(F[, i + nf], type = 'l', col = "orange")
          if(qml) lines(F[, i + 2L * nf], type = 'l', col = "blue")
          if(i == 1L) legend("topleft", c("PCA", "TwoStep", if(qml) "QML"),
                             col = c("red", "orange", "blue"), lty = 1, bty = "n")
        }
      } else {
        oldpar <- par(mfrow = c(nf, 1L))
        on.exit(par(oldpar))
        cnF <- colnames(F)
        for (i in seq_len(nf)) plot(F[, i], type = 'l', main = cnF[i], ylab = "Value",
                                    xlab = if(i == nf) "Time" else "" , ...)
      }
    },
    residual = {
      if(method[1L] == "all") stop("Need to choose a specific method for residual plots")
      oldpar <- par(mar = c(11.5, 4.1, 4.1, 2.1))
      on.exit(par(oldpar))
      boxplot(x$X_imp - tcrossprod(F, x$C), main = "Residuals by input variable", las = 2, ...)
    },
    stop("Unknown plot type: ", type[1L])
  )
}


#' Extract Factor Estimates in a Data Frame
#' @param x an object class 'dfm'.
#' @param method character. The factor estimates to use: any of \code{"qml"}, \code{"twostep"}, \code{"pca"} (multiple can be supplied) or \code{"all"} for all estimates.
#' @param pivot character. The orientation of the frame: \code{"long"}, \code{"wide.factor"} or \code{"wide.method"} or \code{"wide"}.
#' @param time a vector identifying the time dimension, or \code{NULL} to omit a time variable.
#' @param stringsAsFactors make factors from method and factor identifiers. Same as option to \code{\link{as.data.frame.table}}.
#' @param \dots not used.
#'
#' @importFrom collapse ckmatch na_rm seq_row t_list unattrib
#' @importFrom stats setNames
#' @export
as.data.frame.dfm <- function(x, ...,
                              method = "all",
                              pivot = c("long", "wide.factor", "wide.method", "wide", "t.wide"),
                              time = seq_row(x$F_pca),
                              stringsAsFactors = TRUE) {

  estm <- c(PCA = "pca", TwoStep = "twostep", QML = "qml")
  if(length(method) > 1L || method != "all")
     estm <- estm[ckmatch(method, estm, e = "Unknown method:")]
  estlist <- x[paste0("F_", estm)]
  names(estlist) <- names(estm)
  estlist <- na_rm(estlist) # Also removes NULL elements

  nam <- names(estlist)
  m <- length(estlist)
  T <- nrow(estlist[[1L]])
  r <- ncol(estlist[[1L]])

  if(!is.null(time) && length(time) != T) {
    if(length(x$na.rm)) time <- time[-x$na.rm]
    if(length(time) != T) stop(sprintf("time must be a length %s vector or NULL", T))
  }

  res <- switch(pivot[1L],
    long = list(Method = if(stringsAsFactors) setAttrib(rep(1:m, each = T*r), list(levels = nam, class = "factor")) else rep(nam, each = T*r),
                Factor = if(stringsAsFactors) setAttrib(rep(1:r, times = m, each = T), list(levels = paste0("f", 1:r), class = "factor")) else rep(paste0("f", 1:r), times = m, each = T),
                Time = if(length(time)) rep(time, times = m*r) else NULL,
                Value = unlist(estlist, use.names = FALSE)),
    wide.factor = c(list(Method = if(stringsAsFactors) setAttrib(rep(1:m, each = T), list(levels = nam, class = "factor")) else rep(nam, each = T),
                         Time = if(length(time)) rep(time, times = m) else NULL),
                    setNames(lapply(t_list(unattrib(lapply(estlist, mctl))), unlist, FALSE, FALSE), paste0("f", 1:r))),
    wide.method = c(list(Factor = if(stringsAsFactors) setAttrib(rep(1:r, each = T), list(levels = paste0("f", 1:r), class = "factor")) else rep(paste0("f", 1:r), each = T),
                         Time = if(length(time)) rep(time, times = r) else NULL),
                    lapply(estlist, unattrib)),
    # If only one method, do not do combine names e.g. "QML_f1"? -> most of the time people just want a simple frame like this...
    wide = c(list(Time = time), setNames(unlist(lapply(estlist, mctl), FALSE, FALSE), if(length(nam) == 1L) paste0("f", 1:r) else t(outer(nam, 1:r, paste, sep = "_f")))),
    t.wide = c(list(Time = time), setNames(unlist(t_list(lapply(estlist, mctl)), FALSE, FALSE), if(length(nam) == 1L) paste0("f", 1:r) else outer(nam, 1:r, paste, sep = "_f"))),
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
#'
#' @param object an object of class 'dfm'.
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"twostep"} or \code{"pca"}.
#' @param orig.format logical. \code{TRUE} returns residuals/fitted values in a data format similar to \code{X}.
#' @param standardized logical. \code{FALSE} will put residuals/fitted values on the original data scale.
#' @param \dots not used.
#' @importFrom collapse TRA.matrix mctl setAttrib pad
#' @export
residuals.dfm <- function(object,
                          method = switch(object$em.method, none = "twostep", "qml"),
                          orig.format = FALSE,
                          standardized = FALSE, ...) {
  X <- object$X_imp
  F <- switch(method, pca = object$F_pca, twostep = object$F_twostep, qml = object$F_qml,
              stop("Unkown method", method))
  X_pred <- tcrossprod(F, object$C)
  if(!standardized) {
    stats <- attr(X, "stats")
    X_pred <- unscale(X_pred, stats)
    res <- unscale(X, stats) - X_pred
  } else res <- X - X_pred
  if(object$anyNA) res[attr(X, "missing")] <- NA
  if(orig.format) {
    if(length(object$na.rm)) res <- pad(res, object$na.rm, method = "vpos")
    if(attr(X, "is.list")) res <- mctl(res)
    return(setAttrib(res, attr(X, "attributes")))
  }
  return(qM(res))
}

#' @rdname residuals.dfm
#' @export
fitted.dfm <- function(object,
                       method = switch(object$em.method, none = "twostep", "qml"),
                       orig.format = FALSE,
                       standardized = FALSE, ...) {
  X <- object$X_imp
  F <- switch(method, pca = object$F_pca, twostep = object$F_twostep, qml = object$F_qml,
              stop("Unkown method", method))
  res <- tcrossprod(F, object$C)
  if(!standardized) res <- unscale(res, attr(X, "stats"))
  if(object$anyNA) res[attr(X, "missing")] <- NA
  if(orig.format) {
    if(length(object$na.rm)) res <- pad(res, object$na.rm, method = "vpos")
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
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"twostep"} or \code{"pca"}.
#' @param standardized logical. \code{FALSE} will return data forecasts on the original scale.
#' @param resFUN an (optional) function to compute a univariate forecast of the residuals.
#' The function needs to have a second argument providing the forecast horizon (\code{h}) and return a vector or forecasts. See Examples.
#' @param resAC numeric. Threshold for residual autocorrelation to apply \code{resFUN}: only residual series where AC1 > resAC will be forecasted.
#' @param \dots further arguments to \code{resFUN}.
#' @examples
#' dfm <- DFM(diff(Seatbelts[, 1:7], lag = 12), 3, 3)
#' predict(dfm)
#' fcfun <- function(x, h) predict(ar(x), n.ahead = h)$pred
#' predict(dfm, resFUN = fcfun)
#'
#' @export
# TODO: Option for prediction in original format??
predict.dfm <- function(object,
                        h = 10L,
                        method = switch(object$em.method, none = "twostep", "qml"),
                        standardized = TRUE,
                        resFUN = NULL,
                        resAC = 0.1, ...) {

  F <- switch(method, pca = object$F_pca, twostep = object$F_twostep, qml = object$F_qml,
              stop("Unkown method", method))
  nf <- dim(F)[2L]
  C <- object$C
  ny <- dim(C)[1L]
  A <- object$A
  r <- dim(A)[1L]
  p <- dim(A)[2L] / r
  X <- object$X_imp

  F_fc <- matrix(NA_real_, nrow = h, ncol = nf)
  X_fc <- matrix(NA_real_, nrow = h, ncol = ny)
  F_last <- ftail(F, p)   # dimnames(F_last) <- list(c("L2", "L1"), c("f1", "f2"))
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
  if(!is.null(resFUN)) { # TODO: What about missing values??
    if(!is.function(resFUN)) stop("resFUN needs to be a forecasting function with second argument h that produces a numeric h-step ahead forecast of a univariate time series")
    # If X is a multivariate time series object for which the univariate forecasting function could have methods.
    ofl <- !attr(X, "is.list") && length(attr(X, "attributes")[["class"]])
    rsid <- residuals(object, method, orig.format = ofl, standardized = TRUE)
    if(ofl && length(object$na.rm)) rsid <- rsid[-object$na.rm, , drop = FALSE] # drop = FALSE?
    ACF <- AC1(rsid, object$anyNA)
    fcr <- which(abs(ACF) >= abs(resAC)) # TODO: Check length of forecast??
    for (i in fcr) X_fc[, i] <- X_fc[, i] + as.numeric(resFUN(rsid[, i], h, ...))
  } else fcr <- NULL

  if(!standardized) { # Unstandardize factors with the average mean and SD??
    stats <- attr(X, "stats")
    X_fc <- unscale(X_fc, stats)
    X <- unscale(X, stats)
  }

  dimnames(X_fc) <- dimnames(X)
  dimnames(F_fc) <- dimnames(F)

  if(object$anyNA) X[attr(X, "missing")] <- NA

  # model = object, # Better only save essential objects...
  res <- list(X_fcst = X_fc,
              F_fcst = F_fc,
              X = X,
              F = F,
              method = method,
              h = h,
              resid.fc = !is.null(resFUN), # TODO: Rename list elements??
              resid.fc.ind = fcr,
              call = match.call())

  class(res) <- "dfm_forecast"
  return(res)
}
# forecast.dfm <- predict.dfm

# TODO: data frame method.

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
                              factors = 1:ncol(x$F),
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
  F <- x$F[, factors, drop = FALSE]
  r <- ncol(F)
  T <- nrow(F)
  if(ffl) {
    F_fcst <- x$F_fcst[, factors, drop = FALSE]
    if(scale.factors) {
      fcstat <- qsu(F)
      F_fcst <- setop(TRA.matrix(F_fcst, fcstat[, "Mean"], "-"), "/", fcstat[, "SD"], rowwise = TRUE)
      F <- fscale(F)
    }
    if(nyliml) Fr <- frange(F)
    F <- rbind(F, matrix(NA_real_, x$h, r))
  } else Fr <- NULL
  if(dcl) {
    X <- x$X
    n <- ncol(X)
    if(nyliml) {
      Xr <- frange(X, na.rm = TRUE)
      Pr <- frange(if(ffl) c(F_fcst, x$X_fcst) else x$X_fcst)
    }
    X_fcst <- rbind(matrix(NA_real_, T-1L, n), X[T, , drop = FALSE], x$X_fcst)
    X <- rbind(X, matrix(NA_real_, x$h, n))
  } else {
    data.col <- Xr <- NULL
    X <- F[, 1L]
    if(nyliml) Pr <- frange(F_fcst)
  }
  if(ffl) F_fcst <- rbind(matrix(NA_real_, T-1L, r), F[T, , drop = FALSE], F_fcst)
  if(nyliml) {
    ts.plot(X, col = data.col[1L],
            ylim = c(min(Xr[1L], Fr[1L], Pr[1L]), max(Xr[2L], Fr[2L], Pr[1L])),
            main = main, xlab = xlab, ylab = ylab, ...)
  } else ts.plot(X, col = data.col[1L], main = main, xlab = xlab, ylab = ylab, ...)
  if(grid) grid()
  if(dcl) for (i in seq_len(n)) lines(X_fcst[, i], col = data.col[2L], lty = fcst.lty)
  if(ffl) for (i in seq_len(r)) {
    lines(F[, i], col = factor.col[i], lwd = factor.lwd)
    lines(F_fcst[, i], col = factor.col[i], lwd = factor.lwd, lty = fcst.lty)
  }
  if(ffl && legend) legend("topleft", legend.items, col = factor.col,
                           lwd = factor.lwd, lty = 1L, bty = "n")
  if(vline) abline(v = T, col = vline.col, lwd = 1L, lty = vline.lty)
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
#' @param X a \code{T x n} data matrix or frame.
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
#' @details Following Bai and Ng (2002) and De Valk et al. (2019), let NSSR(r) be the normalized sum of squared residuals [= SSR(r) / (n x T)] when r factors are estimated using principal components.
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
#' @references
#' Bai, J., Ng, S. (2002). Determining the Number of Factors in Approximate Factor Models. \emph{Econometrica, 70}(1), 191-221. <doi:10.1111/1468-0262.00273>
#'
#' Onatski, A. (2010). Determining the number of factors from empirical distribution of eigenvalues. \emph{The Review of Economics and Statistics, 92}(4), 1004-1016.
#'
#' De Valk, S., de Mattos, D., & Ferreira, P. (2019). Nowcasting: An R package for predicting economic variables using dynamic factor models. \emph{The R Journal, 11}(1), 230-244.
#' @export
ICr <- function(X, max.r = min(20, ncol(X)-1)) {

  # Converting to matrix and standardizing
  X <- fscale(qM(X))

  if(anyNA(X)) {
    message("Missing values detected: imputing data with tsremimpNA() with default settings")
    X <- tsremimpNA(X)$X_imp
  }

  n <- ncol(X)
  T <- nrow(X)

  # defining rmax and checking if it is a positive integer
  if(max.r < 1L || max.r != as.integer(max.r)) stop("rmax needs to be a positive integer")
  else if(max.r > n) max.r <- n

  # Eigen decomposition
  eigen_decomp = eigen(cov(X), symmetric = TRUE)
  evs = eigen_decomp$vectors
  F_pca = X %*% evs

  # Various constant terms, according to the 3 criteria of Bai and Ng (2002)
  Tn <- T * n
  npTdTn <- (n + T) / Tn
  minnT <- min(n, T)
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
screeplot.ICr <- function(x, type = c("pve", "cum.pve"), show.grid = TRUE, max.r = 30, ...) {
  ev = x$eigenvalues
  n = length(ev)
  pve = (ev / sum(ev)) * 100
  cs_pve = cumsum(pve)

  if(length(ev) > max.r) {
    ev = ev[1:max.r]
    pve = pve[1:max.r]
    cs_pve = cs_pve[1:max.r]
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

