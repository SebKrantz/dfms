#' @name DFM_summary
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
#' @importFrom collapse qsu
#' @export
print.dfm <- function(x,
                      digits = 4L, ...) {

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

#' @rdname DFM_summary
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"twostep"} or \code{"pca"}.
#' @return Summary information following a dynamic factor model estimation.
#' @importFrom stats cov
#' @importFrom collapse pwcov
#' @export
summary.dfm <- function(object,
                        method = "qml", ...) {

  X <- object$X_imp
  F <- object[[method]]
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

#' @rdname DFM_summary
#' @param compact integer. Display a more compact printout: \code{0} prints everything, \code{1} omits the observation matrix [C] and covariance matrix [R], and \code{2} omits all disaggregated information - yielding a summary of only the factor estimates.
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
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"twostep"} or \code{"pca"}.
#' @param type character. The type of plot: \code{"joint"}, \code{"individual"} or \code{"residual"}.
#' @importFrom graphics boxplot
#' @export
plot.dfm <- function(x,
                     method = "qml",
                     type = c("joint", "individual", "residual"), ...) {
  F <- x[[method]]
  nf <- dim(F)[2L]
  switch(type[1L],
    joint = {
      Xr <- range(x$X_imp)
      Fr <- range(F)
      ts.plot(x$X_imp, col = "grey85", ylim = c(min(Xr[1L], Fr[1L]), max(Xr[2L], Fr[2L])))
      cols <- rainbow(nf)
      for (i in seq_len(nf)) lines(F[, i], col = cols[i])
      legend("topleft", paste("Factor", seq_len(nf)), col = cols, lty = 1, bty = "n")
    },
    individual = {
      d <- ceiling(nf / 2L)
      d <- if(d == 1) c(2L, 1L) else c(d, 1L)
      oldpar <- par(mfrow = d)
      on.exit(par(oldpar))
      for (i in seq_len(nf)) plot(F[, i], type = 'l', main = paste0("QML estimated factor ", i),
                                  ylab = "Value", xlab = "Time")
    },
    residual = boxplot(x$X_imp - tcrossprod(F, x$C), main = "Residuals by input variable"),
    stop("Unknown plot method: ", type[1L])
  )
}

#' @name residuals.dfm
#' @aliases fitted.dfm
#' @aliases resid.dfm
#'
#' @title DFM Residuals and Fitted Values
#'
#' @param object an object of class 'dfm'.
#' @param method character. The factor estimates to use: one of \code{"qml"}, \code{"twostep"} or \code{"pca"}.
#' @param orig.format logical. \code{TRUE} returns residuals/fitted values in a data format similar to \code{X}.
#' @param standardized logical. \code{FALSE} will put residuals/fitted values on the original data scale.
#' @importFrom collapse TRA.matrix mctl setAttrib pad
#' @export
residuals.dfm <- function(object,
                          method = "qml",
                          orig.format = FALSE,
                          standardized = FALSE, ...) {
  X <- object$X_imp
  X_pred <- tcrossprod(object[[method]], object$C)
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

resid.dfm <- residuals.dfm

#' @rdname residuals.dfm
#' @export
fitted.dfm <- function(object,
                       method = "qml",
                       orig.format = FALSE,
                       standardized = FALSE, ...) {
  X <- object$X_imp
  res <- tcrossprod(object[[method]], object$C)
  if(!standardized) res <- unscale(res, attr(X, "stats"))
  if(object$anyNA) res[attr(X, "missing")] <- NA
  if(orig.format) {
    if(length(object$na.rm)) res <- pad(res, object$na.rm, method = "vpos")
    if(attr(X, "is.list")) res <- mctl(res)
    return(setAttrib(res, attr(X, "attributes")))
  }
  return(qM(res))
}

#' @name predict.dfm
#' @aliases forecast.dfm
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
#' @param resFUN an (optional) function to compute a univariate forecast of the residuals.
#' The function needs to have a second argument providing the forecast horizon (\code{h}) and return a vector or forecasts. See Examples.
#' @param resAC numeric. Threshold for residual autocorrelation to apply \code{resFUN}: only residual series where AC1 > resAC will be forecasted.
#'
#' @examples
#' dfm <- DFM(diff(Eustockmarkets), 2, 2)
#' predict(dfm)
#' fcfun <- function(x, h) predict(ar(x), n.ahead = h)$pred
#' predict(dfm, resFUN = fcfun)
#'
#' @export

predict.dfm <- function(object,
                        h = 10L,
                        method = "qml",
                        standardized = TRUE,
                        resFUN = NULL,
                        resAC = 0.1, ...) {

  F <- object[[method]]
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

  for (i in seq_len(h)) {
    F_reg <- ftail(F_last, p)
    F_fc[i, ] <- A %*% `dim<-`(t(F_reg)[, spi, drop = FALSE], NULL)
    X_fc[i, ] <- C %*% F_fc[i, ]
    F_last <- rbind(F_last, F_fc[i, ])
  }
  # TODO: What about missing values??
  if(!is.null(resFUN)) {
    if(!is.function(resFUN)) stop("resFUN needs to be a forecasting function with second argument h that produces a numeric h-step ahead forecast of a univariate time series")
    ofl <- !attr(X, "is.list") && length(attr(X, "attributes")[["class"]])
    resid <- residuals(object, method, orig.format = ofl, standardized = TRUE)
    if(ofl && length(object$na.rm)) resid <- resid[-object$na.rm, , drop = FALSE] # drop = FALSE?
    ACF <- AC1(resid, object$anyNA)
    fcr <- which(abs(ACF) >= abs(resAC)) # TODO: Check length of forecast??
    for (i in fcr) X_fc[, i] <- X_fc[, i] + as.numeric(resFUN(resid[, i], h, ...))
  } else fcr <- NULL

  if(!standardized) {
    stats <- attr(X, "stats")
    X_fc <- unscale(X_fc, stats)
    X <- unscale(X, stats)
  }

  dimnames(X_fc)[[2L]] <- dimnames(X)[[2L]]

  if(object$anyNA) X[attr(X, "missing")] <- NA

  # model = object, # Better only save essential objects ??
  res <- list(X_fcst = X_fc,
              F_fcst = F_fc,
              X = X,
              F = F,
              method = method,
              h = h,
              resid.fc = !is.null(resFUN),
              resid.fc.ind = fcr,
              call = match.call())
  class(res) <- "dfm_forecast"
  return(res)
}

forecast.dfm <- predict.dfm

#' @rdname predict.dfm
#' @param digits integer. The number of digits to print out.
#' @export
print.dfm_forecast <- function(x,
                               digits = 4L, ...) {
  h <- x$h
  cat(h, "Step Ahead Forecast from Dynamic Factor Model\n\n")
  cat("Factor Forecasts\n")
  F_fcst <- x$F_fcst
  dimnames(F_fcst) <- list(seq_len(h), paste0("f", seq_len(ncol(F_fcst))))
  print(round(F_fcst, digits))
  cat("\nSeries Forecasts\n")
  X_fcst <- x$X_fcst
  dimnames(X_fcst)[[1L]] <- seq_len(h)
  print(round(X_fcst, digits))
}

#' @rdname predict.dfm
#' @export
# TODO: multiple plot types...
plot.dfm_forecast <- function(x, ...) { # , type = c("joint", "individual", "residual")
  F <- x$F
  Fr <- range(F)
  r <- ncol(F)
  cols <- rainbow(r)
  F <- rbind(F, matrix(NA_real_, x$h, r))
  X <- x$X_imp
  T <- nrow(X)
  n <- ncol(X)
  if(length(W <- attr(X, "missing"))) X[W] <- NA
  Xr <- range(X, na.rm = TRUE)
  Pr <- range(c(x$F_fcst, x$X_fcst))
  X <- rbind(X, matrix(NA_real_, x$h, n))
  F_fcst <- rbind(matrix(NA_real_, T, r), x$F_fcst)
  X_fcst <- rbind(matrix(NA_real_, T, n), x$X_fcst)
  ts.plot(X, col = "grey85", ylim = c(min(Xr[1L], Fr[1L], Pr[1L]), max(Xr[2L], Fr[2L], Pr[1L])))
  for (i in seq_len(n)) lines(X_fcst[, i], col = "grey50", lty = 3)
  for (i in seq_len(r)) {
    lines(F[, i], col = cols[i])
    lines(F_fcst[, i], col = cols[i], lty = 3)
  }
  legend("topleft", paste("Factor", seq_len(r)), col = cols, lty = 1, bty = "n")
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
