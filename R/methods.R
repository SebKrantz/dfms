
#' Print DFM
#' @importFrom collapse qsu
#' @export
print.dfm <- function(x, digits = 4, ...) {

  X <- x$X_imp
  A <- x$A
  r <- dim(A)[1L]
  p <- dim(A)[2L]/r
  cat("Dynamic Factor Model: n = ", dim(X)[2L], ", T = ", dim(X)[1L], ", r = ", r, ", p = ", p, ", %NA = ",
      if(x$anyNA) round(sum(attr(X, "missing"))/prod(dim(X))*100, digits) else 0,"\n", sep = "")
  fnam <- paste0("f", seq_len(r))
  dimnames(A) <- list(fnam, as.vector(t(outer(paste0("L", seq_len(p)), fnam, paste, sep = "."))))
  cat("\nFactor Transition Matrix [A]\n")
  print(round(A, digits))
}

msum <- function(x) {
   stats <- qsu(x)
   med <- fmedian(x)
   res <- if(is.matrix(x)) cbind(stats[, 1:2], Median = med, stats[, -(1:2)]) else
                           c(stats[1:2], Median = med, stats[-(1:2)])
   class(res) <- "qsu"
   res
}

#' Summary information on dynamic factor model estimation
#'
#' @param x An object of \code{dfm} class
#' @return Prints out a summary information following a dynamic factor
#' model estimation. Also can return summary plots.
#' @importFrom stats cov
#' @importFrom collapse pwcov
#' @export
summary.dfm <- function(object, method = "qml", ...) {
  # TODO: Compact summary option: Don't print C and R
  X <- object$X_imp
  nam <- dimnames(X)[[2L]]
  F <- object[[method]]
  A <- object$A
  r <- dim(A)[1L]
  p <- dim(A)[2L] / r
  fnam <- paste0("f", seq_len(r))
  dimnames(F)[[2L]] <- fnam
  dimnames(A) <- list(fnam, as.vector(t(outer(paste0("L", seq_len(p)), fnam, paste, sep = "."))))
  unam <- paste0("u", seq_len(r))
  Q <- object$Q
  dimnames(Q) <- list(unam, unam)
  C <- object$C
  dimnames(C) <- list(nam, fnam)
  res <- X - tcrossprod(F, C)
  anymissing <- object$anyNA
  if(anymissing) res[attr(X, "missing")] <- NA
  rescov <- pwcov(res, use = if(anymissing) "pairwise.complete.obs" else "everything", P = TRUE)
  ACF <- cov(res[-nrow(res), ], res[-1L, ], use = if(anymissing) "pairwise.complete.obs" else "everything")
  ACF <- diag(ACF)/fvar(res)
  R2 <- 1 - diag(rescov[,, 1L])
  summ <- list(info = c(n = dim(X)[2L], T = dim(X)[1L], r = r, p = p,
                        `%NA` = if(anymissing) sum(attr(X, "missing")) / prod(dim(X)) * 100 else 0),
               call = object$call,
               F_stats = msum(F),
               A = A,
               F_cov = pwcov(F, P = TRUE),
               Q = Q,
               C = C,
               R_diag = setNames(diag(object$R), nam),
               res_cov = rescov,
               res_ACF = ACF,
               R2 = R2,
               R2_stats = msum(R2))
  class(summ) <- "dfm_summary"
  return(summ)
}

#' @export
print.dfm_summary <- function(x, digits = 4L, compact = sum(x$info["n"] > 15, x$info["n"] > 40), ...) {

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
#' @param x object of class 'dfm'.
#' @param method character. \code{"qml"}, \code{"twostep"} or \code{"pca"}.
#' @param type character. The type of plot: \code{"joint"}, \code{"individual"} or \code{"residual"}.
#' @importFrom graphics boxplot
#' @export
plot.dfm <- function(x, method = "qml", type = c("joint", "individual", "residual"), ...) {
  F <- x[[method]]
  nf <- dim(F)[2L]
  switch(type[1L],
    joint = {
      ts.plot(x$X_imp, col = "grey85")
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

unscale <- function(x, stats) TRA.matrix(TRA.matrix(x, stats[, "SD"], "*"), stats[, "Mean"], "+")

#' Residuals from DFM
#' @param object object of class 'dfm'.
#' @param method character. \code{"qml"}, \code{"twostep"} or \code{"pca"}.
#' @param orig.format logical. \code{TRUE} returns residuals in a data format similar to \code{X}.
#' @param standardized logical. \code{FALSE} will put residuals on the original data scale.
#' @importFrom collapse TRA.matrix mctl setAttrib pad
#' @export
residuals.dfm <- function(object, method = "qml", orig.format = FALSE, standardized = FALSE, ...) {
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

#' Fitted Values from DFM
#' @inheritParams residuals.dfm
#' @export
fitted.dfm <- function(object, method = "qml", orig.format = FALSE, standardized = FALSE, ...) {
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

ftail <- function(x, p) {n <- dim(x)[1L]; x[(n-p+1L):n, , drop = FALSE]}

#' DFM Forecasts
#' @param object object of class 'dfm'.
#' @param h forecast horizon.
#' @param method character. \code{"qml"}, \code{"twostep"} or \code{"pca"}.
#' @param resid.fc.FUN function to forecast residuals.
#' @export
# TODO: univariate forecast and non-standardized forecast...
predict.dfm <- function(object, h = 10L, method = "qml", resid.fc.FUN = NULL, ...) {

  F <- object[[method]]
  nf <- dim(F)[2L]
  C <- object$C
  ny <- dim(C)[1L]
  A <- object$A
  r <- dim(A)[1L]
  p <- dim(A)[2L] / r

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

  dimnames(X_fc)[[2L]] <- dimnames(object$X_imp)[[2L]]

  # model = object, # Better only save essential objects ??
  res <- list(X_fcst = X_fc,
              F_fcst = F_fc,
              X_imp = object$X_imp,
              F = F,
              method = method,
              h = h,
              resid.fc = !is.null(resid.fc.FUN),
              call = match.call())
  class(res) <- "dfm_fc"
  return(res)
}

forecast.dfm <- predict.dfm

#' @export
print.dfm_fc <- function(x, digits = 4L, ...) {
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

#' @export
# TODO: Proper plot margins and multiple plot types...
plot.dfm_fc <- function(x, ...) { # , type = c("joint", "individual", "residual")
  F <- x$F
  r <- ncol(F)
  cols <- rainbow(r)
  F <- rbind(F, matrix(NA_real_, x$h, r))
  X <- x$X_imp
  T <- nrow(X)
  n <- ncol(X)
  X <- rbind(X, matrix(NA_real_, x$h, n))
  if(length(W <- attr(X, "missing"))) X[W] <- NA
  F_fcst <- rbind(matrix(NA_real_, T, r), x$F_fcst)
  X_fcst <- rbind(matrix(NA_real_, T, n), x$X_fcst)
  ts.plot(X, col = "grey85")
  for (i in seq_len(n)) lines(X_fcst[, i], col = "grey40", lty = 3)
  for (i in seq_len(r)) lines(F[, i], col = cols[i])
  for (i in seq_len(r)) lines(F_fcst[, i], col = cols[i], lty = 3)
  legend("topleft", paste("Factor", seq_len(nf)), col = cols, lty = 1, bty = "n")
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
