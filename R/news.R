dfm_news_build_companion <- function(A, r, p, k) {
  rp <- r * p
  rp_k <- r * (p + k)
  A_comp <- matrix(0, rp_k, rp_k)
  A_comp[1:r, 1:rp] <- A
  if(rp_k > r) A_comp[(r + 1L):rp_k, 1:(rp_k - r)] <- diag(rp_k - r)
  A_comp
}

dfm_news_stationary_cov <- function(A, Q) {
  rp <- nrow(A)
  M <- diag(rp^2) - kronecker(A, A)
  P_0 <- ainv(M) %*% unattrib(Q)
  dim(P_0) <- c(rp, rp)
  P_0
}

dfm_news_kfs <- function(X, A, C, Q, R, r, p, k) {
  A_comp <- dfm_news_build_companion(A, r, p, k)
  rp_k <- nrow(A_comp)
  C_aug <- cbind(C, matrix(0, nrow(C), rp_k - r))
  Q_comp <- matrix(0, rp_k, rp_k)
  Q_comp[1:r, 1:r] <- Q
  F_0 <- rep(0, rp_k)
  P_0 <- dfm_news_stationary_cov(A_comp, Q_comp)
  kfs_res <- SKFS(X, A_comp, C_aug, Q_comp, R, F_0, P_0, FALSE)
  X_sm <- tcrossprod(kfs_res$F_smooth[, 1:r, drop = FALSE], C)
  list(X_sm = X_sm, P = kfs_res$P_smooth)
}

dfm_news_restore_missing <- function(X) {
  W <- attr(X, "missing")
  if(!is.null(W)) X[W] <- NA
  qM(X)
}

dfm_news_stats <- function(X) {
  stats <- attr(X, "stats")
  if(is.null(stats)) {
    n <- ncol(X)
    return(list(Mx = rep(0, n), Wx = rep(1, n)))
  }
  list(Mx = stats[, "Mean"], Wx = stats[, "SD"])
}

dfm_news_unscale_vec <- function(x, Mx, Wx) {
  x * Wx + Mx
}

dfm_news_scale <- function(X, stats) {
  if(is.null(stats)) stop("stats are required to scale X")
  Mx <- stats[, "Mean"]
  Wx <- stats[, "SD"]
  TRA.matrix(TRA.matrix(X, Mx, "-"), Wx, "/")
}
