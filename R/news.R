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

dfm_news_state <- function(dfm, require_full = FALSE) {
  ss_full <- dfm[["ss_full"]]
  if(!is.null(ss_full)) {
    needed <- c("A", "C", "Q", "R", "F_0", "P_0")
    if(any(!needed %in% names(ss_full))) stop("ss_full is missing required state matrices")
    return(ss_full)
  }

  if(isTRUE(require_full)) {
    stop("news() requires dfm objects with full state matrices for MQ or idio.ar1 models; re-estimate with the updated package")
  }

  # Build full state from compact matrices (baseline case)
  A <- dfm$A
  C <- dfm$C
  Q <- dfm$Q
  R <- dfm$R
  r <- nrow(A)
  if(ncol(A) %% r != 0L) stop("Invalid transition matrix dimensions in dfm object")
  p <- ncol(A) / r
  rp <- r * p
  A_state <- matrix(0, rp, rp)
  A_state[1:r, 1:rp] <- A
  if(rp > r) A_state[(r + 1L):rp, 1:(rp - r)] <- diag(rp - r)
  C_state <- cbind(C, matrix(0, nrow(C), rp - r))
  Q_state <- matrix(0, rp, rp)
  Q_state[1:r, 1:r] <- Q
  F_0 <- rep(0, rp)
  P_0 <- dfm_news_stationary_cov(A_state, Q_state)
  list(A = A_state, C = C_state, Q = Q_state, R = R, F_0 = F_0, P_0 = P_0)
}

dfm_news_kfs <- function(X, state, k) {
  A <- state$A
  C <- state$C
  Q <- state$Q
  R <- state$R
  F_0 <- state$F_0
  P_0 <- state$P_0

  if(any(vapply(list(A, C, Q, R, F_0, P_0), is.null, logical(1L)))) {
    stop("state is missing required matrices")
  }

  state_dim <- ncol(C)
  if(nrow(A) != state_dim || ncol(A) != state_dim) stop("state transition matrix has incompatible dimensions")
  if(nrow(Q) != state_dim || ncol(Q) != state_dim) stop("state covariance matrix has incompatible dimensions")
  if(length(F_0) != state_dim) stop("state initial mean has incompatible dimensions")
  if(nrow(P_0) != state_dim || ncol(P_0) != state_dim) stop("state initial covariance has incompatible dimensions")
  if(nrow(C) != ncol(X)) stop("state observation matrix has incompatible dimensions")

  if(k > 0L) {
    aug_dim <- state_dim * (k + 1L)
    A_aug <- matrix(0, aug_dim, aug_dim)
    A_aug[1:state_dim, 1:state_dim] <- A
    A_aug[(state_dim + 1L):aug_dim, 1:(k * state_dim)] <- diag(k * state_dim)

    C_aug <- cbind(C, matrix(0, nrow(C), aug_dim - state_dim))
    Q_aug <- matrix(0, aug_dim, aug_dim)
    Q_aug[1:state_dim, 1:state_dim] <- Q

    F_0 <- c(F_0, rep(0, aug_dim - state_dim))
    P_0 <- rbind(cbind(P_0, matrix(0, state_dim, aug_dim - state_dim)),
                 matrix(0, aug_dim - state_dim, aug_dim))
    diag(P_0)[(state_dim + 1L):aug_dim] <- 1e-8
  } else {
    A_aug <- A
    C_aug <- C
    Q_aug <- Q
  }

  kfs_res <- SKFS(X, A_aug, C_aug, Q_aug, R, F_0, P_0, FALSE)

  F_sm <- kfs_res$F_smooth[, 1:state_dim, drop = FALSE]
  X_sm <- tcrossprod(F_sm, C)
  list(X_sm = X_sm, P = kfs_res$P_smooth, C = C, R = R, state_dim = state_dim)
}

dfm_news_restore_missing <- function(X) {
  W <- attr(X, "missing")
  if(!is.null(W)) X[W] <- NA
  qM(X)
}

dfm_news_stats <- function(X) {
  stats <- unclass(attr(X, "stats"))
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
  class(stats) <- NULL
  Mx <- stats[, "Mean"]
  Wx <- stats[, "SD"]
  TRA.matrix(TRA.matrix(X, Mx, "-"), Wx, "/")
}

resolve_vars <- function(target.vars, n, names) {
  if(is.null(target.vars)) return(seq_len(n))
  if(is.character(target.vars)) {
    if(is.null(names)) stop("target.vars is a name but data have no column names")
    return(as.integer(ckmatch(target.vars, names, e = "Unknown target.vars:")))
  }
  if(!is.numeric(target.vars)) stop("target.vars must be NULL, numeric indices, or names")
  target.vars <- as.integer(target.vars)
  if(any(target.vars < 1L | target.vars > n)) stop("target.vars is out of bounds")
  unique(target.vars)
}
