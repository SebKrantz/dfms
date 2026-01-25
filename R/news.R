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

dfm_news_stationary_cov_ar1 <- function(A, Q, rp_k, rho, R_idio) {
  state_dim <- nrow(A)
  n <- length(rho)

  # Factor block: solve Lyapunov for factor part
  A_fac <- A[1:rp_k, 1:rp_k, drop = FALSE]
  Q_fac <- Q[1:rp_k, 1:rp_k, drop = FALSE]
  P_0_fac <- dfm_news_stationary_cov(A_fac, Q_fac)

  # Error block: stationary AR(1) variance = sigma^2 / (1 - rho^2)
  P_0_err <- diag(diag(R_idio) / (1 - rho^2))

  # Combined (no cross-covariance between factors and errors at t=0)
  P_0 <- matrix(0, state_dim, state_dim)
  P_0[1:rp_k, 1:rp_k] <- P_0_fac
  P_0[(rp_k + 1L):state_dim, (rp_k + 1L):state_dim] <- P_0_err
  P_0
}

dfm_news_compute_X_sm <- function(F_sm, C, quarterly.vars = NULL, e_sm = NULL) {
  r <- ncol(C)
  if(is.null(quarterly.vars)) {
    X_sm <- tcrossprod(F_sm[, 1:r, drop = FALSE], C)
  } else {
    # MQ case: temporal aggregation for quarterly variables
    qind <- ckmatch(quarterly.vars, rownames(C))
    X_sm_m <- tcrossprod(F_sm[, 1:r, drop = FALSE], C[-qind, , drop = FALSE])
    Fa_lags <- flag(F_sm[, 1:r, drop = FALSE], 0:4, fill = 0, stubs = FALSE)
    Cq_lags <- C[qind, rep(1:r, each = 5), drop = FALSE] %r*% rep(c(1, 2, 3, 2, 1), r)
    X_sm_q <- tcrossprod(Fa_lags, Cq_lags)
    X_sm <- cbind(X_sm_m, X_sm_q)
  }
  if(!is.null(e_sm)) X_sm <- X_sm + e_sm
  X_sm
}

dfm_news_kfs <- function(X, A, C, Q, R, r, p, k,
                         idio.ar1 = FALSE, rho = NULL, R_idio = NULL,
                         quarterly.vars = NULL) {
  n <- nrow(C)

  if(idio.ar1 && !is.null(rho)) {
    # AR(1) case: augment state with error terms
    # State = [f_t, f_{t-1}, ..., f_{t-p+1}, f_{t-p}, ..., f_{t-p-k}, e_1, ..., e_n]
    rp <- r * p
    rp_k <- r * (p + k)
    state_dim <- rp_k + n

    # Build augmented A: factors + AR(1) errors
    A_aug <- matrix(0, state_dim, state_dim)
    A_aug[1:r, 1:rp] <- A
    if(rp_k > r) A_aug[(r + 1L):rp_k, 1:(rp_k - r)] <- diag(rp_k - r)
    A_aug[(rp_k + 1L):state_dim, (rp_k + 1L):state_dim] <- diag(rho)

    # Build augmented C: [C | zeros | I_n]
    C_aug <- cbind(C, matrix(0, n, rp_k - r), diag(n))

    # Build augmented Q
    Q_aug <- matrix(0, state_dim, state_dim)
    Q_aug[1:r, 1:r] <- Q
    Q_aug[(rp_k + 1L):state_dim, (rp_k + 1L):state_dim] <- R_idio

    # R for observation is kappa (small) since errors are in state
    R_obs <- 1e-4 * diag(n)

    # Initial conditions
    F_0 <- rep(0, state_dim)
    P_0 <- dfm_news_stationary_cov_ar1(A_aug, Q_aug, rp_k, rho, R_idio)

    kfs_res <- SKFS(X, A_aug, C_aug, Q_aug, R_obs, F_0, P_0, FALSE)

    F_sm <- kfs_res$F_smooth[, 1:r, drop = FALSE]
    e_sm <- kfs_res$F_smooth[, (rp_k + 1L):state_dim, drop = FALSE]
    # Only return factor covariances for P1/P2 computation
    P <- kfs_res$P_smooth[1:rp_k, 1:rp_k, , drop = FALSE]

  } else {
    # Standard case (no AR1)
    A_comp <- dfm_news_build_companion(A, r, p, k)
    rp_k <- nrow(A_comp)
    C_aug <- cbind(C, matrix(0, n, rp_k - r))
    Q_comp <- matrix(0, rp_k, rp_k)
    Q_comp[1:r, 1:r] <- Q
    F_0 <- rep(0, rp_k)
    P_0 <- dfm_news_stationary_cov(A_comp, Q_comp)

    kfs_res <- SKFS(X, A_comp, C_aug, Q_comp, R, F_0, P_0, FALSE)

    F_sm <- kfs_res$F_smooth[, 1:r, drop = FALSE]
    e_sm <- NULL
    P <- kfs_res$P_smooth
  }

  X_sm <- dfm_news_compute_X_sm(F_sm, C, quarterly.vars, e_sm)
  list(X_sm = X_sm, P = P)
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
