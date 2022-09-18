
library(collapse)

# Two-Step factor estimates from monthly BM (2014) data
X <- fscale(diff(qM(BM14_M)))
r <- 5L # 5 Factors
p <- 3L # 3 Lags
n <- ncol(X)

# Initializing the Kalman Filter with PCA results
X_imp <- tsnarmimp(X)                 # Imputing Data
v <- eigen(cov(X_imp))$vectors[, 1:r] # PCA
F_pc <- X_imp %*% v                   # Principal component factor estimates
C <- cbind(v, matrix(0, n, r*p-r))    # Observation matrix
res <- X - tcrossprod(F_pc, v)        # Residuals from static predictions
R <- diag(fvar(res))                  # Observation residual covariance
var <- .VAR(F_pc, p)                  # VAR(p)
A <- rbind(t(var$A), diag(1, r*p-r, r*p))
Q <- matrix(0, r*p, r*p)              # VAR residual matrix
Q[1:r, 1:r] <- cov(var$res)
F_0 <- var$X[1L, ]                    # Initial factor estimate and covariance
P_0 <- ainv(diag((r*p)^2) - kronecker(A,A)) %*% unattrib(Q)
dim(P_0) <- c(r*p, r*p)

## Run standartized data through Kalman filter and smoother once
kfs_res <- SKFS(X, A, C, Q, R, F_0, P_0, FALSE)

## Two-step solution is state mean from the Kalman smoother
F_kal <- kfs_res$F_smooth[, 1:r, drop = FALSE]
colnames(F_kal) <- paste0("f", 1:r)

# See that this is equal to the Two-Step estimate by DFM()
expect_equal(F_kal, DFM(X, r, p, em.method = "none", pos.corr = FALSE)$F_2s)

# Same in two parts
kfs_res2 <- with(SKF(X, A, C, Q, R, F_0, P_0, FALSE), FIS(A, F, F_pred, P, P_pred))
F_kal2 <- kfs_res2$F_smooth[, 1:r, drop = FALSE]
colnames(F_kal2) <- paste0("f", 1:r)

expect_equal(F_kal, F_kal2)

