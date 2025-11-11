# (Fast) Stationary Kalman Filter and Smoother

(Fast) Stationary Kalman Filter and Smoother

## Usage

``` r
SKFS(X, A, C, Q, R, F_0, P_0, loglik = FALSE)
```

## Arguments

- X:

  numeric data matrix (\\T \times n\\).

- A:

  transition matrix (\\rp \times rp\\).

- C:

  observation matrix (\\n \times rp\\).

- Q:

  state covariance (\\rp \times rp\\).

- R:

  observation covariance (\\n \times n\\).

- F_0:

  initial state vector (\\rp \times 1\\).

- P_0:

  initial state covariance (\\rp \times rp\\).

- loglik:

  logical. Compute log-likelihood?

## Value

All results from
[`SKF`](https://sebkrantz.github.io/dfms/reference/SKF.md) and
[`FIS`](https://sebkrantz.github.io/dfms/reference/FIS.md), and
additionally a \\rp \times rp \times T\\ matrix `PPm_smooth`, which is
equal to the estimate of \\Cov(F^{smooth}\_t, F^{smooth}\_{t-1} \| T)\\
and needed for EM iterations. See 'Property 6.3: The Lag-One Covariance
Smoother' in Shumway & Stoffer (2017).

## References

Shumway, R. H., & Stoffer, D. S. (2017). Time Series Analysis and Its
Applications: With R Examples. Springer.

## See also

[`SKF`](https://sebkrantz.github.io/dfms/reference/SKF.md)
[`FIS`](https://sebkrantz.github.io/dfms/reference/FIS.md)
[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
library(collapse)
#> collapse 2.1.4, see ?`collapse-package` or ?`collapse-documentation`
#> 
#> Attaching package: ‘collapse’
#> The following object is masked from ‘package:vars’:
#> 
#>     B
#> The following object is masked from ‘package:stats’:
#> 
#>     D

## Two-Step factor estimates from monthly BM (2014) data
X <- fscale(diff(qM(BM14_M))) # Standardizing as KF has no intercept
r <- 5L # 5 Factors
p <- 3L # 3 Lags
n <- ncol(X)

## Initializing the Kalman Filter with PCA results
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

## Run standartized data through Kalman Filter and Smoother once
kfs_res <- SKFS(X, A, C, Q, R, F_0, P_0, FALSE)

## Two-step solution is state mean from the Kalman Smoother
F_kal <- kfs_res$F_smooth[, 1:r, drop = FALSE]
colnames(F_kal) <- paste0("f", 1:r)

## See that this is equal to the Two-Step estimate by DFM()
all.equal(F_kal, DFM(X, r, p, em.method = "none", pos.corr = FALSE)$F_2s)
#> [1] TRUE

## Same in two steps using SKF() and FIS()
kfs_res2 <- with(SKF(X, A, C, Q, R, F_0, P_0, FALSE), FIS(A, F, F_pred, P, P_pred))
F_kal2 <- kfs_res2$F_smooth[, 1:r, drop = FALSE]
colnames(F_kal2) <- paste0("f", 1:r)
all.equal(F_kal, F_kal2)
#> [1] TRUE

rm(X, r, p, n, X_imp, v, F_pc, C, res, R, var, A, Q, F_0, P_0, kfs_res, F_kal, kfs_res2, F_kal2)
```
