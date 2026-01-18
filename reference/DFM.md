# Estimate a Dynamic Factor Model

Efficient estimation of a Dynamic Factor Model via the EM Algorithm - on
stationary data with time-invariant system matrices and classical
assumptions, while permitting missing data.

## Usage

``` r
DFM(
  X,
  r,
  p = 1L,
  ...,
  idio.ar1 = FALSE,
  quarterly.vars = NULL,
  rQ = c("none", "diagonal", "identity"),
  rR = c("diagonal", "identity", "none"),
  em.method = c("auto", "DGR", "BM", "none"),
  min.iter = 25L,
  max.iter = 100L,
  tol = 1e-04,
  pos.corr = TRUE,
  check.increased = FALSE
)
```

## Arguments

- X:

  a `T x n` numeric data matrix or frame of stationary time series. May
  contain missing values. *Note* that data is internally standardized
  (scaled and centered) before estimation.

- r:

  integer. Number of factors.

- p:

  integer. Number of lags in factor VAR.

- ...:

  (optional) arguments to
  [`tsnarmimp`](https://sebkrantz.github.io/dfms/reference/tsnarmimp.md).
  The default settings impute internal missing values with a cubic
  spline and the edges with the median and a 3-period moving average.

- idio.ar1:

  logical. Model observation errors as AR(1) processes: \\e_t = \rho
  e\_{t-1} + v_t\\. *Note* that this substantially increases computation
  time, and is generally not needed if `n` is large (\>30). See
  theoretical vignette for details.

- quarterly.vars:

  character. Names of quarterly variables in `X` (if any). Monthly
  variables should be to the left of the quarterly variables in the data
  matrix and quarterly observations should be provided every 3rd period.

- rQ:

  character. Restrictions on the state (transition) covariance matrix
  (Q).

- rR:

  character. Restrictions on the observation (measurement) covariance
  matrix (R).

- em.method:

  character. The implementation of the Expectation Maximization
  Algorithm used. The options are:

  |          |     |                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |     |
  |----------|-----|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----|
  | `"auto"` |     | Automatic selection: `"BM"` if `anyNA(X)`, else `"DGR"`.                                                                                                                                                                                                                                                                                                                                                                                                             |     |
  | `"DGR"`  |     | The classical EM implementation of Doz, Giannone and Reichlin (2012). This implementation is efficient and quite robust, missing values are removed on a casewise basis in the Kalman Filter and Smoother, but not explicitly accounted for in EM iterations.                                                                                                                                                                                                        |     |
  | `"BM"`   |     | The modified EM algorithm of Banbura and Modugno (2014) which also accounts for missing data in the EM iterations. Optimal for datasets with systematically missing data e.g. datasets with ragged edges or series at different frequencies.                                                                                                                                                                                                                         |     |
  | `"none"` |     | Performs no EM iterations and just returns the Two-Step estimates from running the data through the Kalman Filter and Smoother once as in Doz, Giannone and Reichlin (2011) (the Kalman Filter is Initialized with system matrices obtained from a regression and VAR on PCA factor estimates). This yields significant performance gains over the iterative methods. Final system matrices are estimated by running a regression and a VAR on the smoothed factors. |     |

- min.iter:

  integer. Minimum number of EM iterations (to ensure a convergence
  path).

- max.iter:

  integer. Maximum number of EM iterations.

- tol:

  numeric. EM convergence tolerance.

- pos.corr:

  logical. Increase the likelihood that factors correlate positively
  with the data, by scaling the eigenvectors such that the principal
  components (used to initialize the Kalman Filter) co-vary positively
  with the row-means of the standardized data.

- check.increased:

  logical. Check if likelihood has increased. Passed to
  [`em_converged`](https://sebkrantz.github.io/dfms/reference/em_converged.md).
  If `TRUE`, the algorithm only terminates if convergence was reached
  with decreasing likelihood.

## Value

A list-like object of class 'dfm' with the following elements:

- `X_imp`:

  \\T \times n\\ matrix with the imputed and standardized (scaled and
  centered) data—after applying
  [`tsnarmimp`](https://sebkrantz.github.io/dfms/reference/tsnarmimp.md).
  It has attributes attached allowing for reconstruction of the original
  data:

  |                |     |                                                                                                                                     |     |
  |----------------|-----|-------------------------------------------------------------------------------------------------------------------------------------|-----|
  | `"stats"`      |     | is a \\n \times 5\\ matrix of summary statistics of class `"qsu"` (see [`qsu`](https://fastverse.org/collapse/reference/qsu.html)). |     |
  | `"missing"`    |     | is a \\T \times n\\ logical matrix indicating missing or infinite values in the original data (which are imputed in `X_imp`).       |     |
  | `"attributes"` |     | contains the [`attributes`](https://rdrr.io/r/base/attributes.html) of the original data input.                                     |     |
  | `"is.list"`    |     | is a logical value indicating whether the original data input was a list / data frame.                                              |     |

- `eigen`:

  `eigen(cov(X_imp))`.

- `F_pca`:

  \\T \times r\\ matrix of principal component factor estimates -
  `X_imp %*% eigen$vectors`.

- `P_0`:

  \\r \times r\\ initial factor covariance matrix estimate based on PCA
  results.

- `F_2s`:

  \\T \times r\\ matrix two-step factor estimates as in Doz, Giannone
  and Reichlin (2011) - obtained from running the data through the
  Kalman Filter and Smoother once, where the Filter is initialized with
  results from PCA.

- `P_2s`:

  \\r \times r \times T\\ covariance matrices of two-step factor
  estimates.

- `F_qml`:

  \\T \times r\\ matrix of quasi-maximum likelihood factor estimates -
  obtained by iteratively Kalman Filtering and Smoothing the factor
  estimates until EM convergence.

- `P_qml`:

  \\r \times r \times T\\ covariance matrices of QML factor estimates.

- `A`:

  \\r \times rp\\ factor transition matrix.

- `C`:

  \\n \times r\\ observation matrix.

- `Q`:

  \\r \times r\\ state (error) covariance matrix.

- `R`:

  \\n \times n\\ observation (error) covariance matrix.

- `e`:

  \\T \times n\\ estimates of observation errors \\\textbf{e}\_t\\. Only
  available if `idio.ar1 = TRUE`.

- `rho`:

  \\n \times 1\\ estimates of AR(1) coefficients (\\\rho\\) in
  observation errors: \\e_t = \rho e\_{t-1} + v_t\\. Only available if
  `idio.ar1 = TRUE`.

- `loglik`:

  vector of log-likelihoods - one for each EM iteration. The final value
  corresponds to the log-likelihood of the reported model.

- `tol`:

  The numeric convergence tolerance used.

- `converged`:

  single logical valued indicating whether the EM algorithm converged
  (within `max.iter` iterations subject to `tol`).

- `anyNA`:

  single logical valued indicating whether there were any (internal)
  missing values in the data (determined after removal of rows with too
  many missing values). If `FALSE`, `X_imp` is simply the original data
  in matrix form, and does not have the `"missing"` attribute attached.

- `rm.rows`:

  vector of any cases (rows) that were removed beforehand (subject to
  `max.missing` and `na.rm.method`). If no cases were removed the slot
  is `NULL`.

- `quarterly.vars`:

  names of the quarterly variables (if any).

- `em.method`:

  The EM method used.

- `call`:

  call object obtained from
  [`match.call()`](https://rdrr.io/r/base/match.call.html).

## Details

This function efficiently estimates a Dynamic Factor Model with the
following classical assumptions:

1.  Linearity

2.  Idiosynchratic measurement (observation) errors (*R* is diagonal)

3.  No direct relationship between series and lagged factors (*ceteris
    paribus* contemporaneous factors)

4.  No relationship between lagged error terms in the either measurement
    or transition equation (no serial correlation), unless explicitly
    modeled as AR(1) processes using `idio.ar1 = TRUE`.

Factors are allowed to evolve in a \\VAR(p)\\ process, and data is
internally standardized (scaled and centered) before estimation
(removing the need of intercept terms). By assumptions 1-4, this
translates into the following dynamic form:

\$\$\textbf{x}\_t = \textbf{C}\_0 \textbf{f}\_t + \textbf{e}\_t \\
\sim\\ N(\textbf{0}, \textbf{R})\$\$ \$\$\textbf{f}\_t = \sum\_{j=1}^p
\textbf{A}\_j \textbf{f}\_{t-j} + \textbf{u}\_t \\ \sim\\ N(\textbf{0},
\textbf{Q}\_0)\$\$

where the first equation is called the measurement or observation
equation and the second equation is called transition, state or process
equation, and

|                   |     |                                                                                                                                                                                                                                               |     |
|-------------------|-----|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----|
| \\n\\             |     | number of series in \\\textbf{x}\_t\\ (\\r\\ and \\p\\ as the arguments to `DFM`).                                                                                                                                                            |     |
| \\\textbf{x}\_t\\ |     | \\n \times 1\\ vector of observed series at time \\t\\: \\(x\_{1t}, \dots, x\_{nt})'\\. Some observations can be missing.                                                                                                                     |     |
| \\\textbf{f}\_t\\ |     | \\r \times 1\\ vector of factors at time \\t\\: \\(f\_{1t}, \dots, f\_{rt})'\\.                                                                                                                                                               |     |
| \\\textbf{C}\_0\\ |     | \\n \times r\\ measurement (observation) matrix.                                                                                                                                                                                              |     |
| \\\textbf{A}\_j\\ |     | \\r \times r\\ state transition matrix at lag \\j\\.                                                                                                                                                                                          |     |
| \\\textbf{Q}\_0\\ |     | \\r \times r\\ state covariance matrix.                                                                                                                                                                                                       |     |
| \\\textbf{R}\\    |     | \\n \times n\\ measurement (observation) covariance matrix. It is diagonal by assumption 2 that \\E\[\textbf{x}\_{it}\|\textbf{x}\_{-i,t},\textbf{x}\_{i,t-1}, \dots, \textbf{f}\_t, \textbf{f}\_{t-1}, \dots\] = \textbf{Cf}\_t \forall i\\. |     |

This model can be estimated using a classical form of the Kalman Filter
and the Expectation Maximization (EM) algorithm, after transforming it
to State-Space (stacked, VAR(1)) form:

\$\$\textbf{x}\_t = \textbf{C} \textbf{F}\_t + \textbf{e}\_t \\ \sim\\
N(\textbf{0}, \textbf{R})\$\$ \$\$\textbf{F}\_t = \textbf{A F}\_{t-1} +
\textbf{u}\_t \\ \sim\\ N(\textbf{0}, \textbf{Q})\$\$

where

|                   |     |                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |     |
|-------------------|-----|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----|
| \\n\\             |     | number of series in \\\textbf{x}\_t\\ (\\r\\ and \\p\\ as the arguments to `DFM`).                                                                                                                                                                                                                                                                                                                                                                              |     |
| \\\textbf{x}\_t\\ |     | \\n \times 1\\ vector of observed series at time \\t\\: \\(x\_{1t}, \dots, x\_{nt})'\\. Some observations can be missing.                                                                                                                                                                                                                                                                                                                                       |     |
| \\\textbf{F}\_t\\ |     | \\rp \times 1\\ vector of stacked factors at time \\t\\: \\(f\_{1t}, \dots, f\_{rt}, f\_{1,t-1}, \dots, f\_{r,t-1}, \dots, f\_{1,t-p}, \dots, f\_{r,t-p})'\\.                                                                                                                                                                                                                                                                                                   |     |
| \\\textbf{C}\\    |     | \\n \times rp\\ observation matrix. Only the first \\n \times r\\ terms are non-zero, by assumption 3 that \\E\[\textbf{x}\_t\|\textbf{F}\_t\] = E\[\textbf{x}\_t\|\textbf{f}\_t\]\\ (no relationship of observed series with lagged factors given contemporaneous factors).                                                                                                                                                                                    |     |
| \\\textbf{A}\\    |     | stacked \\rp \times rp\\ state transition matrix consisting of 3 parts: the top \\r \times rp\\ part provides the dynamic relationships captured by \\(\textbf{A}\_1, \dots, \textbf{A}\_p)\\ in the dynamic form, the terms `A[(r+1):rp, 1:(rp-r)]` constitute an \\(rp-r) \times (rp-r)\\ identity matrix mapping all lagged factors to their known values at times t. The remaining part `A[(rp-r+1):rp, (rp-r+1):rp]` is an \\r \times r\\ matrix of zeros. |     |
| \\\textbf{Q}\\    |     | \\rp \times rp\\ state covariance matrix. The top \\r \times r\\ part gives the contemporaneous relationships, the rest are zeros by assumption 4.                                                                                                                                                                                                                                                                                                              |     |
| \\\textbf{R}\\    |     | \\n \times n\\ observation covariance matrix. It is diagonal by assumption 2 and identical to \\\textbf{R}\\ as stated in the dynamic form.                                                                                                                                                                                                                                                                                                                     |     |

The filter is initialized with PCA estimates on the imputed dataset—see
[`SKFS`](https://sebkrantz.github.io/dfms/reference/SKFS.md) for a
complete code example.

## References

Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for
large approximate dynamic factor models based on Kalman filtering.
*Journal of Econometrics, 164*(1), 188-205.

Doz, C., Giannone, D., & Reichlin, L. (2012). A quasi-maximum likelihood
approach for large, approximate dynamic factor models. *Review of
Economics and Statistics, 94*(4), 1014-1024.

Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of
factor models on datasets with arbitrary pattern of missing data.
*Journal of Applied Econometrics, 29*(1), 133-160.

Stock, J. H., & Watson, M. W. (2016). Dynamic Factor Models,
Factor-Augmented Vector Autoregressions, and Structural Vector
Autoregressions in Macroeconomics. *Handbook of Macroeconomics, 2*,
415–525. https://doi.org/10.1016/bs.hesmac.2016.04.002

## See also

[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
library(magrittr)
library(xts)
library(vars)
#> Loading required package: MASS
#> Loading required package: strucchange
#> Loading required package: sandwich
#> Loading required package: urca
#> Loading required package: lmtest

# BM14 Replication Data. Constructing the database:
BM14 <- merge(BM14_M, BM14_Q)
BM14[, BM14_Models$log_trans] %<>% log()
BM14[, BM14_Models$freq == "M"] %<>% diff()
BM14[, BM14_Models$freq == "Q"] %<>% diff(3)

# \donttest{
### Small Model ---------------------------------------

# IC for number of factors
IC_small <- ICr(BM14[, BM14_Models$small], max.r = 5)
#> Missing values detected: imputing data with tsnarmimp() with default settings
plot(IC_small)

screeplot(IC_small)


# I take 2 factors. Now number of lags
VARselect(IC_small$F_pca[, 1:2])
#> $selection
#> AIC(n)  HQ(n)  SC(n) FPE(n) 
#>      6      5      2      6 
#> 
#> $criteria
#>                 1          2         3          4         5          6
#> AIC(n) -1.5028163 -1.7020957 -1.701093 -1.7021266 -1.782692 -1.7977899
#> HQ(n)  -1.4762557 -1.6578279 -1.639118 -1.6224446 -1.685303 -1.6826936
#> SC(n)  -1.4361151 -1.5909270 -1.545457 -1.5020229 -1.538121 -1.5087511
#> FPE(n)  0.2225028  0.1823018  0.182486  0.1822997  0.168192  0.1656763
#>                 7          8          9         10
#> AIC(n) -1.7892503 -1.7887748 -1.7793434 -1.7708400
#> HQ(n)  -1.6564469 -1.6382643 -1.6111258 -1.5849153
#> SC(n)  -1.4557440 -1.4108011 -1.3569022 -1.3039313
#> FPE(n)  0.1671036  0.1671913  0.1687862  0.1702408
#> 

# Estimating the model with 2 factors and 3 lags
dfm_small <- DFM(BM14[, BM14_Models$small], r = 2, p = 3,
    quarterly.vars = BM14_Models %$% series[freq == "Q" & small])
#> Converged after 26 iterations.

# Inspecting the model
summary(dfm_small)
#> Mixed Frequency Dynamic Factor Model
#> n = 14, nm = 10, nq = 4, T = 356, r = 2, p = 3
#> %NA = 38.3628, %NAm = 26.3202
#> 
#> Call:  DFM(X = BM14[, BM14_Models$small], r = 2, p = 3, quarterly.vars = BM14_Models %$%      series[freq == "Q" & small])
#> 
#> Summary Statistics of Factors [F]
#>       N     Mean  Median      SD      Min     Max
#> f1  356  -0.1605  0.1274  2.0755  -9.4638  3.1838
#> f2  356   0.1049  0.0463  1.3905  -4.2535  5.3709
#> 
#> Factor Transition Matrix [A]
#>     L1.f1   L1.f2   L2.f1   L2.f2   L3.f1    L3.f2
#> f1 1.2824 -0.1946 -0.2100  0.2179 -0.1402 -0.05373
#> f2 0.1709  0.5829  0.2514 -0.2866 -0.6082  0.41312
#> 
#> Factor Covariance Matrix [cov(F)]
#>           f1        f2
#> f1   4.3077   -0.5272*
#> f2  -0.5272*   1.9334 
#> 
#> Factor Transition Error Covariance Matrix [Q]
#>        u1     u2
#> u1 0.3026 0.2191
#> u2 0.2191 0.4764
#> 
#> Observation Matrix [C]
#>                             f1      f2
#> ip_tot_cstr             0.2516  0.1905
#> new_cars                0.0265  0.0290
#> orders                  0.1725  0.1659
#> ret_turnover_defl       0.0533  0.0067
#> ecs_ec_sent_ind         0.2650  0.4645
#> pms_pmi                 0.1296  0.4786
#> urx                    -0.3799  0.1838
#> extra_ea_trade_exp_val  0.0601  0.0865
#> euro325                 0.1307  0.2383
#> raw_mat                 0.0958  0.1625
#> gdp                     0.0434  0.0217
#> empl                    0.0437 -0.0249
#> capacity                0.0392  0.0116
#> gdp_us                  0.0268  0.0437
#> 
#> Observation Error Covariance Matrix [diag(R) - Restricted]
#>            ip_tot_cstr               new_cars                 orders 
#>                 0.6269                 0.9909                 0.7924 
#>      ret_turnover_defl        ecs_ec_sent_ind                pms_pmi 
#>                 0.9845                 0.3332                 0.3706 
#>                    urx extra_ea_trade_exp_val                euro325 
#>                 0.1681                 0.9688                 0.8183 
#>                raw_mat                    gdp                   empl 
#>                 0.9096                 0.3599                 0.1326 
#>               capacity                 gdp_us 
#>                 0.4659                 0.6048 
#> 
#> Observation Residual Covariance Matrix [cov(resid(DFM))]
#>                        ip_tot_cstr  new_cars    orders ret_turnover_defl
#> ip_tot_cstr                0.6016   -0.0155    0.1950*          -0.0191 
#> new_cars                  -0.0155    0.9946   -0.0080            0.1262 
#> orders                     0.1950*  -0.0080    0.7821            0.0324 
#> ret_turnover_defl         -0.0191    0.1262    0.0324            0.9864 
#> ecs_ec_sent_ind           -0.0209   -0.0282   -0.0861*          -0.0101 
#> pms_pmi                   -0.0687   -0.0472   -0.0187            0.0022 
#> urx                        0.0113   -0.0158    0.0547*          -0.0098 
#> extra_ea_trade_exp_val     0.0697    0.0963    0.0742           -0.0190 
#> euro325                   -0.0373    0.0556   -0.0904            0.0107 
#> raw_mat                    0.0921    0.1067    0.0067            0.0323 
#> gdp                        0.0720   -0.0844    0.0107           -0.0055 
#> empl                      -0.0031   -0.1102   -0.0405           -0.0695 
#> capacity                  -0.0066   -0.0474   -0.0520           -0.1209 
#> gdp_us                    -0.0740   -0.0115    0.0609           -0.1050 
#>                        ecs_ec_sent_ind   pms_pmi       urx
#> ip_tot_cstr                   -0.0209   -0.0687    0.0113 
#> new_cars                      -0.0282   -0.0472   -0.0158 
#> orders                        -0.0861*  -0.0187    0.0547*
#> ret_turnover_defl             -0.0101    0.0022   -0.0098 
#> ecs_ec_sent_ind                0.2434   -0.0836*  -0.0114 
#> pms_pmi                       -0.0836*   0.2982    0.0106 
#> urx                           -0.0114    0.0106    0.1409 
#> extra_ea_trade_exp_val        -0.0697*  -0.0419    0.0077 
#> euro325                       -0.0358   -0.0396   -0.0014 
#> raw_mat                       -0.0433    0.0339    0.0248 
#> gdp                           -0.0258    0.0604    0.0200 
#> empl                           0.0444    0.1233    0.0351 
#> capacity                      -0.0473    0.0574    0.0304 
#> gdp_us                        -0.0519    0.0210    0.0555 
#>                        extra_ea_trade_exp_val   euro325   raw_mat       gdp
#> ip_tot_cstr                           0.0697   -0.0373    0.0921    0.0720 
#> new_cars                              0.0963    0.0556    0.1067   -0.0844 
#> orders                                0.0742   -0.0904    0.0067    0.0107 
#> ret_turnover_defl                    -0.0190    0.0107    0.0323   -0.0055 
#> ecs_ec_sent_ind                      -0.0697*  -0.0358   -0.0433   -0.0258 
#> pms_pmi                              -0.0419   -0.0396    0.0339    0.0604 
#> urx                                   0.0077   -0.0014    0.0248    0.0200 
#> extra_ea_trade_exp_val                0.9673    0.0112   -0.0329    0.0482 
#> euro325                               0.0112    0.7982    0.0097    0.0341 
#> raw_mat                              -0.0329    0.0097    0.8984    0.0022 
#> gdp                                   0.0482    0.0341    0.0022    0.8794 
#> empl                                 -0.0869    0.0905    0.0372    0.4779*
#> capacity                             -0.0428    0.0477   -0.1930*   0.5392*
#> gdp_us                               -0.0780    0.1076    0.0037    0.3432*
#>                             empl  capacity    gdp_us
#> ip_tot_cstr             -0.0031   -0.0066   -0.0740 
#> new_cars                -0.1102   -0.0474   -0.0115 
#> orders                  -0.0405   -0.0520    0.0609 
#> ret_turnover_defl       -0.0695   -0.1209   -0.1050 
#> ecs_ec_sent_ind          0.0444   -0.0473   -0.0519 
#> pms_pmi                  0.1233    0.0574    0.0210 
#> urx                      0.0351    0.0304    0.0555 
#> extra_ea_trade_exp_val  -0.0869   -0.0428   -0.0780 
#> euro325                  0.0905    0.0477    0.1076 
#> raw_mat                  0.0372   -0.1930*   0.0037 
#> gdp                      0.4779*   0.5392*   0.3432*
#> empl                     0.8218    0.4565*   0.1683*
#> capacity                 0.4565*   0.9166    0.2969*
#> gdp_us                   0.1683*   0.2969*   0.9371 
#> 
#> Residual AR(1) Serial Correlation
#>            ip_tot_cstr               new_cars                 orders 
#>               -0.38749               -0.40829               -0.41500 
#>      ret_turnover_defl        ecs_ec_sent_ind                pms_pmi 
#>               -0.49326               -0.09350                0.08832 
#>                    urx extra_ea_trade_exp_val                euro325 
#>               -0.09171               -0.52651                0.24433 
#>                raw_mat                    gdp                   empl 
#>                0.27610                     NA                     NA 
#>               capacity                 gdp_us 
#>                     NA                     NA 
#> 
#> Summary of Residual AR(1) Serial Correlations
#>    N     Mean   Median      SD      Min     Max
#>   10  -0.1807  -0.2405  0.3062  -0.5265  0.2761
#> 
#> Goodness of Fit: R-Squared
#>            ip_tot_cstr               new_cars                 orders 
#>                 0.3984                 0.0054                 0.2179 
#>      ret_turnover_defl        ecs_ec_sent_ind                pms_pmi 
#>                 0.0136                 0.7566                 0.7018 
#>                    urx extra_ea_trade_exp_val                euro325 
#>                 0.8591                 0.0327                 0.2018 
#>                raw_mat                    gdp                   empl 
#>                 0.1016                 0.1206                 0.1782 
#>               capacity                 gdp_us 
#>                 0.0834                 0.0629 
#> 
#> Summary of Individual R-Squared's
#>    N    Mean  Median      SD     Min     Max
#>   14  0.2667  0.1494  0.2939  0.0054  0.8591
plot(dfm_small)  # Factors and data

plot(dfm_small, method = "all", type = "individual") # Factor estimates

plot(dfm_small, type = "residual") # Residuals from factor predictions


# 10 periods ahead forecast
plot(predict(dfm_small), xlim = c(300, 370))



### Medium-Sized Model ---------------------------------

# IC for number of factors
IC_medium <- ICr(BM14[, BM14_Models$medium])
#> Missing values detected: imputing data with tsnarmimp() with default settings
plot(IC_medium)

screeplot(IC_medium)


# I take 3 factors. Now number of lags
VARselect(IC_medium$F_pca[, 1:3])
#> $selection
#> AIC(n)  HQ(n)  SC(n) FPE(n) 
#>      7      2      1      7 
#> 
#> $criteria
#>               1        2        3        4        5        6        7        8
#> AIC(n) 1.539267 1.492781 1.498671 1.509014 1.486111 1.479450 1.468524 1.490378
#> HQ(n)  1.592388 1.585744 1.631474 1.681659 1.698597 1.731776 1.760691 1.822386
#> SC(n)  1.672669 1.726236 1.832177 1.942573 2.019721 2.113112 2.202238 2.324143
#> FPE(n) 4.661186 4.449528 4.475952 4.522752 4.420750 4.391987 4.345059 4.442132
#>               9       10
#> AIC(n) 1.500199 1.534401
#> HQ(n)  1.872048 1.946092
#> SC(n)  2.434016 2.568271
#> FPE(n) 4.487352 4.645257
#> 

# Estimating the model with 3 factors and 3 lags
dfm_medium <- DFM(BM14[, BM14_Models$medium], r = 3, p = 3,
    quarterly.vars = BM14_Models %$% series[freq == "Q" & medium])
#> Converged after 26 iterations.

# Inspecting the model
summary(dfm_medium)
#> Mixed Frequency Dynamic Factor Model
#> n = 48, nm = 39, nq = 9, T = 356, r = 3, p = 3
#> %NA = 32.6135, %NAm = 24.0781
#> 
#> Call:  DFM(X = BM14[, BM14_Models$medium], r = 3, p = 3, quarterly.vars = BM14_Models %$%      series[freq == "Q" & medium])
#> 
#> Summary Statistics of Factors [F]
#>       N     Mean   Median      SD       Min     Max
#> f1  356  -0.0985   0.4918  3.1283  -19.3582  5.3549
#> f2  356   0.0341   -0.224  2.2831    -6.462  8.3411
#> f3  356  -0.0172  -0.1399  1.5726   -6.7371  4.3577
#> 
#> Factor Transition Matrix [A]
#>       L1.f1    L1.f2  L1.f3   L2.f1    L2.f2    L2.f3   L3.f1  L3.f2    L3.f3
#> f1  0.93895 -0.19728 0.1715 0.05318  0.06486  0.01496 -0.1298 0.0859 -0.02337
#> f2 -0.30923  0.75786 0.1641 0.38513 -0.38768 -0.02307 -0.2750 0.4086  0.02054
#> f3  0.06881  0.03244 0.4077 0.36873 -0.32608 -0.17889 -0.4565 0.3461  0.11805
#> 
#> Factor Covariance Matrix [cov(F)]
#>           f1        f2        f3
#> f1   9.7864   -0.6567    0.3212 
#> f2  -0.6567    5.2124    0.0848 
#> f3   0.3212    0.0848    2.4729 
#> 
#> Factor Transition Error Covariance Matrix [Q]
#>         u1      u2      u3
#> u1  1.8586  1.3143 -0.3865
#> u2  1.3143  1.8104 -0.5786
#> u3 -0.3865 -0.5786  1.8320
#> 
#> Summary of Residual AR(1) Serial Correlations
#>    N     Mean   Median      SD      Min     Max
#>   39  -0.0901  -0.0699  0.2595  -0.5508  0.3205
#> 
#> Summary of Individual R-Squared's
#>    N    Mean  Median      SD     Min     Max
#>   48  0.2894   0.155  0.2852  0.0086  0.9362
plot(dfm_medium)  # Factors and data

plot(dfm_medium, method = "all", type = "individual") # Factor estimates

plot(dfm_medium, type = "residual") # Residuals from factor predictions


# 10 periods ahead forecast
plot(predict(dfm_medium), xlim = c(300, 370))



### Large Model ---------------------------------

# IC for number of factors
IC_large <- ICr(BM14)
#> Missing values detected: imputing data with tsnarmimp() with default settings
plot(IC_large)

screeplot(IC_large)


# I take 6 factors. Now number of lags
VARselect(IC_large$F_pca[, 1:6])
#> $selection
#> AIC(n)  HQ(n)  SC(n) FPE(n) 
#>      6      1      1      6 
#> 
#> $criteria
#>                 1          2          3          4          5          6
#> AIC(n)   6.020768   5.926947   5.883307   5.889824   5.897327   5.845883
#> HQ(n)    6.206692   6.272235   6.387960   6.553841   6.720708   6.828628
#> SC(n)    6.487677   6.794063   7.150631   7.557355   7.965066   8.313830
#> FPE(n) 411.908362 375.087393 359.232802 361.888929 365.116976 347.515300
#>                 7          8          9         10
#> AIC(n)   5.908082   5.963699   5.982969   6.111831
#> HQ(n)    7.050191   7.265172   7.443806   7.732032
#> SC(n)    8.776235   9.232061   9.651538  10.180607
#> FPE(n) 370.862536 393.546431 403.139879 461.355079
#> 

# Estimating the model with 6 factors and 3 lags
dfm_large <- DFM(BM14, r = 6, p = 3,
    quarterly.vars = BM14_Models %$% series[freq == "Q"])
#> Converged after 39 iterations.

# Inspecting the model
summary(dfm_large)
#> Mixed Frequency Dynamic Factor Model
#> n = 101, nm = 92, nq = 9, T = 356, r = 6, p = 3
#> %NA = 29.7363, %NAm = 25.8366
#> 
#> Call:  DFM(X = BM14, r = 6, p = 3, quarterly.vars = BM14_Models %$%      series[freq == "Q"])
#> 
#> Summary Statistics of Factors [F]
#>       N     Mean   Median      SD       Min      Max
#> f1  356   -0.107   0.6573  4.8387  -23.1251  14.7986
#> f2  356  -0.2216   0.1847  3.4756  -13.9892  17.5939
#> f3  356  -0.0085    0.038   2.538    -9.952    7.023
#> f4  356   0.1225   0.2815  3.1901  -17.2132  11.1281
#> f5  356  -0.0816  -0.0857  2.6977  -11.8225  11.6105
#> f6  356  -0.0021  -0.0044  2.2638   -8.0332  14.3884
#> 
#> Factor Transition Matrix [A]
#>      L1.f1    L1.f2     L1.f3     L1.f4    L1.f5     L1.f6      L2.f1    L2.f2
#> f1  0.4146 -0.31547  0.535594 -0.409677  0.24935 -0.042001 -0.0013023  0.10535
#> f2 -0.1181  0.42170 -0.132113 -0.044417  0.11608  0.126305  0.0865712 -0.02060
#> f3  0.2563  0.07508  0.387405  0.069762 -0.13764 -0.102988 -0.1223568  0.01331
#> f4 -0.2962 -0.12730  0.098611 -0.104417  0.48309 -0.037316  0.0466138 -0.10279
#> f5  0.3225  0.08940 -0.008011  0.344206  0.09136  0.009044 -0.0134199  0.07329
#> f6 -0.1849  0.10457 -0.088462 -0.009701 -0.07510  0.362577  0.0007091 -0.06256
#>       L2.f3    L2.f4     L2.f5    L2.f6    L3.f1    L3.f2    L3.f3    L3.f4
#> f1 -0.18944 -0.17725 -0.036033  0.48920  0.32166 -0.06934  0.06825  0.05621
#> f2  0.04108 -0.06438  0.172166 -0.08944  0.15363  0.22166  0.11864 -0.13273
#> f3 -0.27139  0.19914 -0.010737  0.02744 -0.04957 -0.08776  0.10569  0.12061
#> f4  0.04184 -0.17465 -0.003223  0.16731  0.07214 -0.11843 -0.03447  0.12219
#> f5 -0.08275 -0.01601 -0.107952 -0.10232 -0.17649  0.02642 -0.00477  0.01682
#> f6  0.03254 -0.04805  0.086758 -0.10416  0.09119  0.06214  0.13715 -0.12928
#>       L3.f5    L3.f6
#> f1  0.15754  0.05201
#> f2  0.01799  0.09386
#> f3 -0.11424  0.01713
#> f4 -0.06460  0.07272
#> f5  0.03395 -0.14320
#> f6  0.05089  0.02547
#> 
#> Factor Covariance Matrix [cov(F)]
#>           f1        f2        f3        f4        f5        f6
#> f1  23.4128    0.7031   -0.2198    3.7881*  -2.5178*  -1.9928*
#> f2   0.7031   12.0798   -1.6485*  -0.4811   -0.6676    2.4008*
#> f3  -0.2198   -1.6485*   6.4415   -1.4691*   1.4553*  -0.9491*
#> f4   3.7881*  -0.4811   -1.4691*  10.1768   -4.6435*  -0.4141 
#> f5  -2.5178*  -0.6676    1.4553*  -4.6435*   7.2774   -0.4585 
#> f6  -1.9928*   2.4008*  -0.9491*  -0.4141   -0.4585    5.1250 
#> 
#> Factor Transition Error Covariance Matrix [Q]
#>         u1      u2      u3      u4      u5      u6
#> u1 10.2115  0.2060 -1.4420  2.4061 -2.2408  0.2968
#> u2  0.2060  5.8639 -0.0814  0.8044 -0.6562  0.1438
#> u3 -1.4420 -0.0814  4.3219 -0.5280  0.4627  0.1162
#> u4  2.4061  0.8044 -0.5280  4.9253 -1.6502  0.0570
#> u5 -2.2408 -0.6562  0.4627 -1.6502  4.4696 -0.0469
#> u6  0.2968  0.1438  0.1162  0.0570 -0.0469  3.1868
#> 
#> Summary of Residual AR(1) Serial Correlations
#>    N     Mean   Median      SD      Min     Max
#>   92  -0.0673  -0.0953  0.2657  -0.4957  0.6614
#> 
#> Summary of Individual R-Squared's
#>     N    Mean  Median      SD     Min     Max
#>   101  0.4297  0.3892  0.3071  0.0121  0.9987
plot(dfm_large)  # Factors and data

# plot(dfm_large, method = "all", type = "individual") # Factor estimates
plot(dfm_large, type = "residual") # Residuals from factor predictions


# 10 periods ahead forecast
plot(predict(dfm_large), xlim = c(300, 370))



### Mixed-Frequency Model with AR(1) Idiosyncratic Errors ---------

# Estimate model with AR(1) observation errors
# This models e(t) = rho * e(t-1) + v(t) for each series
dfm_large_ar1 <- DFM(BM14, r = 6, p = 3, idio.ar1 = TRUE,
    quarterly.vars = BM14_Models %$% series[freq == "Q"])
#> Converged after 27 iterations.

# Model summary shows AR(1) coefficients
summary(dfm_large_ar1)
#> Mixed Frequency Dynamic Factor Model
#> n = 101, nm = 92, nq = 9, T = 356, r = 6, p = 3
#> %NA = 29.7363, %NAm = 25.8366
#>    with AR(1) errors: mean(abs(rho)) = 0.272 
#> 
#> Call:  DFM(X = BM14, r = 6, p = 3, idio.ar1 = TRUE, quarterly.vars = BM14_Models %$%      series[freq == "Q"])
#> 
#> Summary Statistics of Factors [F]
#>       N     Mean   Median      SD       Min      Max
#> f1  356  -0.0932   0.4356  4.1258  -26.3399   7.2689
#> f2  356   0.0141   0.3348  2.8144  -13.6388  14.9466
#> f3  356  -0.0904   0.0623  2.4063  -10.6054   7.7731
#> f4  356   0.0655   0.1567  2.4883   -9.7491     9.14
#> f5  356  -0.0532  -0.0793  2.2006   -7.0683   7.9831
#> f6  356   0.0785   0.1214  1.7718    -9.055   4.8855
#> 
#> Factor Transition Matrix [A]
#>       L1.f1    L1.f2     L1.f3    L1.f4     L1.f5    L1.f6    L2.f1     L2.f2
#> f1  0.36294 -0.33739  0.390802 -0.16640  0.015708  0.17678  0.16782 -0.034706
#> f2 -0.04303  0.20884 -0.072546  0.15323  0.009384 -0.03609  0.08371  0.035102
#> f3  0.16766  0.12552  0.333441  0.03715 -0.056107 -0.13851 -0.07633 -0.049216
#> f4 -0.14238 -0.02027  0.104313 -0.03705  0.409972  0.05844  0.02029 -0.025738
#> f5  0.18039  0.01489 -0.035480  0.28023  0.120216 -0.06277 -0.01386  0.017727
#> f6 -0.06093  0.04987  0.009273  0.05994 -0.011549  0.19555 -0.04598  0.006984
#>        L2.f3     L2.f4     L2.f5    L2.f6    L3.f1    L3.f2    L3.f3     L3.f4
#> f1 -0.061362 -0.084291  0.066044  0.13487  0.29337 -0.02864  0.07579 -0.025627
#> f2 -0.079912 -0.011331  0.004003  0.06906  0.22750  0.28135  0.10763 -0.008841
#> f3 -0.150995  0.187179  0.072732 -0.17937 -0.06504 -0.08119  0.10897  0.058212
#> f4  0.009792 -0.147390  0.052710  0.04956 -0.08077 -0.11240 -0.10943  0.115191
#> f5 -0.051458 -0.089720 -0.112780  0.06085 -0.05774  0.04967  0.06603  0.012034
#> f6 -0.068374  0.007569 -0.049427  0.04002  0.06585  0.06722  0.05816 -0.026570
#>       L3.f5    L3.f6
#> f1  0.02247 -0.05281
#> f2 -0.01189  0.13155
#> f3 -0.11803 -0.04513
#> f4 -0.13004  0.06949
#> f5  0.07094 -0.15219
#> f6  0.04678  0.07416
#> 
#> Factor Covariance Matrix [cov(F)]
#>           f1        f2        f3        f4        f5        f6
#> f1  17.0224    0.9532   -0.1667   -0.0775   -0.3385   -0.2222 
#> f2   0.9532    7.9206   -0.2488   -1.0206*   0.1217   -0.3593 
#> f3  -0.1667   -0.2488    5.7905    0.0333   -0.0485   -0.0480 
#> f4  -0.0775   -1.0206*   0.0333    6.1916   -1.3488*   0.3836 
#> f5  -0.3385    0.1217   -0.0485   -1.3488*   4.8426    0.1443 
#> f6  -0.2222   -0.3593   -0.0480    0.3836    0.1443    3.1392 
#> 
#> Factor Transition Error Covariance Matrix [Q]
#>         u1      u2      u3      u4      u5      u6
#> u1  6.6526 -0.5783 -1.1405  1.3603 -1.3792  0.6104
#> u2 -0.5783  5.3191  0.0332 -0.1322 -0.2439 -0.5950
#> u3 -1.1405  0.0332  4.4231  0.0626 -0.0724  0.1648
#> u4  1.3603 -0.1322  0.0626  4.7055 -0.9944  0.4509
#> u5 -1.3792 -0.2439 -0.0724 -0.9944  4.1330  0.1573
#> u6  0.6104 -0.5950  0.1648  0.4509  0.1573  2.9458
#> 
#> Summary of Residual AR(1) Serial Correlations
#>     N     Mean   Median      SD      Min     Max
#>   101  -0.0581  -0.1326  0.3254  -0.6311  0.9075
#> 
#> Summary of Individual R-Squared's
#>     N    Mean  Median      SD     Min     Max
#>   101  0.6087  0.6793  0.2909  0.0039  0.9975

# Access AR(1) coefficients (rho) for each series
head(dfm_large_ar1$rho)
#>       ip_total    ip_tot_cstr ip_tot_cstr_en      ip_constr    ip_im_goods 
#>     -0.4565217     -0.4189399     -0.6310923     -0.3832824     -0.1150151 
#>     ip_capital 
#>     -0.3351250 

# Access estimated observation errors
head(dfm_large_ar1$e)
#>           ip_total   ip_tot_cstr ip_tot_cstr_en     ip_constr    ip_im_goods
#> [1,]  4.090700e-41  2.773400e-43  -9.284405e-25  1.068275e-49  1.229385e-109
#> [2,] -3.169609e-41 -1.163997e-43   7.234494e-25 -1.071486e-49 -1.413799e-110
#> [3,]  4.299249e-41  4.926788e-44  -6.744676e-25  2.137958e-49  1.610407e-111
#> [4,] -8.210496e-41 -2.184188e-44   7.709324e-25 -5.325980e-49 -4.893789e-113
#> [5,]  1.743392e-40  1.201861e-44  -1.033646e-24  1.379910e-48 -1.179288e-111
#> [6,] -3.793705e-40 -1.188133e-44   1.519262e-24 -3.596541e-48  1.043790e-110
#>         ip_capital     ip_d_cstr    ip_nd_cons         ip_en       ip_en_2
#> [1,] -1.944303e-56  1.888720e-84  4.816108e-66 -9.555744e-56 -1.944248e-66
#> [2,]  8.178211e-57 -5.369288e-85 -3.094479e-66  1.295765e-55  2.187136e-66
#> [3,] -7.701158e-57  9.752601e-85  7.051130e-66 -3.266343e-55 -6.486240e-66
#> [4,]  1.738260e-56 -4.665112e-84 -2.394640e-65  9.346773e-55  2.281459e-65
#> [5,] -4.999318e-56  2.390846e-83  8.478264e-65 -2.718771e-54 -8.145441e-65
#> [6,]  1.485491e-55 -1.228627e-82 -3.011926e-64  7.923738e-54  2.911579e-64
#>           ip_manuf     ip_metals  ip_chemicals   ip_electric  ip_machinery
#> [1,]  3.811841e-28  1.061574e-61  2.808696e-57 -1.379006e-49 -9.142109e-37
#> [2,] -2.798538e-28 -4.290610e-62 -2.626191e-57  1.129823e-49  5.296934e-37
#> [3,]  2.611274e-28  4.675522e-62  5.939577e-57 -1.982737e-49 -4.251618e-37
#> [4,] -3.194775e-28 -1.237244e-61 -1.715952e-56  4.769731e-49  5.453626e-37
#> [5,]  4.721270e-28  3.941849e-61  5.122161e-56 -1.220941e-48 -9.538318e-37
#> [6,] -7.641334e-28 -1.281107e-60 -1.534681e-55  3.155888e-48  1.866478e-36
#>           ip_paper    ip_plastic      new_cars        orders ret_turnover_defl
#> [1,] -6.441647e-94 -6.682650e-44 -5.749264e-45 -1.469969e-66       -3.86140700
#> [2,]  2.642492e-93  5.150991e-44  4.387543e-45  9.461044e-67        2.01426869
#> [3,] -1.558201e-92 -7.472124e-44 -6.504336e-45 -1.156951e-66       -0.97525169
#> [4,]  9.303855e-92  1.538219e-43  1.377787e-44  2.266244e-66        0.04166604
#> [5,] -5.557195e-91 -3.479769e-43 -3.197471e-44 -5.135422e-66       -0.51489075
#> [6,]  3.319347e-90  8.024085e-43  7.552159e-44  1.199259e-65        1.81347305
#>      ecs_ec_sent_ind  ecs_ind_conf ecs_ind_order_book ecs_ind_stocks
#> [1,]    7.019681e-40  2.458348e-58      -2.324592e-53  -2.327511e-29
#> [2,]   -2.529104e-39 -3.460365e-58       3.951806e-53   2.307864e-29
#> [3,]    1.099037e-38  2.899652e-57      -2.799512e-52  -5.547478e-29
#> [4,]   -4.828063e-38 -2.601194e-56       2.108370e-51   1.662147e-28
#> [5,]    2.122165e-37  2.335500e-55      -1.589624e-50  -5.116898e-28
#> [6,]   -9.328206e-37 -2.096969e-54       1.198534e-49   1.579794e-27
#>      ecs_ind_prod_exp ecs_ind_prod_rec_m ecs_ind_x_orders ecs_ind_empl_exp
#> [1,]    -1.963344e-42      -7.053766e-17     2.581352e-56    -5.788957e-51
#> [2,]     1.781038e-42       5.036705e-17    -1.823120e-56     3.407008e-50
#> [3,]    -7.377403e-42      -5.294718e-17     1.386374e-55    -2.242466e-49
#> [4,]     3.691013e-41       7.944348e-17    -1.232320e-54     1.480006e-48
#> [5,]    -1.861997e-40      -1.418243e-16     1.097726e-53    -9.768506e-48
#> [6,]     9.396237e-40       2.682671e-16    -9.778589e-53     6.447533e-47
#>      ecs_cons_conf ecs_cons_sit_over_next_12 ecs_cons_exp_unempl
#> [1,] -6.746136e-74             -8.441318e-59       7.186846e-111
#> [2,]  1.093414e-72              7.454530e-58      -8.399174e-112
#> [3,] -1.792053e-71             -6.769903e-57       5.646915e-110
#> [4,]  2.937213e-70              6.150267e-56      -4.278867e-108
#> [5,] -4.814155e-69             -5.587369e-55       3.242321e-106
#> [6,]  7.890504e-68              5.075990e-54      -2.456875e-104
#>      ecs_cons_gen_last_12m ecs_cstr_conf ecs_cstr_order_books ecs_cstr_empl_exp
#> [1,]         -1.558279e-59 -3.864650e-39        -3.206632e-32     -3.142529e-34
#> [2,]         -1.801058e-59  4.273319e-39         2.448859e-32      2.393847e-34
#> [3,]         -1.614725e-58 -1.614561e-38        -5.921301e-32     -6.334039e-34
#> [4,]         -1.569366e-57  7.133020e-38         1.962234e-31      2.268084e-33
#> [5,]         -1.526639e-56 -3.178655e-37        -6.721948e-31     -8.345304e-33
#> [6,]         -1.485088e-55  1.417108e-36         2.309332e-30      3.076863e-32
#>      ecs_cstr_prod_recent ecs_ret_tr_conf ecs_ret_tr_bus_sit ecs_ret_tr_stocks
#> [1,]        -2.074378e-53    3.631713e-39       1.110684e-40     -1.086935e-23
#> [2,]         8.861033e-53   -1.719347e-39      -6.909123e-41      9.565073e-24
#> [3,]        -6.762395e-52    4.616078e-39       2.365973e-40     -1.649752e-23
#> [4,]         5.230494e-51   -2.042420e-38      -1.121462e-39      3.763642e-23
#> [5,]        -4.046531e-50    9.335982e-38       5.406575e-39     -9.118502e-23
#> [6,]         3.130580e-49   -4.274274e-37      -2.608431e-38      2.232554e-22
#>      ecs_ret_tr_exp_bus ecs_ret_tr_empl  ecs_serv_conf ecs_serv_empl_exp
#> [1,]       1.768988e-43   -8.626704e-38 -1.746649e-150    -7.221803e-182
#> [2,]      -2.671458e-43    1.183339e-37  2.045560e-150     8.941669e-183
#> [3,]       1.254029e-42   -3.936167e-37 -1.202534e-149    -1.207501e-183
#> [4,]      -6.449884e-42    1.477913e-36  7.891663e-149     9.738593e-184
#> [5,]       3.329385e-41   -5.599817e-36 -5.192910e-148    -6.789456e-183
#> [6,]      -1.718839e-40    2.123122e-35  3.417277e-147     5.477853e-182
#>      pms_comp_output  pms_comp_empl       pms_pmi pms_manuf_empl
#> [1,]               0 -3.446477e-217 2.425447e-271 -1.875658e-234
#> [2,]               0  2.708394e-218 4.646938e-270  9.584543e-236
#> [3,]               0 -2.128376e-219 8.899469e-269 -4.897667e-237
#> [4,]               0  1.672572e-220 1.704358e-267  2.502690e-238
#> [5,]               0 -1.314381e-221 3.264056e-266 -1.278866e-239
#> [6,]               0  1.032898e-222 6.251069e-265  6.534958e-241
#>      pms_manuf_output pms_manuf_product   pms_serv_out  pms_serv_empl
#> [1,]   -1.385982e-159    -8.746639e-158 -5.438809e-313 -1.475161e-142
#> [2,]    1.957084e-160    -1.311423e-158  1.862736e-314  3.138414e-143
#> [3,]   -2.763512e-161    -1.966276e-159 -6.379680e-316 -6.677037e-144
#> [4,]    3.902232e-162    -2.948127e-160  2.185075e-317  1.420752e-144
#> [5,]   -5.510168e-163    -4.420261e-161 -7.775506e-319 -3.032455e-145
#> [6,]    7.780664e-164    -6.627498e-162  8.787303e-319  6.912067e-146
#>      pms_serv_new_bus pms_serv_product          urx    empl_total  empl_tot_xc
#> [1,]                0   -5.532815e-122 1.527306e-50 -2.189655e-08 1.664538e-06
#> [2,]                0   -1.507222e-122 6.916613e-51 -2.105954e-08 1.666183e-06
#> [3,]                0   -5.383310e-123 3.140323e-51 -2.054119e-08 1.689056e-06
#> [4,]                0   -6.611946e-123 1.443543e-51 -2.033363e-08 1.733449e-06
#> [5,]                0   -2.124985e-122 7.026725e-52 -2.043374e-08 1.799928e-06
#> [6,]                0   -7.898323e-122 4.271091e-52 -2.084302e-08 1.889340e-06
#>          empl_cstr   empl_manuf extra_ea_trade_exp_val intra_ea_trade_exp_val
#> [1,] -3.909774e-20 1.130804e-05             -0.1173080              0.5300049
#> [2,] -3.037622e-20 1.125250e-05             -0.6082641              0.5613434
#> [3,] -2.429536e-20 1.130307e-05              0.2089724              0.6935698
#> [4,] -2.032653e-20 1.146023e-05             -0.3888678             -0.8632375
#> [5,] -1.812472e-20 1.172545e-05             -0.7122444              0.9965301
#> [6,] -1.749851e-20 1.210123e-05              0.1964240             -0.8188774
#>      extra_ea_trade_imp_val intra_ea_trade_imp_val      us_ip     us_urx
#> [1,]              0.7441252              0.3905597  0.4464705 -0.5972064
#> [2,]              2.2986685             -2.2884988 -0.6994935 -0.3270729
#> [3,]             -2.5212749              1.2639157 -1.5727228  1.6812686
#> [4,]              1.4758380             -0.1257012 -1.9704646  1.5982896
#> [5,]              0.1805347              0.1019460 -0.7847065 -0.6068660
#> [6,]              0.1648355              0.3371970 -0.5442462  0.4663372
#>          us_empl us_retail_sales us_ip_manuf_exp us_cons_exp   us_r3_m
#> [1,]  0.68805923   -3.879492e-85      1.09995435 -0.16821408  1.841994
#> [2,] -1.05923627    8.964934e-85      0.21804294 -0.96714338  3.912328
#> [3,] -0.32778663   -3.271557e-84     -2.65396143  0.31640642 -2.906317
#> [4,]  0.87144459    1.245807e-83     -2.95870846 -0.02023559 -7.023246
#> [5,]  0.01071255   -4.758257e-83      0.16277231  0.99325909 -1.871839
#> [6,]  0.64722814    1.817750e-82      0.06523476 -0.58291072  2.277664
#>      us_r10_year          m3        loans      ir_long       ir_short
#> [1,]   4.8907708  0.07194975 6.259551e-78 3.374236e-65 -1.710178e-112
#> [2,]   0.1963391  0.94025272 1.492565e-78 8.668876e-66 -3.188068e-113
#> [3,]  -3.0176803 -0.63796022 3.566169e-79 4.511839e-66 -5.943110e-114
#> [4,]  -2.1176326 -0.08400808 8.822942e-80 1.124106e-65 -1.107899e-114
#> [5,]  -0.2053264  0.37246930 3.448220e-80 4.509294e-65 -2.065324e-115
#> [6,]   1.8400670  1.31812770 6.462160e-80 1.877460e-64 -3.850596e-116
#>          ir_1_year     ir_2_year     ir_5_year          eer       eer_cpi
#> [1,]  1.220272e-98 8.549943e-151 3.681757e-156 0.0004557372  6.827481e-09
#> [2,]  2.745420e-99 1.796414e-151 7.405038e-157 0.0002974932 -4.009777e-10
#> [3,] 6.176765e-100 3.808120e-152 1.618525e-157 0.0002342173  4.175719e-11
#> [4,] 1.389675e-100 9.676741e-153 9.959707e-158 0.0002457102 -3.143740e-10
#> [5,] 3.126552e-101 1.002608e-152 3.551094e-157 0.0003356407  5.343853e-09
#> [6,] 7.034274e-102 4.016726e-152 1.743610e-156 0.0005327172 -9.123236e-08
#>           eer_ppi     exr_usd    exr_gbp    rxr_yen        euro50       euro325
#> [1,] 2.978315e-69  0.73454196 -1.0565458  2.4954660 -2.529723e-60 -6.362085e-59
#> [2,] 6.337630e-69  0.37396618 -0.2080339 -0.4297539 -3.060252e-61 -8.108712e-60
#> [3,] 2.205470e-68 -0.20597328 -1.6671123  0.4294326 -3.702041e-62 -1.033485e-60
#> [4,] 8.077629e-68  0.53680789 -1.5342581 -3.8573492 -4.478426e-63 -1.317215e-61
#> [5,] 2.970038e-67  0.07440597 -0.7810578 -2.2184222 -5.417632e-64 -1.678839e-62
#> [6,] 1.092360e-66  0.24086153 -0.4188286  1.8411043 -6.553806e-65 -2.139743e-63
#>            sp500       dow_j  raw_mat_en raw_mat_oil raw_mat_gold
#> [1,] 0.431919190 -0.17589384  0.81784534 -0.34892116  -0.06331777
#> [2,] 0.010685447  0.10735651 -0.25481550 -0.01625208  -3.53588829
#> [3,] 0.109158959 -0.11497711 -0.22658581  0.06371555  -2.19498973
#> [4,] 0.007235308  0.03244591  0.08718499  0.16768271  -1.29570277
#> [5,] 0.166750379 -0.18766044 -0.69483895  0.06920533   3.06659121
#> [6,] 0.015507310  0.02006015  0.31392001 -0.07390485   1.48255729
#>      raw_mat_oil_fwd     raw_mat         gdp   priv_cons     invest      export
#> [1,]    1.322800e-28 -0.40362743  0.04292267 -0.07682953 0.04696923 -0.13995419
#> [2,]    4.144101e-29 -0.12983426  0.05129271 -0.05790409 0.06563161 -0.08840784
#> [3,]    1.299791e-29  0.37714437  0.09987802 -0.15427628 0.11585092 -0.22419183
#> [4,]    4.125192e-30  0.16670567  0.04691497 -0.01023578 0.08187876 -0.10553460
#> [5,]    1.463590e-30 -0.06393378  0.01649269  0.05402172 0.06276925 -0.06993188
#> [6,]    1.005651e-30 -0.08081010 -0.02636112  0.20745746 0.04525212  0.03071138
#>           import       empl    prductivity      capacity       gdp_us
#> [1,] -0.01986328 0.25685434 -2.190174e-184 -2.715323e-43 -0.087264721
#> [2,] -0.01055342 0.18337707  1.464337e-185  1.253871e-42 -0.104537037
#> [3,] -0.04353313 0.13800700 -9.790471e-187 -5.812130e-42 -0.194121324
#> [4,]  0.04643295 0.09459807  6.545849e-188  2.694602e-41 -0.112607167
#> [5,]  0.09398922 0.05497276 -4.376515e-189 -1.249274e-40 -0.065167393
#> [6,]  0.21910457 0.02206330  2.926111e-190  5.791894e-40 -0.001692607

# Compare with model without AR(1) errors
dfm_large_ar1$loglik  # Log-likelihood path
#>  [1] -37569.02 -37273.48 -37113.12 -36999.05 -36919.58 -36867.57 -36833.99
#>  [8] -36810.86 -36793.18 -36778.36 -36765.15 -36752.93 -36741.45 -36730.58
#> [15] -36720.29 -36710.54 -36701.32 -36692.67 -36684.63 -36677.25 -36670.55
#> [22] -36664.51 -36659.09 -36654.22 -36649.84 -36645.88 -36642.28
# }
```
