# Package index

- [`dfms-package`](https://sebkrantz.github.io/dfms/reference/dfms-package.md)
  [`dfms`](https://sebkrantz.github.io/dfms/reference/dfms-package.md) :
  Dynamic Factor Models

## Information Criteria

Choose the number of factors and the lag-order of the factor VAR.

- [`ICr()`](https://sebkrantz.github.io/dfms/reference/ICr.md)
  [`print(`*`<ICr>`*`)`](https://sebkrantz.github.io/dfms/reference/ICr.md)
  [`plot(`*`<ICr>`*`)`](https://sebkrantz.github.io/dfms/reference/ICr.md)
  [`screeplot(`*`<ICr>`*`)`](https://sebkrantz.github.io/dfms/reference/ICr.md)
  : Information Criteria to Determine the Number of Factors (r)

## Fit a Dynamic Factor Model

DFM estimation via the EM algorithm and PCA, and various methods inspect
the model and extract results.

- [`DFM()`](https://sebkrantz.github.io/dfms/reference/DFM.md) :
  Estimate a Dynamic Factor Model
- [`print(`*`<dfm>`*`)`](https://sebkrantz.github.io/dfms/reference/summary.dfm.md)
  [`coef(`*`<dfm>`*`)`](https://sebkrantz.github.io/dfms/reference/summary.dfm.md)
  [`logLik(`*`<dfm>`*`)`](https://sebkrantz.github.io/dfms/reference/summary.dfm.md)
  [`summary(`*`<dfm>`*`)`](https://sebkrantz.github.io/dfms/reference/summary.dfm.md)
  [`print(`*`<dfm_summary>`*`)`](https://sebkrantz.github.io/dfms/reference/summary.dfm.md)
  : DFM Summary Methods
- [`plot(`*`<dfm>`*`)`](https://sebkrantz.github.io/dfms/reference/plot.dfm.md)
  [`screeplot(`*`<dfm>`*`)`](https://sebkrantz.github.io/dfms/reference/plot.dfm.md)
  : Plot DFM
- [`as.data.frame(`*`<dfm>`*`)`](https://sebkrantz.github.io/dfms/reference/as.data.frame.dfm.md)
  : Extract Factor Estimates in a Data Frame
- [`residuals(`*`<dfm>`*`)`](https://sebkrantz.github.io/dfms/reference/residuals.dfm.md)
  [`fitted(`*`<dfm>`*`)`](https://sebkrantz.github.io/dfms/reference/residuals.dfm.md)
  : DFM Residuals and Fitted Values

## Forecasting

Forecast both the factors and the data, and methods to visualize
forecasts and extract results.

- [`predict(`*`<dfm>`*`)`](https://sebkrantz.github.io/dfms/reference/predict.dfm.md)
  [`print(`*`<dfm_forecast>`*`)`](https://sebkrantz.github.io/dfms/reference/predict.dfm.md)
  [`plot(`*`<dfm_forecast>`*`)`](https://sebkrantz.github.io/dfms/reference/predict.dfm.md)
  [`as.data.frame(`*`<dfm_forecast>`*`)`](https://sebkrantz.github.io/dfms/reference/predict.dfm.md)
  : DFM Forecasts

## Fast Stationary Kalman Filtering and Smoothing

Optimized Armadillo C++ implementations of the stationary Kalman Filter
and Smoother.

- [`SKF()`](https://sebkrantz.github.io/dfms/reference/SKF.md) : (Fast)
  Stationary Kalman Filter
- [`FIS()`](https://sebkrantz.github.io/dfms/reference/FIS.md) : (Fast)
  Fixed-Interval Smoother (Kalman Smoother)
- [`SKFS()`](https://sebkrantz.github.io/dfms/reference/SKFS.md) :
  (Fast) Stationary Kalman Filter and Smoother

## Helper Functions

Fast VAR, matrix inverses, imputation/removal of missing values in
multivariate time series, and convergence check for EM algorithm.

- [`.VAR()`](https://sebkrantz.github.io/dfms/reference/dot-VAR.md) :
  (Fast) Barebones Vector-Autoregression
- [`tsnarmimp()`](https://sebkrantz.github.io/dfms/reference/tsnarmimp.md)
  : Remove and Impute Missing Values in a Multivariate Time Series
- [`ainv()`](https://sebkrantz.github.io/dfms/reference/ainv.md)
  [`apinv()`](https://sebkrantz.github.io/dfms/reference/ainv.md) :
  Armadillo's Inverse Functions
- [`em_converged()`](https://sebkrantz.github.io/dfms/reference/em_converged.md)
  : Convergence Test for EM-Algorithm

## Data

Euro area macroeconomic data from Banbura and Modugno (2014), and 3 DFM
specifications considered in their paper.

- [`BM14_Models`](https://sebkrantz.github.io/dfms/reference/BM14_Models.md)
  [`BM14_M`](https://sebkrantz.github.io/dfms/reference/BM14_Models.md)
  [`BM14_Q`](https://sebkrantz.github.io/dfms/reference/BM14_Models.md)
  : Euro Area Macroeconomic Data from Banbura and Modugno 2014
