# DFM
<!-- badges: start -->
[![R-CMD-check](https://github.com/SebKrantz/DFM/workflows/R-CMD-check/badge.svg)](https://github.com/SebKrantz/DFM/actions)
<!-- badges: end -->
<!-- **NOTE**: Package is under development. Feel free to contribute. -->

*DFM* provides efficient estimation of Dynamic Factor Models via the EM Algorithm. Estimation can be done in 3 different ways following:

- Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large approximate dynamic factor models based on Kalman filtering. *Journal of Econometrics, 164*(1), 188-205.

- Doz, C., Giannone, D., & Reichlin, L. (2012). A quasi-maximum likelihood approach for large, approximate dynamic factor models. *Review of economics and statistics, 94*(4), 1014-1024.

- Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor models on datasets with arbitrary pattern of missing data. *Journal of Applied Econometrics, 29*(1), 133-160.

All 3 estimation methods support missing data, with various preprocessing options. The default is `em.method = "DGR"` following Doz, Giannone & Reichlin (2012). Using `em.method = "none"` generates Two-Step estimates following Doz, Giannone & Reichlin (2011). This is extremely efficient on bigger datasets. PCA and Two-Step estimates are also reported in EM based methods. Choosing `em.method = "BM"` adjusts the EM iterations suitably to allow for arbitrary patterns of missing data following Banbura & Modugno (2014), and is the preferred (though least efficient) method to deal with missing data (particularly if there are ragged-edges or series at different frequencies). 

### Comparison with Other R Packages

The implementation is based on efficient C++ code, making *DFM* orders of magnitude faster than packages such as [*MARSS*](<https://CRAN.R-project.org/package=MARSS>) that can be used to fit dynamic factor models, or packages like [*nowcasting*](<https://github.com/nmecsys/nowcasting>) and [*nowcastDFM*](<https://github.com/dhopp1/nowcastDFM>), which fit dynamic factor models specific to mixed-frequency nowcasting applications. The latter two packages additionally support blocking of variables into different groups for which factors are to be estimated, and EM adjustments for variables at different frequencies, whereas *DFM* with `em.method = "BM"` does allow mixed-frequency data but performs no adjustments for the frequency of the data. *DFM* currently also does not allow residual autocorrelation in the estimation, although the addition of this feature is planned. 

The package is intended to provide a robust and computationally efficient baseline implementation of Dynamic Factor Models for R, allowing straightforward application to various contexts such as time series dimensionality reduction and multivariate forecasting. The package is not meant to fit very general forms of the state space model such as provided by [*MARSS*](<https://CRAN.R-project.org/package=MARSS>), or advanced specifications of Dynamic Factor Models tailored to mixed-frequency nowcasting applications such as [*nowcasting*](<https://github.com/nmecsys/nowcasting>) and [*nowcastDFM*](<https://github.com/dhopp1/nowcastDFM>). Such software could however benefit from the functions and methods provided in *DFM*, most notably *DFM* exports stationary Kalman Filters and Smoothers also used in nowcasting applications, that are noticeably faster than the more general implementations provided by the [*FKF*](<https://CRAN.R-project.org/package=FKF>) package. 

<!-- Estimation with *DFM* also requires stationary data of a single frequency, and assumes time-invariant system matrices and classical assumptions (i.e. the 'exact factor model', assuming away residual autocorrelation in the observation equation). -->

### Usage Example 
```
library(DFM)

# Fit DFM with 6 factors and 3 lags in the transition equation
mod = DFM(diff(BM14_M), r = 6, p = 3, em.method = "BM")

# 'dfm' methods
summary(mod)
plot(mod)

# Forecasting 20 periods ahead
fc = predict(mod, n.ahead = 20)

# 'dfm_forecast' methods
print(fc)
plot(fc)
```
