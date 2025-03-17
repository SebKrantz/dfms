# **dfms**: Dynamic Factor Models for R
<!-- badges: start -->
[![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/556_status.svg)](https://github.com/ropensci/software-review/issues/556)
[![R-CMD-check](https://github.com/SebKrantz/dfms/workflows/R-CMD-check/badge.svg)](https://github.com/SebKrantz/dfms/actions)
[![dfms status badge](https://sebkrantz.r-universe.dev/badges/dfms)](https://sebkrantz.r-universe.dev)
[![CRAN status](https://www.r-pkg.org/badges/version/dfms)](https://cran.r-project.org/package=dfms) 
[![cran checks](https://badges.cranchecks.info/worst/dfms.svg)](https://cran.r-project.org/web/checks/check_results_dfms.html)
![downloads per month](https://cranlogs.r-pkg.org/badges/dfms?color=blue)
![downloads](https://cranlogs.r-pkg.org/badges/grand-total/dfms?color=blue)
[![Codecov test coverage](https://codecov.io/gh/SebKrantz/dfms/branch/main/graph/badge.svg)](https://app.codecov.io/gh/SebKrantz/dfms?branch=main)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![dependencies](https://tinyverse.netlify.app/badge/dfms)](https://CRAN.R-project.org/package=dfms)
<!-- [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) -->
<!-- badges: end

<!-- **NOTE**: This package is under [rOpenSci Statistical Software Peer Review](https://stats-devguide.ropensci.org/). Peer review might result in changes to the API. -->
<!--
The package is fully functional though, and you are very welcome to install it using `remotes::install_github("SebKrantz/dfms")` and give feedback. -->

*dfms* provides efficient estimation of Dynamic Factor Models via the EM Algorithm. Estimation can be done in 3 different ways following:

- Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large approximate dynamic factor models based on Kalman filtering. *Journal of Econometrics, 164*(1), 188-205. <doi:10.1016/j.jeconom.2011.02.012> 

- Doz, C., Giannone, D., & Reichlin, L. (2012). A quasi-maximum likelihood approach for large, approximate dynamic factor models. *Review of Economics and Statistics, 94*(4), 1014-1024. <doi:10.1162/REST_a_00225>

- Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor models on datasets with arbitrary pattern of missing data. *Journal of Applied Econometrics, 29*(1), 133-160. <doi:10.1002/jae.2306>

The default is `em.method = "auto"`, which chooses `"BM"` following Banbura & Modugno (2014) with missing data or mixed frequency, and `"DGR"` following Doz, Giannone & Reichlin (2012) otherwise. Using `em.method = "none"` generates Two-Step estimates following Doz, Giannone & Reichlin (2011). This is extremely efficient on bigger datasets. PCA and Two-Step estimates are also reported in EM-estimation. All methods support missing data, but `em.method = "DGR"` does not model them in EM iterations.

### Comparison with Other R Packages

*dfms* is intended to provide a simple, numerically robust, and computationally efficient baseline implementation of (linear Gaussian) Dynamic Factor Models for R, allowing straightforward application to various contexts such as time series dimensionality reduction and forecasting. The implementation is based on efficient C++ code, making *dfms* orders of magnitude faster than packages that can be used to fit dynamic factor models such as [*MARSS*](<https://CRAN.R-project.org/package=MARSS>), or [*nowcasting*](<https://github.com/nmecsys/nowcasting>) and [*nowcastDFM*](<https://github.com/dhopp1/nowcastDFM>) geared to mixed-frequency nowcasting applications - supporting blocking of variables into different groups for which factors are to be estimated and evaluation of news content. For large-scale nowcasting models the [`DynamicFactorMQ`](https://www.statsmodels.org/dev/generated/statsmodels.tsa.statespace.dynamic_factor_mq.DynamicFactorMQ.html) class in the `statsmodels` Python library is probably the most robust implementation - see the [example](http://www.chadfulton.com/topics/statespace_large_dynamic_factor_models.html) by Chad Fulton.   <!-- , and EM adjustments for variables at different frequencies. *dfms* with `em.method = "BM"` does allow mixed-frequency data but performs no specific adjustments for the frequency of the data^[All series are weighted equally, and the prevalence of missing values in lower-frequency series downweights them. To remedy this lower frequency series could be included multiple times in the dataset e.g. include a quarterly series 3 times in a monthly dataset.]. *dfms* currently also does not allow residual autocorrelation in the estimation (i.e. it cannot estimate *approximate factor models*), but the addition of this feature is planned. -->
The *dfms* package is not intended to fit more general forms of the state space model like [*MARSS*](<https://CRAN.R-project.org/package=MARSS>).

<!-- or advanced specifications of Dynamic Factor Models tailored to mixed-frequency nowcasting applications such as [*nowcasting*](<https://github.com/nmecsys/nowcasting>) and [*nowcastDFM*](<https://github.com/dhopp1/nowcastDFM>). Such software could however benefit from the functions and methods provided in *dfms*, most notably *dfms* exports stationary Kalman Filters and Smoothers used in nowcasting applications, that are noticeably faster than the more general implementations provided by the [*FKF*](<https://CRAN.R-project.org/package=FKF>) package. -->

<!-- Estimation with *dfms* also requires stationary data of a single frequency, and assumes time-invariant system matrices and classical assumptions (i.e. the 'exact factor model', assuming away residual autocorrelation in the observation equation). -->

### Installation 

```r
# CRAN
install.packages("dfms")

# Development Version
install.packages('dfms', repos = c('https://sebkrantz.r-universe.dev', 'https://cloud.r-project.org'))

```
### Usage Example 
```r
library(dfms)

# Fit DFM with 6 factors and 3 lags in the transition equation
mod = DFM(diff(BM14_M), r = 6, p = 3) 

# 'dfm' methods
summary(mod)
plot(mod)
as.data.frame(mod)

# Forecasting 20 periods ahead
fc = predict(mod, h = 20)

# 'dfm_forecast' methods
print(fc)
plot(fc)
as.data.frame(fc)
```
