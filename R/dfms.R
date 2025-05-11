#' Dynamic Factor Models
#'
#' *dfms* provides efficient estimation of Dynamic Factor Models via the EM Algorithm.
#'
#' Estimation can be done in 3 different ways following:
#'
#'  - Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large approximate dynamic factor models based on Kalman filtering. *Journal of Econometrics, 164*(1), 188-205. <doi:10.1016/j.jeconom.2011.02.012>
#'
#'  - Doz, C., Giannone, D., & Reichlin, L. (2012). A quasi-maximum likelihood approach for large, approximate dynamic factor models. *Review of Economics and Statistics, 94*(4), 1014-1024. <doi:10.1162/REST_a_00225>
#'
#'  - Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor models on datasets with arbitrary pattern of missing data. *Journal of Applied Econometrics, 29*(1), 133-160. <doi:10.1002/jae.2306>
#'
#'  The default is `em.method = "auto"`, which chooses `"BM"` following Banbura & Modugno (2014) with missing data or mixed frequency, and `"DGR"` following Doz, Giannone & Reichlin (2012) otherwise. Using `em.method = "none"` generates Two-Step estimates following Doz, Giannone & Reichlin (2011). This is extremely efficient on bigger datasets. PCA and Two-Step estimates are also reported in EM-estimation. All methods support missing data, but `em.method = "DGR"` does not model them in EM iterations.
#'
#' @section Package Contents:
#'
#' **Functions to Specify/Estimate Model and Key Methods**
#'
#'  \code{\link[=ICr]{ICr()}} --- Information Criteria\cr
#'
#'  - \code{\link[=plot.ICr]{plot(<ICr>)}}\cr
#'  - \code{\link[=screeplot.ICr]{screeplot(<ICr>)}}\cr
#'
#'  \code{\link[=DFM]{DFM()}} --- Estimate the Model\cr
#'
#'  - \code{\link[=summary.dfm]{summary(<dfm>)}}\cr
#'  - \code{\link[=plot.dfm]{plot(<dfm>)}}\cr
#'  - \code{\link[=as.data.frame.dfm]{as.data.frame(<dfm>)}}\cr
#'  - \code{\link[=residuals.dfm]{residuals(<dfm>)}}\cr
#'  - \code{\link[=fitted.dfm]{fitted(<dfm>)}}
#'
#'  \code{\link[=predict.dfm]{predict(<dfm>)}} --- Generate Forecasts\cr
#'
#'  - \code{\link[=plot.dfm_forecast]{plot(<dfm_forecast>)}}\cr
#'  - \code{\link[=as.data.frame.dfm_forecast]{as.data.frame(<dfm_forecast>)}}\cr
#'
#' **Auxiliary Functions**
#'
#'  \code{\link[=.VAR]{.VAR()}} --- Estimate Vector Autoregression\cr
#'  \code{\link[=SKF]{SKF()}} --- Stationary Kalman Filter\cr
#'  \code{\link[=FIS]{FIS()}} --- Fixed Interval Smoother\cr
#'  \code{\link[=SKFS]{SKFS()}} --- Stationary Kalman Filter + Smoother\cr
#'  \code{\link[=tsnarmimp]{tsnarmimp()}} --- Remove and Impute Missing Values in a Multivariate Time Series\cr
#'  \code{\link[=ainv]{ainv()}} --- Rcpp Armadillo's Inverse Function\cr
#'  \code{\link[=apinv]{apinv()}} --- Rcpp Armadillo's Pseudo-Inverse Function\cr
#'
#' **Data**
#'
#'  \code{\link{BM14_M}} --- Monthly Series by Banbura and Modugno (2014)\cr
#'  \code{\link{BM14_Q}} --- Quarterly Series by Banbura and Modugno (2014)\cr
#'  \code{\link{BM14_Models}} --- Series Metadata + Small/Medium/Large Model Specifications\cr
#'
#' @docType package
#' @name dfms-package
#' @aliases dfms
#'
NULL
