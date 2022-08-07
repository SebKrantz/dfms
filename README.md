# DFM
Dynamic Factor Models for R

<!-- badges: start -->
[![R-CMD-check](https://github.com/SebKrantz/DFM/workflows/R-CMD-check/badge.svg)](https://github.com/SebKrantz/DFM/actions)
<!-- badges: end -->

**NOTE**: Package is under development. Feel free to contribute. 

*DFM* provides efficient estimation of Dynamic Factor Models via the EM Algorithm. Estimation can be done in 3 different ways following:

- Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large approximate dynamic factor models based on Kalman filtering. *Journal of Econometrics, 164*(1), 188-205.

- Doz, C., Giannone, D., & Reichlin, L. (2012). A quasi-maximum likelihood approach for large, approximate dynamic factor models. *Review of economics and statistics, 94*(4), 1014-1024.

- Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor models on datasets with arbitrary pattern of missing data. *Journal of Applied Econometrics, 29*(1), 133-160.

All 3 estimation methods supports missing data, with various preprocessing options. The default is `em.method = "DGR"` following Doz, Giannone & Reichlin (2012). Using `em.method = "none"` generates two-step estimates following Doz, Giannone & Reichlin (2011). This is extremely efficient on bigger datasets. PCA and two-tep estimates are also reported in EM based methods. Choosing `em.method = "BM"` adjusts the EM iterations suitably to allow for arbitrary patterns of missing data following Banbura & Modugno (2014), and is the preferred (though least efficient) method to deal with missing data. 

The implementation is based on efficient C++ code but less general compared to other packages such as *nowcasting* and *nowcastDFM* (i.e. it currently allows no blocking or mixed-frequency models). Estimation with *DFM* also requires stationary data of a single frequency, and assumes time-invariant system matrices and classical assumptions (i.e. the 'exact factor model', assuming away residual autocorrelation in the observation equation).

