library(collapse)
library(magrittr)
library(xts)

#' @srrstats {G5.0} *Where applicable or practicable, tests should use standard data sets with known properties (for example, the [NIST Standard Reference Datasets](https://www.itl.nist.gov/div898/strd/), or data sets provided by other widely-used R packages).*
#' @srrstats {G5.1} *Data sets created within, and used to test, a package should be exported (or otherwise made generally available) so that users can confirm tests and run examples.*
#' @srrstats {G5.2} *Appropriate error and warning behaviour of all functions should be explicitly demonstrated through tests. In particular,*
#' @srrstats {G5.2a} *Every message produced within R code by `stop()`, `warning()`, `message()`, or equivalent should be unique*
#' @srrstats {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstats {G5.3} *For functions which are expected to return objects containing no missing (`NA`) or undefined (`NaN`, `Inf`) values, the absence of any such values in return objects should be explicitly tested.*
#' @srrstats {G5.4} **Correctness tests** *to test that statistical algorithms produce expected results to some fixed test data sets (potentially through comparisons using binding frameworks such as [RStata](https://github.com/lbraglia/RStata)).*
#' -> I have translated the authors original matlab code into R and run tests with that code (see misc/ directory in the repo). These comparisons yielded that my implementation is equivalent to the original Matlab code.
#'  testing can be improved upon in general, an idea would be to run tests against the r-translations of those Matlab codes in the testthat framework, but I have not come round to do that yet and I would first like to gather some more substantive feedback on the software.


# BM14 Replication Data. Constructing the database:
BM14 = merge(BM14_M, BM14_Q)
BM14[, BM14_Models$log_trans] %<>% log()
BM14[, BM14_Models$freq == "M"] %<>% diff()
BM14[, BM14_Models$freq == "Q"] %<>% diff(3)

### Missing value removal
for(narm in c("LE", "all")) {
  for(m in c("median", "rnorm", "median.ma", "median.ma.spline"))
    expect_false(anyNA(tsnarmimp(BM14, na.rm.method = narm, na.impute = m)))
}


### Small Model ---------------------------------------

ics = ICr(BM14[, BM14_Models$small], max.r = 5)
expect_output(print(ics))

# Estimating the model with 2 factors and 3 lags
dfm_small = DFM(BM14[, BM14_Models$small], 2, 3, em.method = "BM", check.increased = TRUE)

expect_equal(dfm_small$A[1, 1], 1.21068973, tolerance = 1e-5)
expect_equal(dfm_small$A[2, 6], 0.391438990, tolerance = 1e-5)

# Rudimentry testing various methods
expect_output(print(dfm_small))
expect_visible(summary(dfm_small))
expect_true(is.matrix(fitted(dfm_small)))
expect_true(is.matrix(residuals(dfm_small)))
expect_true(ncol(as.data.frame(dfm_small)) == 4L)
expect_true(ncol(as.data.frame(dfm_small, pivot = "wide.method")) == 5L)
expect_true(ncol(as.data.frame(dfm_small, pivot = "wide.factor")) == 4L)
expect_true(ncol(as.data.frame(dfm_small, pivot = "wide")) == 7L)
expect_true(ncol(as.data.frame(dfm_small, method = "qml", pivot = "wide")) == 3L)
expect_true(ncol(as.data.frame(dfm_small, pivot = "t.wide")) == 7L)
expect_true(is.list(predict(dfm_small)))
expect_true(is.list(predict(dfm_small, standardized = FALSE)))
expect_true(is.list(predict(dfm_small, resFUN = function(x, h) predict(ar(na_rm(x)), n.ahead = h)$pred)))

expect_output(print(predict(dfm_small)))
expect_true(ncol(as.data.frame(predict(dfm_small))) == 4L)
expect_true(ncol(as.data.frame(predict(dfm_small), pivot = "wide")) == 4L)
expect_true(ncol(as.data.frame(predict(dfm_small), use = "data", pivot = "wide")) == 16L)
expect_true(ncol(as.data.frame(predict(dfm_small), use = "both", pivot = "wide")) == 18L)

# Other missing value options
expect_visible(tsnarmimp(BM14[, BM14_Models$small], na.rm.method = "all"))
expect_visible(tsnarmimp(BM14[, BM14_Models$small], na.impute = "median.ma"))
expect_visible(tsnarmimp(BM14[, BM14_Models$small], na.impute = "median"))
expect_visible(tsnarmimp(BM14[, BM14_Models$small], na.impute = "rnorm"))

### Medium-Sized Model ---------------------------------

# Estimating the model with 3 factors and 3 lags
dfm_medium = DFM(BM14[, BM14_Models$medium], 3, 3, em.method = "BM")

expect_equal(dfm_medium$A[1, 1], 0.74619087, tolerance = 1e-5)
expect_equal(dfm_medium$A[3, 9], 0.15380781, tolerance = 1e-5)

### Large Model ---------------------------------

# Estimating the model with 6 factors and 3 lags
dfm_large = DFM(BM14, 6, 3, em.method = "BM")

expect_equal(dfm_large$A[1, 1], 0.48915420, tolerance = 1e-5)
expect_equal(dfm_large$A[6, 18], 0.10027110, tolerance = 1e-5)

# Twostep estimates
dfm_large_2s = DFM(BM14, 6, 3, em.method = "none")
expect_equal(dfm_large_2s$A[1, 1], 0.6778009, tolerance = 1e-5)
expect_equal(dfm_large_2s$A[6, 18], 0.1314988, tolerance = 1e-5)

# Now testing DGR 2012 method
Xcols = colSums(is.na(BM14_M[70:350, ])) == 0
X_cc = BM14_M[70:350, Xcols]

expect_identical(ICr(diff(X_cc), max.r = 10)$r.star, c(IC1 = 7L, IC2 = 4L, IC3 = 7L))
mod = DFM(diff(X_cc), 4, 3)

expect_equal(mod$F_qml[10,], c(f1 = 0.4147068, f2 = 3.2364928, f3 = 1.6050372, f4 = -3.1349201), tolerance = 1e-5)
expect_equal(mod$A[1, 1], 0.5770928, tolerance = 1e-5)
expect_equal(mod$A[4, 12], -0.02980421, tolerance = 1e-5)

# BM should give almost the same...
mod_BM = DFM(diff(X_cc), 4, 3, em.method = "BM")

expect_equal(mod$F_twostep, mod_BM$F_twostep, tolerance = 1e-2)
expect_equal(mod$F_qml, mod_BM$F_qml, tolerance = 1e-2)
expect_equal(mod$A, mod_BM$A, tolerance = 1e-1)

