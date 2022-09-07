library(magrittr)
library(xts)

# BM14 Replication Data. Constructing the database:
BM14 = merge(BM14_M, BM14_Q)
BM14[, BM14_Models$log_trans] %<>% log()
BM14[, BM14_Models$freq == "M"] %<>% diff()
BM14[, BM14_Models$freq == "Q"] %<>% diff(3)

### Small Model ---------------------------------------

# Estimating the model with 2 factors and 3 lags
dfm_small = DFM(BM14[, BM14_Models$small], 2, 3, em.method = "BM")

expect_equal(dfm_small$A[1, 1], 1.21068973, tolerance = 1e-5)
expect_equal(dfm_small$A[2, 6], 0.391438990, tolerance = 1e-5)

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

