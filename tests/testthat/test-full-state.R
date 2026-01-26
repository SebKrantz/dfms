library(collapse)
library(xts)

test_that("use.full.state preserves output shapes", {
  set.seed(2024)
  X <- matrix(rnorm(120), nrow = 24)
  colnames(X) <- paste0("v", seq_len(ncol(X)))
  X[20, 1] <- NA
  X[21, 2] <- NA

  mod <- DFM(X, r = 1, p = 1, em.method = "none")

  fit_full <- fitted(mod, use.full.state = TRUE)
  fit_compact <- fitted(mod, use.full.state = FALSE)
  expect_equal(dim(fit_full), dim(fit_compact))
  expect_equal(colnames(fit_full), colnames(fit_compact))

  res_full <- residuals(mod, use.full.state = TRUE)
  res_compact <- residuals(mod, use.full.state = FALSE)
  expect_equal(dim(res_full), dim(res_compact))
  expect_equal(colnames(res_full), colnames(res_compact))

  fc_full <- predict(mod, h = 3, use.full.state = TRUE)
  fc_compact <- predict(mod, h = 3, use.full.state = FALSE)
  expect_equal(dim(fc_full$X_fcst), dim(fc_compact$X_fcst))
  expect_equal(dim(fc_full$F_fcst), dim(fc_compact$F_fcst))
})

test_that("use.full.state shapes match with idio.ar1", {
  set.seed(2025)
  X <- matrix(rnorm(120), nrow = 24)
  colnames(X) <- paste0("v", seq_len(ncol(X)))
  X[20, 1] <- NA
  X[21, 2] <- NA

  mod <- DFM(X, r = 1, p = 1, idio.ar1 = TRUE, em.method = "none")

  fit_full <- fitted(mod, use.full.state = TRUE)
  fit_compact <- fitted(mod, use.full.state = FALSE)
  expect_equal(dim(fit_full), dim(fit_compact))

  res_full <- residuals(mod, use.full.state = TRUE)
  res_compact <- residuals(mod, use.full.state = FALSE)
  expect_equal(dim(res_full), dim(res_compact))

  fc_full <- predict(mod, h = 3, use.full.state = TRUE)
  fc_compact <- predict(mod, h = 3, use.full.state = FALSE)
  expect_equal(dim(fc_full$X_fcst), dim(fc_compact$X_fcst))
  expect_equal(dim(fc_full$F_fcst), dim(fc_compact$F_fcst))
})

compare_full_state <- function(mod, h = 3, tol = 1e-6, expect_equal_values = TRUE, exclude.vars = NULL) {
  fit_full <- fitted(mod, use.full.state = TRUE)
  fit_compact <- fitted(mod, use.full.state = FALSE)
  expect_equal(dim(fit_full), dim(fit_compact))
  if(expect_equal_values) {
    keep <- if(length(exclude.vars)) setdiff(colnames(fit_full), exclude.vars) else colnames(fit_full)
    expect_equal(fit_full[, keep, drop = FALSE], fit_compact[, keep, drop = FALSE], tolerance = tol)
  }

  res_full <- residuals(mod, use.full.state = TRUE)
  res_compact <- residuals(mod, use.full.state = FALSE)
  # fmedian(abs(res_full)) / fmedian(abs(res_compact))
  expect_equal(dim(res_full), dim(res_compact))
  if(expect_equal_values) {
    keep <- if(length(exclude.vars)) setdiff(colnames(res_full), exclude.vars) else colnames(res_full)
    expect_equal(res_full[, keep, drop = FALSE], res_compact[, keep, drop = FALSE], tolerance = tol)
  }

  fc_full <- predict(mod, h = h, use.full.state = TRUE)
  fc_compact <- predict(mod, h = h, use.full.state = FALSE)
  expect_equal(dim(fc_full$X_fcst), dim(fc_compact$X_fcst))
  expect_equal(dim(fc_full$F_fcst), dim(fc_compact$F_fcst))
  if(expect_equal_values) {
    keep <- if(length(exclude.vars)) setdiff(colnames(fc_full$X_fcst), exclude.vars) else colnames(fc_full$X_fcst)
    expect_equal(fc_full$X_fcst[, keep, drop = FALSE], fc_compact$X_fcst[, keep, drop = FALSE], tolerance = tol)
    expect_equal(fc_full$F_fcst, fc_compact$F_fcst, tolerance = tol)
  }
}

test_that("use.full.state is numerically close for BM14 small monthly", {
  skip_on_cran()
  X <- qM(BM14_M[, BM14_Models$small[BM14_Models$freq == "M"]])
  mod <- DFM(X, r = 2, p = 2, em.method = "none")
  compare_full_state(mod, h = 3, tol = 1e-6, expect_equal_values = TRUE)
})

test_that("use.full.state preserves shapes for BM14 small monthly with idio.ar1", {
  skip_on_cran()
  X <- qM(BM14_M[, BM14_Models$small[BM14_Models$freq == "M"]])
  mod <- DFM(X, r = 2, p = 2, idio.ar1 = TRUE, em.method = "none")
  compare_full_state(mod, h = 3, tol = 1e-6, expect_equal_values = FALSE)
})

test_that("use.full.state is numerically close for BM14 small MQ", {
  skip_on_cran()
  BM14 <- merge.xts(BM14_M, BM14_Q)
  BM14[, BM14_Models$log_trans] <- log(BM14[, BM14_Models$log_trans])
  BM14[, BM14_Models$freq == "M"] <- fdiff(BM14[, BM14_Models$freq == "M"])
  BM14[, BM14_Models$freq == "Q"] <- fdiff(BM14[, BM14_Models$freq == "Q"], 3)
  X <- qM(BM14[, BM14_Models$small])
  colnames(X) <- BM14_Models$series[BM14_Models$small]
  quarterly.vars <- intersect(colnames(X), colnames(BM14_Q))
  mod <- DFM(X, r = 2, p = 2, quarterly.vars = quarterly.vars)
  compare_full_state(mod, h = 3, tol = 1e-6, expect_equal_values = TRUE, exclude.vars = quarterly.vars)
})

test_that("use.full.state preserves shapes for BM14 small MQ with idio.ar1", {
  skip_on_cran()
  BM14 <- merge.xts(BM14_M, BM14_Q)
  BM14[, BM14_Models$log_trans] <- log(BM14[, BM14_Models$log_trans])
  BM14[, BM14_Models$freq == "M"] <- fdiff(BM14[, BM14_Models$freq == "M"])
  BM14[, BM14_Models$freq == "Q"] <- fdiff(BM14[, BM14_Models$freq == "Q"], 3)
  X <- qM(BM14[, BM14_Models$small])
  colnames(X) <- BM14_Models$series[BM14_Models$small]
  quarterly.vars <- BM14_Models$series[BM14_Models$small & BM14_Models$freq == "Q"]
  mod <- DFM(X, r = 2, p = 2, quarterly.vars = quarterly.vars, idio.ar1 = TRUE)
  compare_full_state(mod, h = 3, tol = 1e-6, expect_equal_values = FALSE)
})
