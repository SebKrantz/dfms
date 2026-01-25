library(collapse)
library(xts)

test_that("news returns expected shapes", {
  set.seed(123)
  X <- matrix(rnorm(120), nrow = 24)
  colnames(X) <- paste0("v", seq_len(ncol(X)))
  X_old <- X
  X_new <- X
  X_old[20, 1] <- NA
  X_new[20, 1] <- 0.9
  X_old[18, 2] <- NA
  X_new[18, 2] <- X[18, 2]

  dfm_old <- DFM(X_old, r = 1, p = 1, em.method = "none")
  dfm_new <- DFM(X_new, r = 1, p = 1, em.method = "none")

  res <- news(dfm_old, dfm_new, t.fcst = 20, target.vars = 1)
  expect_s3_class(res, "dfm.news")

  expect_true(is.list(res))
  expect_length(res$singlenews, ncol(X))
  expect_equal(names(res$singlenews), colnames(X))
  expect_true(is.numeric(res$y_old))
  expect_true(is.numeric(res$y_new))
})

test_that("news uses simple-case shortcut", {
  set.seed(42)
  X <- matrix(rnorm(80), nrow = 20)
  colnames(X) <- paste0("v", seq_len(ncol(X)))
  X_old <- X
  X_new <- X
  X_old[15, 2] <- NA
  X_new[15, 2] <- X[15, 2]

  dfm_old <- DFM(X_old, r = 1, p = 1, em.method = "none")
  dfm_new <- DFM(X_new, r = 1, p = 1, em.method = "none")

  res <- news(dfm_old, dfm_new, t.fcst = 10, target.vars = 1)

  res_data <- news(dfm_old, X_new, t.fcst = 10, target.vars = 1)
  expect_s3_class(res_data, "dfm.news")
  res_all <- news(dfm_old, X_new, t.fcst = 10)
  expect_s3_class(res_all, "dfm.news_list")
  expect_length(res_data$gain, 0)
  expect_true(is.null(res_data$gain_scaled) || length(res_data$gain_scaled) == ncol(X))

  expect_length(res$gain, 0)
  expect_true(is.null(names(res$gain)) || length(names(res$gain)) == 0L)
})

test_that("news works with MQ small model for monthly target", {
  skip_on_cran()

  # Construct BM14 database
  BM14 <- merge(BM14_M, BM14_Q)
  BM14[, BM14_Models$log_trans] <- log(BM14[, BM14_Models$log_trans])
  BM14[, BM14_Models$freq == "M"] <- fdiff(BM14[, BM14_Models$freq == "M"])
  BM14[, BM14_Models$freq == "Q"] <- fdiff(BM14[, BM14_Models$freq == "Q"], 3)

  # Small model data
  X_small <- qM(BM14[, BM14_Models$small])
  colnames(X_small) <- BM14_Models$series[BM14_Models$small]
  quarterly.vars <- BM14_Models$series[BM14_Models$small & BM14_Models$freq == "Q"]

  # Create vintages: X_old has fewer monthly observations (simulating older vintage)
  X_old <- X_small
  X_new <- X_small

  # Create releases from observed variables (new_cars, pms_pmi are observed near end)
  X_old[355, "new_cars"] <- NA
  X_old[356, "new_cars"] <- NA
  X_old[354, "pms_pmi"] <- NA
  X_old[355, "pms_pmi"] <- NA

  # Fit models with mixed-frequency
  dfm_old <- DFM(X_old, r = 2, p = 2, quarterly.vars = quarterly.vars)
  dfm_new <- DFM(X_new, r = 2, p = 2, quarterly.vars = quarterly.vars)

  # Target: 'orders' at t=355 (X_imp coordinates) - naturally missing at end of sample
  # This is a proper nowcast scenario where target is NOT observed
  res_m <- news(dfm_old, dfm_new, t.fcst = 355, target.vars = "orders")

  expect_s3_class(res_m, "dfm.news")
  expect_true(is.numeric(res_m$y_old))
  expect_true(is.numeric(res_m$y_new))

  # For monthly targets, revision should equal sum of singlenews exactly
  revision_m <- unname(res_m$y_new - res_m$y_old)
  sum_news_m <- sum(res_m$singlenews)
  expect_equal(revision_m, sum_news_m, tolerance = 1e-10)

  # Check that released variables contributed news
  expect_true(all(res_m$singlenews[c("new_cars", "pms_pmi")] != 0))

  # Check gains exist for released variables
  expect_true(length(res_m$gain) > 0)

  # # Check that gain produces results.
  # all.equal(res_m$singlenews[names(res_m$gain)], unclass((res_m$forecasts[names(res_m$gain)] - res_m$actual[names(res_m$gain)]) * res_m$gain))

})

test_that("news works with MQ medium model for monthly target", {
  skip_on_cran()

  # Construct BM14 database
  BM14 <- merge(BM14_M, BM14_Q)
  BM14[, BM14_Models$log_trans] <- log(BM14[, BM14_Models$log_trans])
  BM14[, BM14_Models$freq == "M"] <- fdiff(BM14[, BM14_Models$freq == "M"])
  BM14[, BM14_Models$freq == "Q"] <- fdiff(BM14[, BM14_Models$freq == "Q"], 3)

  # Medium model data
  X_medium <- qM(BM14[, BM14_Models$medium])
  colnames(X_medium) <- BM14_Models$series[BM14_Models$medium]
  quarterly.vars <- BM14_Models$series[BM14_Models$medium & BM14_Models$freq == "Q"]

  # Create vintages
  X_old <- X_medium
  X_new <- X_medium

  # Create releases from observed variables
  X_old[355, "new_cars"] <- NA
  X_old[356, "new_cars"] <- NA
  X_old[354, "pms_pmi"] <- NA
  X_old[355, "pms_pmi"] <- NA
  X_old[355, "ip_capital"] <- NA
  X_old[356, "us_ip"] <- NA

  # Fit models with mixed-frequency
  dfm_old <- DFM(X_old, r = 3, p = 2, quarterly.vars = quarterly.vars)
  dfm_new <- DFM(X_new, r = 3, p = 2, quarterly.vars = quarterly.vars)

  # Target: 'orders' at t=355 - naturally missing at end of sample
  res_m <- news(dfm_old, dfm_new, t.fcst = 355, target.vars = "orders")

  expect_s3_class(res_m, "dfm.news")
  expect_true(is.numeric(res_m$y_old))
  expect_true(is.numeric(res_m$y_new))

  # For monthly targets, revision should equal sum of singlenews exactly
  revision_m <- unname(res_m$y_new - res_m$y_old)
  sum_news_m <- sum(res_m$singlenews)
  expect_equal(revision_m, sum_news_m, tolerance = 1e-10)

  # Check gains exist for released variables
  expect_true(length(res_m$gain) > 0)
})

test_that("news works with idio.ar1 model", {
  set.seed(202)
  X <- matrix(rnorm(120), nrow = 24)
  colnames(X) <- paste0("v", seq_len(ncol(X)))
  X_old <- X
  X_new <- X
  X_old[20, 1] <- NA
  X_new[20, 1] <- 0.7

  dfm_old <- DFM(X_old, r = 1, p = 2, idio.ar1 = TRUE, em.method = "none")
  dfm_new <- DFM(X_new, r = 1, p = 2, idio.ar1 = TRUE, em.method = "none")

  res <- news(dfm_old, dfm_new, t.fcst = 20, target.vars = 1)
  expect_s3_class(res, "dfm.news")

  revision <- unname(res$y_new - res$y_old)
  sum_news <- sum(res$singlenews)
  expect_equal(revision, sum_news, tolerance = 1e-8)
})

test_that("news works with MQ + idio.ar1 model", {
  skip_on_cran()

  BM14 <- merge(BM14_M, BM14_Q)
  BM14[, BM14_Models$log_trans] <- log(BM14[, BM14_Models$log_trans])
  BM14[, BM14_Models$freq == "M"] <- fdiff(BM14[, BM14_Models$freq == "M"])
  BM14[, BM14_Models$freq == "Q"] <- fdiff(BM14[, BM14_Models$freq == "Q"], 3)

  X_small <- qM(BM14[, BM14_Models$small])
  colnames(X_small) <- BM14_Models$series[BM14_Models$small]
  quarterly.vars <- BM14_Models$series[BM14_Models$small & BM14_Models$freq == "Q"]

  X_old <- X_small
  X_new <- X_small
  X_old[355, "new_cars"] <- NA
  X_old[356, "new_cars"] <- NA
  X_old[354, "pms_pmi"] <- NA
  X_old[355, "pms_pmi"] <- NA

  dfm_old <- DFM(X_old, r = 2, p = 2, quarterly.vars = quarterly.vars, idio.ar1 = TRUE)
  dfm_new <- DFM(X_new, r = 2, p = 2, quarterly.vars = quarterly.vars, idio.ar1 = TRUE)

  res_m <- news(dfm_old, dfm_new, t.fcst = 355, target.vars = "orders")

  expect_s3_class(res_m, "dfm.news")
  revision_m <- unname(res_m$y_new - res_m$y_old)
  sum_news_m <- sum(res_m$singlenews)
  expect_equal(revision_m, sum_news_m, tolerance = 1e-10)
})
