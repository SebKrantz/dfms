library(collapse)
library(xts)

test_that("news returns expected shapes", {
  set.seed(123)
  X <- matrix(rnorm(120), nrow = 24)
  colnames(X) <- paste0("v", seq_len(ncol(X)))
  X_old <- X
  X_new <- X
  X_old[20, 1] <- NA
  X_new[20, 1] <- NA
  X_old[18, 2] <- NA
  X_new[18, 2] <- X[18, 2]

  dfm_old <- DFM(X_old, r = 1, p = 1, em.method = "none")
  dfm_new <- DFM(X_new, r = 1, p = 1, em.method = "none")
  groups <- rep(c("g1", "g2"), length.out = ncol(X))

  res <- news(dfm_old, dfm_new, t.fcst = 20, target.vars = 1, groups = groups)
  expect_s3_class(res, "dfm.news")

  expect_true(is.list(res))
  expect_length(res$singlenews, ncol(X))
  expect_length(res$groupnews, length(unique(groups)))
  expect_equal(names(res$singlenews), colnames(X))
  expect_equal(names(res$groupnews), unique(groups))
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

test_that("news works with MQ small model for GDP nowcast", {
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

  # Remove some recent monthly observations from X_old
  X_old[356, "ip_tot_cstr"] <- NA
  X_old[354, "orders"] <- NA
  X_old[355, "orders"] <- NA
  X_old[355, "extra_ea_trade_exp_val"] <- NA
  X_old[356, "new_cars"] <- NA
  X_old[357, "new_cars"] <- NA

  # Fit models with mixed-frequency
  dfm_old <- DFM(X_old, r = 2, p = 2, quarterly.vars = quarterly.vars)
  dfm_new <- DFM(X_new, r = 2, p = 2, quarterly.vars = quarterly.vars)

  # Nowcast GDP (variable 11) at row 354 (last observed quarter)
  gdp_idx <- which(colnames(X_new) == "gdp")
  res <- news(dfm_old, dfm_new, t.fcst = 354, target.vars = gdp_idx)

  expect_s3_class(res, "dfm.news")
  expect_true(is.numeric(res$y_old))
  expect_true(is.numeric(res$y_new))

  # Key check: revision should equal sum of singlenews
  revision <- res$y_new - res$y_old
  sum_news <- sum(res$singlenews)
  expect_equal(revision, sum_news, tolerance = 1e-10)

  # GDP should not be in the released variables (we're nowcasting it)
  expect_true(res$singlenews["gdp"] == 0 || is.na(res$singlenews["gdp"]))

  # Check that monthly variables contributed news
  monthly_news <- res$singlenews[c("ip_tot_cstr", "orders", "new_cars", "extra_ea_trade_exp_val")]
  expect_true(any(monthly_news != 0))
})

test_that("news works with MQ medium model for GDP nowcast", {
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

  # Remove some recent monthly observations from X_old
  X_old[356, "ip_tot_cstr"] <- NA
  X_old[354, "orders"] <- NA
  X_old[355, "orders"] <- NA
  X_old[355, "extra_ea_trade_exp_val"] <- NA
  X_old[356, "new_cars"] <- NA
  X_old[357, "new_cars"] <- NA
  X_old[356, "ip_capital"] <- NA
  X_old[356, "us_ip"] <- NA

  # Fit models with mixed-frequency
  dfm_old <- DFM(X_old, r = 3, p = 2, quarterly.vars = quarterly.vars)
  dfm_new <- DFM(X_new, r = 3, p = 2, quarterly.vars = quarterly.vars)

  # Nowcast GDP at row 354
  gdp_idx <- which(colnames(X_new) == "gdp")
  res <- news(dfm_old, dfm_new, t.fcst = 354, target.vars = gdp_idx)

  expect_s3_class(res, "dfm.news")
  expect_true(is.numeric(res$y_old))
  expect_true(is.numeric(res$y_new))

  # Key check: revision should equal sum of singlenews
  revision <- res$y_new - res$y_old
  sum_news <- sum(res$singlenews)
  expect_equal(revision, sum_news, tolerance = 1e-10)

  # GDP should not be in the released variables
  expect_true(res$singlenews["gdp"] == 0 || is.na(res$singlenews["gdp"]))

  # Check gains exist for released variables
  expect_true(length(res$gain) > 0)
})
