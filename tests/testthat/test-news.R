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
  expect_true(is.data.frame(res$news_df))
  expect_equal(res$news_df$series, colnames(X))
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
  expect_true(all(is.na(res_data$news_df$actual)))
  expect_true(all(is.na(res_data$news_df$forecast)))
  revision <- unname(res$y_new - res$y_old)
  sum_impact <- sum(res$news_df$impact)
  expect_equal(sum_impact, revision)
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
  # X_old[355, "new_cars"] <- NA
  # X_old[356, "new_cars"] <- NA
  X_old[357, "new_cars"] <- NA
  # X_old[354, "pms_pmi"] <- NA
  # X_old[356, "pms_pmi"] <- NA
  X_old[357, "pms_pmi"] <- NA

  # Fit models with mixed-frequency
  dfm_old <- DFM(X_old, r = 2, p = 2, quarterly.vars = quarterly.vars)
  dfm_new <- DFM(X_new, r = 2, p = 2, quarterly.vars = quarterly.vars)

  # This is a proper nowcast scenario where target is NOT observed
  res_m <- news(dfm_old, dfm_new, t.fcst = 356, target.vars = "orders")

  expect_s3_class(res_m, "dfm.news")
  expect_true(is.numeric(res_m$y_old))
  expect_true(is.numeric(res_m$y_new))

  # For monthly targets, revision should equal sum of impacts
  revision_m <- unname(res_m$y_new - res_m$y_old)
  sum_impact <- sum(res_m$news_df$impact)
  expect_equal(revision_m, sum_impact, tolerance = 1e-10)

  # Standardized scale should match more tightly
  res_m_std <- news(dfm_old, dfm_new, t.fcst = 355, target.vars = "orders", standardized = TRUE)
  revision_m_std <- unname(res_m_std$y_new - res_m_std$y_old)
  sum_impact_std <- sum(res_m_std$news_df$impact)
  expect_lt(abs(revision_m_std - sum_impact_std), 1e-10)

  # Check that released variables contributed news
  idx_rel <- match(c("new_cars", "pms_pmi"), res_m$news_df$series)
  expect_true(all(res_m$news_df$news[idx_rel] != 0))

  # Check gains exist for released variables
  expect_true(any(res_m$news_df$gain != 0))

  # Check that gain produces results.
  rel_idx <- which(!is.na(res_m$news_df$actual))
  expect_equal(res_m$news_df$impact[rel_idx], res_m$news_df$news[rel_idx] * res_m$news_df$gain[rel_idx])

})

test_that("news works with MQ small model for quarterly target", {
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
  X_old <- X_small[-1, ]
  X_new <- X_small[-1, ]

  # Create releases from observed variables
  X_old[355, "ip_tot_cstr"] <- NA
  X_old[356, "new_cars"] <- NA
  X_old[356, "pms_pmi"] <- NA
  X_old[356, "euro325"] <- NA
  X_old[356, "capacity"] <- NA

  # Fit models with mixed-frequency
  dfm_old <- DFM(X_old, r = 2, p = 2, quarterly.vars = quarterly.vars, max.missing = 1)
  dfm_new <- DFM(X_new, r = 2, p = 2, quarterly.vars = quarterly.vars, max.missing = 1)

  # Target: 'orders' at t=355 (X_imp coordinates) - naturally missing at end of sample
  # This is a proper nowcast scenario where target is NOT observed
  res_m <- news(dfm_old, dfm_new, t.fcst = 355, target.vars = "gdp")

  expect_s3_class(res_m, "dfm.news")
  expect_true(is.numeric(res_m$y_old))
  expect_true(is.numeric(res_m$y_new))

  # revision should equal sum of impacts
  revision_m <- unname(res_m$y_new - res_m$y_old)
  sum_impact <- sum(res_m$news_df$impact)
  expect_equal(revision_m, sum_impact, tolerance = 1e-10)

  # Standardized scale should match more tightly
  res_m_std <- news(dfm_old, dfm_new, t.fcst = 355, target.vars = "gdp", standardized = TRUE)
  revision_m_std <- unname(res_m_std$y_new - res_m_std$y_old)
  sum_impact_std <- sum(res_m_std$news_df$impact)
  expect_lt(abs(revision_m_std - sum_impact_std), 1e-10)

  # Check that released variables contributed news
  idx_rel <- match(c("ip_tot_cstr", "new_cars", "pms_pmi", "euro325", "capacity"), res_m$news_df$series)
  expect_true(all(res_m$news_df$news[idx_rel] != 0))

  # Check gains exist for released variables
  expect_true(any(res_m$news_df$gain != 0))

  # # Check that gain produces results.
  rel_idx <- which(!is.na(res_m$news_df$actual))
  expect_equal(res_m$news_df$impact[rel_idx], res_m$news_df$news[rel_idx] * res_m$news_df$gain[rel_idx])

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

  # For monthly targets, revision should equal sum of impacts
  revision_m <- unname(res_m$y_new - res_m$y_old)
  sum_impact <- sum(res_m$news_df$impact)
  expect_equal(revision_m, sum_impact, tolerance = 1e-4)

  # Standardized scale should match more tightly
  res_m_std <- news(dfm_old, dfm_new, t.fcst = 355, target.vars = "orders", standardized = TRUE)
  revision_m_std <- unname(res_m_std$y_new - res_m_std$y_old)
  sum_impact_std <- sum(res_m_std$news_df$impact)
  expect_lt(abs(revision_m_std - sum_impact_std), 1e-4)

  # Check gains exist for released variables
  expect_true(any(res_m$news_df$gain != 0))

  # Check that gain produces results.
  rel_idx <- which(!is.na(res_m$news_df$actual))
  expect_equal(res_m$news_df$impact[rel_idx], res_m$news_df$news[rel_idx] * res_m$news_df$gain[rel_idx])

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
  sum_impact <- sum(res$news_df$impact)
  expect_equal(revision, sum_impact, tolerance = 1e-8)


  rel_idx <- which(!is.na(res$news_df$actual))
  expect_equal(res$news_df$impact[rel_idx], res$news_df$news[rel_idx] * res$news_df$gain[rel_idx])
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

  X_old <- X_small[-1, ]
  X_new <- X_small[-1, ]

  # Create releases from observed variables
  X_old[355, "ip_tot_cstr"] <- NA
  X_old[356, "new_cars"] <- NA
  X_old[356, "pms_pmi"] <- NA
  X_old[356, "euro325"] <- NA
  X_old[356, "capacity"] <- NA

  dfm_old <- DFM(X_old, r = 2, p = 2, quarterly.vars = quarterly.vars, idio.ar1 = TRUE, max.missing = 1)
  dfm_new <- DFM(X_new, r = 2, p = 2, quarterly.vars = quarterly.vars, idio.ar1 = TRUE, max.missing = 1)

  res_m <- news(dfm_old, dfm_new, t.fcst = 356, target.vars = "orders")

  expect_s3_class(res_m, "dfm.news")
  revision_m <- unname(res_m$y_new - res_m$y_old)
  sum_impact <- sum(res_m$news_df$impact)
  expect_equal(revision_m, sum_impact, tolerance = 1e-6)

  rel_idx <- which(!is.na(res_m$news_df$actual))
  expect_equal(res_m$news_df$impact[rel_idx], res_m$news_df$news[rel_idx] * res_m$news_df$gain[rel_idx])
})

test_that("news handles errors properly", {
  set.seed(101)
  X <- matrix(rnorm(80), nrow = 20)
  colnames(X) <- paste0("v", seq_len(ncol(X)))
  X_old <- X
  X_new <- X
  X_old[15, 2] <- NA
  X_new[15, 2] <- X[15, 2]

  dfm_old <- DFM(X_old, r = 1, p = 1, em.method = "none")

  expect_error(news(dfm_old, X_new, t.fcst = 25), "t.fcst is out of bounds")
  expect_error(news(dfm_old, X_new, t.fcst = 10, target.vars = "nonexistent"))
})

test_that("news_list works as expected", {
  set.seed(303)
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

  res_list <- news(dfm_old, dfm_new, t.fcst = 20)
  expect_s3_class(res_list, "dfm.news_list")
  expect_equal(length(res_list), ncol(X))

  for(res in res_list$news) {
    expect_s3_class(res, "dfm.news")
  }
})
