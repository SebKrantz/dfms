library(collapse)

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
