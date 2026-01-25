###########################################
# Test R vs Julia DFM Implementations
# Comparing EM_DFM_SS_restrMQ and news()
###########################################

library(dfms)
library(xts)
library(magrittr)

# BM14 Replication Data. Constructing the database:
BM14 <- merge(BM14_M, BM14_Q)
BM14[, BM14_Models$log_trans] %<>% log()
BM14[, BM14_Models$freq == "M"] %<>% diff()
BM14[, BM14_Models$freq == "Q"] %<>% diff(3)

# Select 12 variables (10 monthly + 2 quarterly at the end)
X <- BM14[, c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99)]
quarterly.vars <- colnames(X)[11:12]  # invest, productivity

cat("Data dimensions:", dim(X), "\n")
cat("Quarterly variables:", quarterly.vars, "\n")
cat("Column names:", colnames(X), "\n\n")

###########################################
# Test 1: MQ Model (EM_DFM_SS_restrMQ)
###########################################
cat("=== Test 1: MQ Model (r=2, p=3) ===\n")

dfm_mq <- DFM(X, r = 2, p = 3,
              quarterly.vars = quarterly.vars,
              em.method = "BM")

cat("Converged:", dfm_mq$converged, "\n")
cat("Iterations:", length(dfm_mq$loglik), "\n")
cat("Factor dimensions:", dim(dfm_mq$F_pca), "\n\n")

# Save parameters for Julia comparison
cat("--- Parameters for Julia comparison ---\n")
cat("C matrix dimensions:", dim(dfm_mq$C), "\n")
cat("A matrix dimensions:", dim(dfm_mq$A), "\n")
cat("Q matrix dimensions:", dim(dfm_mq$Q), "\n")
cat("R matrix dimensions:", dim(dfm_mq$R), "\n\n")

# Sample of C matrix (first 5 rows)
cat("C matrix (first 5 rows, first 2 cols):\n")
print(round(dfm_mq$C[1:5, 1:2], 6))

###########################################
# Test 2: News Decomposition (MQ model)
###########################################
cat("\n=== Test 2: News Decomposition (MQ) ===\n")

# Create old and new vintages
X_new <- X
X_old <- X

# Simulate data releases: remove recent monthly observations
X_old[350, 1] <- NA  # ip_total
X_old[351, 1] <- NA

cat("Releases at t=350,351 for variable 1 (ip_total)\n")

# Fit models on both vintages
dfm_old <- DFM(X_old, r = 2, p = 3, quarterly.vars = quarterly.vars, em.method = "BM")
dfm_new <- DFM(X_new, r = 2, p = 3, quarterly.vars = quarterly.vars, em.method = "BM")

# Run news decomposition for a monthly target
t_fcst <- 355
target_var <- "ret_turnover_defl"  # Variable 3 (monthly)

cat("Target time:", t_fcst, "\n")
cat("Target variable:", target_var, "\n\n")

res_news <- news(dfm_old, dfm_new, t.fcst = t_fcst, target.vars = target_var)

cat("--- News Results ---\n")
cat("Old forecast:", res_news$y_old, "\n")
cat("New forecast:", res_news$y_new, "\n")
cat("Revision:", res_news$y_new - res_news$y_old, "\n")
cat("Sum of news:", sum(res_news$singlenews), "\n")
cat("Difference (should be ~0):", (res_news$y_new - res_news$y_old) - sum(res_news$singlenews), "\n\n")

cat("Non-zero news contributions:\n")
nz <- which(res_news$singlenews != 0)
if(length(nz) > 0) {
  for(i in nz) {
    cat(sprintf("  Variable %d (%s): %.10f\n", i, colnames(X)[i], res_news$singlenews[i]))
  }
}

###########################################
# Test 3: MQ + AR1 Model (if available)
###########################################
cat("\n=== Test 3: MQ + AR1 Idiosyncratic Model ===\n")

dfm_mq_ar1 <- tryCatch({
  DFM(X, r = 2, p = 3,
      quarterly.vars = quarterly.vars,
      idio.ar1 = TRUE,
      em.method = "BM")
}, error = function(e) {
  cat("Error fitting MQ+AR1 model:", e$message, "\n")
  NULL
})

if(!is.null(dfm_mq_ar1)) {
  cat("Converged:", dfm_mq_ar1$converged, "\n")
  cat("Iterations:", dfm_mq_ar1$em.iterations, "\n")

  # Check rho (AR1 coefficients)
  cat("AR1 coefficients (rho):\n")
  print(round(dfm_mq_ar1$rho, 4))

  # News with AR1 model
  cat("\n--- News with AR1 model ---\n")
  dfm_old_ar1 <- DFM(X_old, r = 2, p = 3, quarterly.vars = quarterly.vars,
                     idio.ar1 = TRUE, em.method = "BM")
  dfm_new_ar1 <- DFM(X_new, r = 2, p = 3, quarterly.vars = quarterly.vars,
                     idio.ar1 = TRUE, em.method = "BM")

  res_news_ar1 <- news(dfm_old_ar1, dfm_new_ar1, t.fcst = t_fcst, target.vars = target_var)

  cat("Old forecast:", res_news_ar1$y_old, "\n")
  cat("New forecast:", res_news_ar1$y_new, "\n")
  cat("Revision:", res_news_ar1$y_new - res_news_ar1$y_old, "\n")
  cat("Sum of news:", sum(res_news_ar1$singlenews), "\n")
  cat("Difference (should be ~0):", (res_news_ar1$y_new - res_news_ar1$y_old) - sum(res_news_ar1$singlenews), "\n")
}

###########################################
# Export data for Julia comparison
###########################################
cat("\n=== Exporting for Julia Comparison ===\n")

# Save the parameters from MQ model for Julia to load
saveRDS(list(
  C = dfm_mq$C,
  A = dfm_mq$A,
  Q = dfm_mq$Q,
  R = dfm_mq$R,
  F = dfm_mq$F_pca,
  Mx = dfm_mq$Mx,
  Wx = dfm_mq$Wx,
  r = dfm_mq$r,
  p = dfm_mq$p,
  quarterly.vars = quarterly.vars,
  news_result = list(
    y_old = res_news$y_old,
    y_new = res_news$y_new,
    singlenews = res_news$singlenews
  )
), "misc/BM2014/R_dfm_mq_results.rds")

cat("Results saved to misc/BM2014/R_dfm_mq_results.rds\n")
cat("\nDone!\n")
