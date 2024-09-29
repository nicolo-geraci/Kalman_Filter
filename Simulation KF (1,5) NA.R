rm(list = ls())
#install.packages("missMDA")
library(missMDA)
set.seed(2)  # Set seed for reproducibility

#--------------------- READ ME

# This script implements a Kalman Filter Maximum Likelihood estimation based 
# on 5 observed (simulated) series and one latent state following an AR(1) dynamic.
# The KF iterations are initialized via Principal Component Analysis. Early estimates of the auto-regressive 
# coefficient governing the latent dynamic is estimated fitting an AR(1) model on the scores of the most informative
# Principal Component.
# In its current version, the Kalman Filter MLE function skips the update step contribution from a specific
# observed series when meets a missing value. The PCA function, instead, deals with dataset where missing values are  
# present using iterative NA values imputation.

# Map of the code: 1. DATA GENERATION: both of the Latent State and the five observed series
#                    1.1. Introduce NA values in the dataset.
#                  2. KALMAN FILTER ESTIMATION
#                    2. 1. Define KF MLE function
#                    2. 2. Recover early estimates of the parameters to be estimated (PCA)
#                    2. 3. Run optimization
#                    2. 4. Print results
#                  3. KALMAN SMOOTHING
#                  4. FINAL RESULTS


# -------------------- 1. Data generation

# Generate Latent State Series though an AR(1) process
{
  n <- 150
  phi <- 0.96
  sigma_eta <- 0.5  # Standard deviation of the process noise
  X <- numeric(n)
  X[1] <- rnorm(1)  # Initial state
  
  for (t in 2:n) {
    X[t] <- phi * X[t-1] + rnorm(1, sd = sigma_eta)
  }
}

# Generate Five Noisy Observed Series from the Latent State
{
  sigma_epsilon <- c(1.5, 1.5, 1.2, 1.0, 0.8)   # Standard deviations of the noise for each series
  alpha <- c(2.0, 1.5, -1.0, 0.8, 0.6)          # Coefficients for each series
  
  Y1 <- numeric(n)
  Y2 <- numeric(n)
  Y3 <- numeric(n)
  Y4 <- numeric(n)
  Y5 <- numeric(n)
    correlation_noise <- 0  # this switcher can be used to control correlation between noises of the observed series 
  
  for (t in 1:n) {
    epsilon1 <- rnorm(1, sd = sigma_epsilon[1])
    epsilon2 <- correlation_noise * epsilon1 + sqrt(1 - correlation_noise^2) * rnorm(1, sd = sigma_epsilon[2])
    epsilon3 <- correlation_noise * epsilon1 + sqrt(1 - correlation_noise^2) * rnorm(1, sd = sigma_epsilon[3])
    epsilon4 <- correlation_noise * epsilon1 + sqrt(1 - correlation_noise^2) * rnorm(1, sd = sigma_epsilon[4])
    epsilon5 <- correlation_noise * epsilon1 + sqrt(1 - correlation_noise^2) * rnorm(1, sd = sigma_epsilon[5])
    
    Y1[t] <- alpha[1] * X[t] + epsilon1
    Y2[t] <- alpha[2] * X[t] + epsilon2
    Y3[t] <- alpha[3] * X[t] + epsilon3
    Y4[t] <- alpha[4] * X[t] + epsilon4
    Y5[t] <- alpha[5] * X[t] + epsilon5
  }
}

# Compare pair-wise correlation 
{
cat("Correlation between Y1 and Y2:", cor(Y1, Y2), "\n")
cat("Correlation between Y1 and Y3:", cor(Y1, Y3), "\n")
cat("Correlation between Y1 and Y4:", cor(Y1, Y4), "\n")
cat("Correlation between Y1 and Y5:", cor(Y1, Y5), "\n")
cat("Correlation between Y2 and Y3:", cor(Y2, Y3), "\n")
cat("Correlation between Y2 and Y4:", cor(Y2, Y4), "\n")
cat("Correlation between Y2 and Y5:", cor(Y2, Y5), "\n")
cat("Correlation between Y3 and Y4:", cor(Y3, Y4), "\n")
cat("Correlation between Y3 and Y5:", cor(Y3, Y5), "\n")
cat("Correlation between Y4 and Y5:", cor(Y4, Y5), "\n")
}
# PLOT
#par(mfrow = c(1, 1))
{plot(Y1, type = "l", col = "blue", main = "Observed Series Y1 to Y5", ylab = "Series", xlab = "Time", sub = "Latent state in black")
lines(Y2, col = "red")
lines(Y3, col = "green")
lines(Y4, col = "purple")
lines(Y5, col = "orange")
lines(X, lwd = 3, col = "black")}



# -------------------- 1.1 Introduce missing values
set.seed(2)  # Set seed for reproducibility

# Define the fraction of missing values
missing_fraction <- 0.05  # 5% missing values

# Number of missing values to insert
num_missing <- round(missing_fraction * length(Y1))

# Randomly select indices for missing values
{
missing_indices_Y1 <- sample(1:length(Y1), num_missing, replace = FALSE)
missing_indices_Y2 <- sample(1:length(Y2), num_missing, replace = FALSE)
missing_indices_Y3 <- sample(1:length(Y3), num_missing, replace = FALSE)
missing_indices_Y4 <- sample(1:length(Y4), num_missing, replace = FALSE)
missing_indices_Y5 <- sample(1:length(Y5), num_missing, replace = FALSE)
}
# Insert missing values
{
Y1[missing_indices_Y1] <- NA
Y2[missing_indices_Y2] <- NA
Y3[missing_indices_Y3] <- NA
Y4[missing_indices_Y4] <- NA
Y5[missing_indices_Y5] <- NA
}
# Verify the insertion of missing values
{
cat("Number of missing values in Y1:", sum(is.na(Y1)), "\n")
cat("Number of missing values in Y2:", sum(is.na(Y2)), "\n")
cat("Number of missing values in Y3:", sum(is.na(Y3)), "\n")
cat("Number of missing values in Y4:", sum(is.na(Y4)), "\n")
cat("Number of missing values in Y5:", sum(is.na(Y5)), "\n")
}


# -------------------- 2. KFM Estimates

# ---------------------- 2.1 Define KF MLE function

kalman_log_likelihood <- function(par, Y1, Y2, Y3, Y4, Y5) {
  # Parameters to be estimated
  phi <- par[1]
  alpha1 <- par[2]
  alpha2 <- par[3]
  alpha3 <- par[4]
  alpha4 <- par[5]
  alpha5 <- par[6]
  sigma_eta <- exp(par[7])  # variance, hence exp to ensure positivity
  sigma_epsilon1 <- exp(par[8])  # variance for the first observed series
  sigma_epsilon2 <- exp(par[9])  # ...
  sigma_epsilon3 <- exp(par[10])  
  sigma_epsilon4 <- exp(par[11])  
  sigma_epsilon5 <- exp(par[12])  # variance for the fifth observed series
  
  # Length of the time series
  n <- length(Y1)
  
  # Initialize Kalman filter variables
  X_pred <- numeric(n)
  P_pred <- numeric(n)
  X_upd <- numeric(n)
  P_upd <- numeric(n)
  
  # Initial state estimates
  X_upd[1] <- pca1[1]  # Alternatively, Assume prior mean of state is 0 ...
  P_upd[1] <- var(pca1)  # ...and set large initial variance
  
  log_likelihood <- 0
  
  for (t in 2:n) {
    # Predict
    X_pred[t] <- phi * X_upd[t-1]
    P_pred[t] <- phi^2 * P_upd[t-1] + sigma_eta^2
    
    # List to store valid updates
    X_upd_list <- c()
    P_upd_list <- c()
    K_list <- c()
    
    # Update for the first observation equation
    if (!is.na(Y1[t])) {
      K1 <- P_pred[t] / (P_pred[t] + sigma_epsilon1^2)
      X_upd1 <- X_pred[t] + K1 * (Y1[t] - alpha1 * X_pred[t])
      P_upd1 <- (1 - K1) * P_pred[t]
      X_upd_list <- c(X_upd_list, X_upd1)
      P_upd_list <- c(P_upd_list, P_upd1)
      K_list <- c(K_list, K1)
      
      # Contribution to the log-likelihood
      log_likelihood <- log_likelihood -
        0.5 * (log(2 * pi * (P_pred[t] + sigma_epsilon1^2)) + (Y1[t] - alpha1 * X_pred[t])^2 / (P_pred[t] + sigma_epsilon1^2))
    }
    
    # Update for the second observation equation
    if (!is.na(Y2[t])) {
      K2 <- P_pred[t] / (P_pred[t] + sigma_epsilon2^2)
      X_upd2 <- X_pred[t] + K2 * (Y2[t] - alpha2 * X_pred[t])
      P_upd2 <- (1 - K2) * P_pred[t]
      X_upd_list <- c(X_upd_list, X_upd2)
      P_upd_list <- c(P_upd_list, P_upd2)
      K_list <- c(K_list, K2)
      
      # Contribution to the log-likelihood
      log_likelihood <- log_likelihood -
        0.5 * (log(2 * pi * (P_pred[t] + sigma_epsilon2^2)) + (Y2[t] - alpha2 * X_pred[t])^2 / (P_pred[t] + sigma_epsilon2^2))
    }
    
    # Update for the third observation equation
    if (!is.na(Y3[t])) {
      K3 <- P_pred[t] / (P_pred[t] + sigma_epsilon3^2)
      X_upd3 <- X_pred[t] + K3 * (Y3[t] - alpha3 * X_pred[t])
      P_upd3 <- (1 - K3) * P_pred[t]
      X_upd_list <- c(X_upd_list, X_upd3)
      P_upd_list <- c(P_upd_list, P_upd3)
      K_list <- c(K_list, K3)
      
      # Contribution to the log-likelihood
      log_likelihood <- log_likelihood -
        0.5 * (log(2 * pi * (P_pred[t] + sigma_epsilon3^2)) + (Y3[t] - alpha3 * X_pred[t])^2 / (P_pred[t] + sigma_epsilon3^2))
    }
    
    # Update for the fourth observation equation
    if (!is.na(Y4[t])) {
      K4 <- P_pred[t] / (P_pred[t] + sigma_epsilon4^2)
      X_upd4 <- X_pred[t] + K4 * (Y4[t] - alpha4 * X_pred[t])
      P_upd4 <- (1 - K4) * P_pred[t]
      X_upd_list <- c(X_upd_list, X_upd4)
      P_upd_list <- c(P_upd_list, P_upd4)
      K_list <- c(K_list, K4)
      
      # Contribution to the log-likelihood
      log_likelihood <- log_likelihood -
        0.5 * (log(2 * pi * (P_pred[t] + sigma_epsilon4^2)) + (Y4[t] - alpha4 * X_pred[t])^2 / (P_pred[t] + sigma_epsilon4^2))
    }
    
    # Update for the fifth observation equation
    if (!is.na(Y5[t])) {
      K5 <- P_pred[t] / (P_pred[t] + sigma_epsilon5^2)
      X_upd5 <- X_pred[t] + K5 * (Y5[t] - alpha5 * X_pred[t])
      P_upd5 <- (1 - K5) * P_pred[t]
      X_upd_list <- c(X_upd_list, X_upd5)
      P_upd_list <- c(P_upd_list, P_upd5)
      K_list <- c(K_list, K5)
      
      # Contribution to the log-likelihood
      log_likelihood <- log_likelihood -
        0.5 * (log(2 * pi * (P_pred[t] + sigma_epsilon5^2)) + (Y5[t] - alpha5 * X_pred[t])^2 / (P_pred[t] + sigma_epsilon5^2))
    }
    
    # Combine updates if any are present
    if (length(X_upd_list) > 0) {
      combined_gain <- sum(K_list)
      X_upd[t] <- sum(K_list * X_upd_list) / combined_gain
      P_upd[t] <- sum(K_list * P_upd_list) / combined_gain
    } else {
      X_upd[t] <- X_pred[t]
      P_upd[t] <- P_pred[t]
    }
  }
  
  assign("X_filtered", X_upd, envir = .GlobalEnv)
  assign("P_filtered", P_upd, envir = .GlobalEnv) # these variables can be used to run KF smoothing algorithm
  assign("P_pred", P_pred, envir = .GlobalEnv)
  
  return(-log_likelihood)  # return negative log-likelihood for minimization
}



# ---------------------- 2.2 Build initial guesses for the parameters

# 2.2.1 Perform PCA on observed series to have early estimate of \phi

{df <- data.frame(Y1, Y2, Y3, Y4, Y5) # Merge observed series into a dataframe
ncp <- estim_ncpPCA(df, scale = TRUE) # Estimate the number of principal components (ncp) to retain
data_imputed_missMDA <- imputePCA(df, ncp = ncp$ncp)$completeObs # Impute missing values with iterative PCA-based imputation
pca_result_missMDA <- prcomp(data_imputed_missMDA, scale. = TRUE) # Perform PCA on the imputed data
summary(pca_result_missMDA) # Summary of PCA results
# View the principal components
#pca_result_missMDA$x

pca1 <- pca_result_missMDA$x[,1]
ar_model <- ar(pca1, order.max = 1, method = "mle") #fit AR(1) on PC1 scores
ar_pca <- ar_model$ar   # AR coefficient
# # Plot variance explained by each principal component
# plot(pca_result_missMDA, type = "lines", main = "Scree Plot")
}


initial_params <- c(ar_pca, 1.0, 1.0, 1.0, 1.0, 1.0, log(1), log(1), log(1), log(1), log(1), log(1))

# ---------------------- 2.3 RUN OPTIMIZATION

mle_results <- optim(initial_params, kalman_log_likelihood, Y1 = Y1, Y2 = Y2, Y3 = Y3, Y4 = Y4, Y5 = Y5, hessian = TRUE, method = "BFGS")

# Extract estimated parameters
{
  estimated_phi <- mle_results$par[1]
  estimated_alpha1 <- mle_results$par[2]
  estimated_alpha2 <- mle_results$par[3]
  estimated_alpha3 <- mle_results$par[4]
  estimated_alpha4 <- mle_results$par[5]
  estimated_alpha5 <- mle_results$par[6]
  estimated_sigma_eta <- exp(mle_results$par[7])
  estimated_sigma_epsilon1 <- exp(mle_results$par[8])
  estimated_sigma_epsilon2 <- exp(mle_results$par[9])
  estimated_sigma_epsilon3 <- exp(mle_results$par[10])
  estimated_sigma_epsilon4 <- exp(mle_results$par[11])
  estimated_sigma_epsilon5 <- exp(mle_results$par[12])
}

# Print results
{
cat("Estimated phi:", estimated_phi, "versus true value:", phi, "\n")
cat("Estimated alpha1:", estimated_alpha1,"versus true value:", alpha[1], "\n")
cat("Estimated alpha2:", estimated_alpha2,"versus true value:", alpha[2], "\n")
cat("Estimated alpha3:", estimated_alpha3,"versus true value:", alpha[3], "\n")
cat("Estimated alpha4:", estimated_alpha4,"versus true value:", alpha[4], "\n")
cat("Estimated alpha5:", estimated_alpha5,"versus true value:", alpha[5], "\n")

cat("Estimated sigma_eta:", estimated_sigma_eta,"versus true value:", sigma_eta, "\n")
cat("Estimated sigma_epsilon1:", estimated_sigma_epsilon1,"versus true value:", sigma_epsilon[1], "\n")
cat("Estimated sigma_epsilon2:", estimated_sigma_epsilon2,"versus true value:", sigma_epsilon[2], "\n")
cat("Estimated sigma_epsilon3:", estimated_sigma_epsilon3,"versus true value:", sigma_epsilon[3], "\n")
cat("Estimated sigma_epsilon4:", estimated_sigma_epsilon4,"versus true value:", sigma_epsilon[4], "\n")
cat("Estimated sigma_epsilon5:", estimated_sigma_epsilon5,"versus true value:", sigma_epsilon[5], "\n")
}



# Plot True vs Filtered data
plot(X, type = "l", col = "black", xlab = "Time", ylab = "Values", main = "True vs Filtered data" )
lines(X_filtered, col = "red")
#legend()



# ---------------------- 3. RUN KALMAN SMOOTHING 
# (Rauch-Tung-Striebel smoother)
kalman_smoothing <- function(phi, X_filtered, P_filtered, P_pred) {
  n <- length(X_filtered)
  
  # Initialize smoothed estimates
  X_smooth <- numeric(n)
  P_smooth <- numeric(n)
  
  # Start from the last filtered estimates
  X_smooth[n] <- X_filtered[n]
  P_smooth[n] <- P_filtered[n]
  
  # Backward pass
  for (t in (n-1):1) {
    # Smoothing gain
    J <- P_filtered[t] * phi / P_pred[t+1]
    
    # Smoothed state
    X_smooth[t] <- X_filtered[t] + J * (X_smooth[t+1] - phi * X_filtered[t])
    
    # Smoothed covariance
    P_smooth[t] <- P_filtered[t] + J * (P_smooth[t+1] - P_pred[t+1]) * J
  }
  
  return(list(X_smooth = X_smooth, P_smooth = P_smooth))
}

smoothing_results <- kalman_smoothing(estimated_phi, X_filtered, P_filtered, P_pred)

# Smoothed states
{
X_smooth <- smoothing_results$X_smooth
P_smooth <- smoothing_results$P_smooth
}



# ------------------- 4. FINAL RESULTS

plot(X, type ="l",col = "black", lwd = 2, xlab="Time", main = "True vs Estimated data", ylab ="Values")
lines(X_filtered, col="red", lwd = 0.7)
lines(X_smooth, col = "blue", lwd=0.7)


# Compute mean forecast error for both Filtered and Smoothed Data

forecast_error_filt <- mean((X-X_filtered)^2)
forecast_error_smooth <-mean((X-X_smooth)^2)

cat("mean squared forecast error of filtered series:", forecast_error_filt, "\n")
cat("mean squared forecast error of smoothed series:", forecast_error_smooth, "\n")


