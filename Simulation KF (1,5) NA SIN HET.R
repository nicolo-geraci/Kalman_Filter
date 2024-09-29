rm(list = ls())
#install.packages("missMDA")
library(missMDA)
set.seed(2)  # Set seed for reproducibility

#--------------------- READ ME

# This script implements a Kalman Filter Maximum Likelihood estimation based 
# on 5 observed (simulated) series and one latent state following a sinusoidal dynamic defined over a period t = 40.
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

# Further Considerations: Under this specification the model can recognize the sinusoidal pattern
# of the latent variable with the correct frequency. However, the precision of the estimates
# regarding the oscillations is greatly influenced by both the number of missing values and
# the noise of the different series.
# The script currently imposes 20% missing values in each series and a good correlation between
# the noises of the observed series. For instance, if the correlation between contemporary noises
# is set to zero, the MSFE would increase by approx. 80%.

# -------------------- 1. Data generation

# Generate Latent State Series though a sinusoidal function
n <- 150

# Sinusoidal function parameters
{
A <- 1     # Amplitude (controls the peak value)
T <- 50    # Period (controls the length of each cycle)
phi <- 0   # Phase shift (optional)
sigma_eta <- 0.5

# Generate the time series following a sinusoidal function
t <- 1:n
X_sin <- 2.3 * sin(2 * pi * t / T + phi) + rnorm(n, sd = sigma_eta)

# Plot the sinusoidal time series
plot(X_sin, type = "l", col = "blue", main = "Sinusoidal Time Series", ylab = "Value", xlab = "Time")
}
X <- X_sin
# Generate Five Noisy Observed Series from the Latent State using a sinusoidal volatility
{
  alpha <- c(2.0, 1.5, -1.0, 0.8, 0.6)          # Coefficients for each series
  
  Y1 <- numeric(n)
  Y2 <- numeric(n)
  Y3 <- numeric(n)
  Y4 <- numeric(n)
  Y5 <- numeric(n)
    correlation_noise <- 0  # this switcher can be used to control correlation between noises of the observed series 
}
# Sinusoidal function for time-varying variance
    {
      A <- 1.0  # Baseline volatility high
      B <- 1.1  # Amplitude high (controls fluctuation in volatility)
      T <- 60   # Period of the sinusoidal wave (controls frequency)
      phi_shift <- 0  # Phase shift
      C <- 0.5  # Baseline volatility low
      D <- 0.9  # Amplitude low
    }
    # Generate the sinusoidal volatility function
{
    sinusoidal_volatility_high <- A + B * abs(sin((2 * pi * (1:n) / T) + phi_shift))
    sinusoidal_volatility_low <- C + D * abs(sin((2 * pi * (1:n) / T) + phi_shift))
    mean_sin_variance_h = mean(sinusoidal_volatility_high)
    mean_sin_variance_l = mean(sinusoidal_volatility_low)
    print(mean_sin_variance_h)
    print(mean_sin_variance_l)
    plot(sinusoidal_volatility_high, type ="l", xlab = "Time", ylim = c(0.5,2.5),
         ylab = "Variance", main = "Example time-varying volatility", col = "blue")
    lines(sinusoidal_volatility_low, col = "red")
    # these are the periodically evolving volatilities from which innovation terms will be drawn
    # during the generation of artificial series
    # decide weather use them as standard dev or variances!
    sigma_epsilon <- c(mean_sin_variance_h, mean_sin_variance_h, mean_sin_variance_h, mean_sin_variance_l, mean_sin_variance_l)  # Standard deviations of the noise for each series
}
  for (t in 1:n) {
    epsilon1 <- rnorm(1, sd = sinusoidal_volatility_high[t])
    epsilon2 <- correlation_noise * epsilon1 + sqrt(1 - correlation_noise^2) * rnorm(1, sd = sinusoidal_volatility_high[t])
    epsilon3 <- correlation_noise * epsilon1 + sqrt(1 - correlation_noise^2) * rnorm(1, sd = sinusoidal_volatility_high[t])
    epsilon4 <- correlation_noise * epsilon1 + sqrt(1 - correlation_noise^2) * rnorm(1, sd = sinusoidal_volatility_low[t])
    epsilon5 <- correlation_noise * epsilon1 + sqrt(1 - correlation_noise^2) * rnorm(1, sd = sinusoidal_volatility_low[t])
    
    Y1[t] <- alpha[1] * X[t] + epsilon1
    Y2[t] <- alpha[2] * X[t] + epsilon2
    Y3[t] <- alpha[3] * X[t] + epsilon3
    Y4[t] <- alpha[4] * X[t] + epsilon4
    Y5[t] <- alpha[5] * X[t] + epsilon5
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

dfx <- data.frame(Y1, Y2, Y3, Y4, Y5) #store original balanced df for post-estimation diagnostics

# -------------------- 1.1 Introduce missing values
set.seed(2)  # Set seed for reproducibility

# Define the fraction of missing values
missing_fraction <- 0.20  # 20% missing values

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

kalman_log_likelihood <- function(par, Y1, Y2, Y3, Y4, Y5, return_all = FALSE) {
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
  
  
  # Return a list containing the log-likelihood, filtered states, and uncertainties
  if (return_all) {
    return(list(log_likelihood = -log_likelihood, X_filtered = X_upd, P_filtered = P_upd, P_pred = P_pred))
  }
  
  # Otherwise, return just the negative log-likelihood for optim
  return(-log_likelihood)
}



# ---------------------- 2.2 Build initial guesses for the parameters

# 2.2.1 Perform PCA on observed series to have early estimate of \phi

{df <- data.frame(Y1, Y2, Y3, Y4, Y5) # Merge observed series into a dataframe
ncp <- estim_ncpPCA(df, scale = TRUE) # Estimate the number of principal components (ncp) to retain
data_imputed_missMDA <- imputePCA(df, ncp = ncp$ncp)$completeObs # Impute missing values with iterative PCA-based imputation
pca_result_missMDA <- prcomp(data_imputed_missMDA, scale. = TRUE) # Perform PCA on the imputed data
summary(pca_result_missMDA) # Summary of PCA results
df_imputed <- data.frame(data_imputed_missMDA)
# View the principal components
#pca_result_missMDA$x

pca1 <- pca_result_missMDA$x[,1]
ar_model <- ar(pca1, order.max = 1, method = "mle") #fit AR(1) on PC1 scores
ar_pca <- ar_model$ar   # AR coefficient
# # Plot variance explained by each principal component
# plot(pca_result_missMDA, type = "lines", main = "Scree Plot")
}


#------ next lines can use regression between PC1 and the series to initialize values of ALPHA...but issue of Rotation Matrix!!
# loading_Y1_PC1 <- pca_result_missMDA$rotation[1, 1]
# loading_PCA <- pca_result_missMDA$rotation
# print(loading_Y1_PC1)
# rm1 <- lm(Y1 ~ pca1, data = df_imputed)
# summary(rm1)
# fitted_values1 <- predict(rm1)
# plot(pca1, Y1, main = "Regression of Y1 on PC1", xlab = "PC1", ylab = "Y1", col = "blue")
# lines(pca1, fitted_values1, col = "red", lwd = 2)  # Add regression line



initial_params <- c(ar_pca, 1.0, 1.0, -1.0, 1.0, 1.0, log(1), log(1), log(1), log(1), log(1), log(1))

# ---------------------- 2.3 RUN OPTIMIZATION

mle_results <- optim(initial_params, kalman_log_likelihood, Y1 = Y1, Y2 = Y2, Y3 = Y3, Y4 = Y4, Y5 = Y5, hessian = TRUE, method = "BFGS")
optimized_params <- mle_results$par
final_result <- kalman_log_likelihood(optimized_params, Y1 = Y1, Y2 = Y2, Y3 = Y3, Y4 = Y4, Y5 = Y5, return_all = TRUE)
# Extract the results
X_filtered <- final_result$X_filtered
P_filtered <- final_result$P_filtered
P_pred <- final_result$P_pred

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
X_NA_filtered <- X_filtered #so that further below we can compare them with X_filtered from balanced DF


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





############################## COMPARISON



# now compare results with the use of different df: original balanced panel and PCA imputed

### 1. try estimating KF on original balanced data
pca_result_bal <- prcomp(dfx, scale. = TRUE) # Perform PCA on the original complete data
summary(pca_result_bal) # Summary of PCA results
pca1 <- pca_result_bal$x[,1]
ar_model <- ar(pca1, order.max = 1, method = "mle") #fit AR(1) on PC1 scores
ar_pca <- ar_model$ar
initial_params <- c(ar_pca, 1.0, 1.0, -1.0, 1.0, 1.0, log(1), log(1), log(1), log(1), log(1), log(1))
Yx1 <- dfx$Y1
Yx2 <- dfx$Y2
Yx3 <- dfx$Y3
Yx4 <- dfx$Y4
Yx5 <- dfx$Y5
mle_results_B <- optim(initial_params, kalman_log_likelihood, Y1 = Yx1, Y2 = Yx2, Y3 = Yx3, Y4 = Yx4, Y5 = Yx5, hessian = TRUE, method = "BFGS")
optimized_params_B <- mle_results_B$par
final_resultB <- kalman_log_likelihood(optimized_params_B, Y1 = Yx1, Y2 = Yx2, Y3 = Yx3, Y4 = Yx4, Y5 = Yx5, return_all = TRUE)
# Extract the results
X_filtered_B <- final_resultB$X_filtered
P_filtered_B <- final_resultB$P_filtered
P_pred_B <- final_resultB$P_pred

plot(X, type ="l",col = "black", lwd = 1, xlab="Time", main = "True vs Estimated data", sub = "results from original balanced df", ylab ="Values")
lines(X_filtered_B, col = "red")
lines(X_filtered, col = "green")



### 2. try estimating KF on PCA imputed balanced data

pca_result_missMDA <- prcomp(data_imputed_missMDA, scale. = TRUE) # Perform PCA on the imputed data
summary(pca_result_missMDA) # Summary of PCA results
#df_imputed <- data.frame(data_imputed_missMDA)
pca1 <- pca_result_missMDA$x[,1]
ar_model <- ar(pca1, order.max = 1, method = "mle") #fit AR(1) on PC1 scores
ar_pca <- ar_model$ar
initial_params <- c(ar_pca, 1.0, 1.0, -1.0, 1.0, 1.0, log(1), log(1), log(1), log(1), log(1), log(1))

Yi1 <- df_imputed$Y1
Yi2 <- df_imputed$Y2
Yi3 <- df_imputed$Y3
Yi4 <- df_imputed$Y4
Yi5 <- df_imputed$Y5

mle_results_I <- optim(initial_params, kalman_log_likelihood, Y1 = Yi1, Y2 = Yi2, Y3 = Yi3, Y4 = Yi4, Y5 = Yi5, hessian = TRUE, method = "BFGS")
optimized_params_I <- mle_results_I$par
final_resultI <- kalman_log_likelihood(optimized_params_I, Y1 = Yi1, Y2 = Yi2, Y3 = Yi3, Y4 = Yi4, Y5 = Yi5, return_all = TRUE)
# Extract the results
X_filtered_I <- final_resultI$X_filtered
P_filtered_I <- final_resultI$P_filtered
P_pred_I <- final_resultI$P_pred

plot(X, type ="l",col = "black", lwd = 1, xlab="Time", main = "Final Comparison", sub = "results from original and PCA balanced df", ylab ="Values")
lines(X_filtered_B, col = "red")
lines(X_filtered_I, col = "orange")
lines(X_filtered, col = "green")

forecast_error_filt_NA <- mean((X-X_filtered)^2)
forecast_error_filt_B <- mean((X-X_filtered_B)^2)
forecast_error_filt_I <- mean((X-X_filtered_I)^2)

cat("mean squared forecast error of filtered series with NA values:", forecast_error_filt_NA, "\n")
cat("mean squared forecast error of filtered series with original balanced df:", forecast_error_filt_B, "\n")
cat("mean squared forecast error of filtered series from PCA imputed df:", forecast_error_filt_I, "\n")


# in sostanza, le stime effettuate con NA, ORIGINAL e PCA_IMPUTED appaiono confrontabili...sebbene visivamente sembra che le stime
# con NA seguano meglio di tutte il latent state vero.
# provare con meno correlazione fra gli errori delle serie osservate!
