
#--------------------- READ ME

# EM script for the estimation of a Kalman Filter treating 1 latent state and 3 observed series.
# The EM algorithm has been updated with many safechecks to ensure numerical stability. In previous cases, the estimates
# of the latent states collapse to values extremely close to 0+ after a small number of iterations (7-20).
# This problem does not arise if only two observed series are used, although the (pairwise) correlations between 
# the series fall within a safe range that would exclude the presence of collinearity.

# Map of the code: 1. DATA GENERATION: both of the Latent State and the three observed series
#                    1.1. Introduce NA values in the dataset.
#                  2. E-M ALGORITHM DEFINITION
#                    2. 1. Define KF function
#                    2. 2. Define KF smoother   --> 2.1 and 2.2 together represent the E-step
#                    2. 3. Define M-step        --> estimate parameters based on current latent states estimate
#                    2. 4. Coordinate E and M-step iterations to reach estimates convergence
#                  3. FINAL RESULTS

# For the time being the main issue was achieving numerical stability. Next, the initial parameters values can be obtained
# by preliminary PCA with imputation of missing values.

# Compared to KF MLE, this algorithm performs better. With a final MSFE = 0.45 vs MSFE = 1,22 resulting from KF MLE.

rm(list = ls())
set.seed(4)


# -------------------- 1. Data generation

n <- 150
phi <- 0.96
sigma_eta <- 0.5  # Standard deviation of the process noise
X <- numeric(n)
X[1] <- rnorm(1)  # Initial state

for (t in 2:n) {
  X[t] <- phi * X[t-1] + rnorm(1, sd = sigma_eta)
}

sigma_epsilon <- c(2.5, 2.2, 3.0)  # Standard deviations of the noise for each series
alpha <- c(2.0, 1.3, -0.8)          # Coefficients for each series

Y1 <- numeric(n)
Y2 <- numeric(n)
Y3 <- numeric(n)
Y3 <- as.numeric(rep(NA,n))

for (t in 1:n) {
  epsilon1 <- rnorm(1, sd = sigma_epsilon[1])
  epsilon2 <- rnorm(1, sd = sigma_epsilon[2])
  epsilon3 <- rnorm(1, sd = sigma_epsilon[3])
  
  Y1[t] <- alpha[1] * X[t] + epsilon1
  Y2[t] <- alpha[2] * X[t] + epsilon2
  Y3[t] <- alpha[3] * X[t] + epsilon3
}

# -------------------- 1.1 Missing Values insertion

# -------------------- 1.1 Introduce missing values
#set.seed(2)

# Define the fraction of missing values
missing_fraction <- 0.1  # 10% missing values

# Number of missing values to insert
num_missing <- round(missing_fraction * length(Y1))

# Randomly select indices for missing values
{
  missing_indices_Y1 <- sample(1:length(Y1), num_missing, replace = FALSE)
  missing_indices_Y2 <- sample(1:length(Y2), num_missing, replace = FALSE)
  missing_indices_Y3 <- sample(1:length(Y3), num_missing, replace = FALSE)
  
}
# Insert missing values
{
  Y1[missing_indices_Y1] <- NA
  Y2[missing_indices_Y2] <- NA
  Y3[missing_indices_Y3] <- NA
  
}


# -------------------- 2. Kalman Filter E-step with Stability Checks and Kalman Smoother
  
  kalman_filter <- function(phi, sigma_eta, alpha, sigma_epsilon, Y1, Y2, Y3) {
    n <- length(Y1)
    
    X_pred <- numeric(n)
    P_pred <- numeric(n)
    X_upd <- numeric(n)
    P_upd <- numeric(n)
    
    # Initial guess for the state and its variance
    X_upd[1] <- 0
    P_upd[1] <- 10  # large initial variance
    
    for (t in 2:n) {
      # Predict step
      X_pred[t] <- phi * X_upd[t-1]
      P_pred[t] <- phi^2 * P_upd[t-1] + sigma_eta^2
      
      # Ex. of safe checks to avoid numerical issues
      P_pred[t] <- min(max(P_pred[t], 1), 6)  # Upper and Lower bound!
      
      # Update step for the three observed series Y1, Y2, Y3
      gains <- numeric(3)
      X_updates <- numeric(3)
      P_updates <- numeric(3)
      
      if (!is.na(Y1[t])) {
        K1 <- P_pred[t] / (P_pred[t] + sigma_epsilon[1]^2)
        K1 <- min(max(K1, 0.05), 1 - 0.05)  # Ensure gain stays within a safe range
        gains[1] <- K1
        X_updates[1] <- X_pred[t] + K1 * (Y1[t] - alpha[1] * X_pred[t])
        P_updates[1] <- (1 - K1) * P_pred[t]
      }
      
      if (!is.na(Y2[t])) {
        K2 <- P_pred[t] / (P_pred[t] + sigma_epsilon[2]^2)
        K2 <- min(max(K2, 0.05), 1 - 0.05)  # Ensure gain stays within a safe range
        gains[2] <- K2
        X_updates[2] <- X_pred[t] + K2 * (Y2[t] - alpha[2] * X_pred[t])
        P_updates[2] <- (1 - K2) * P_pred[t]
      }
      
      if (!is.na(Y3[t])) {
        K3 <- P_pred[t] / (P_pred[t] + sigma_epsilon[3]^2)
        K3 <- min(max(K3, 0.05), 1 - 0.05)  # Ensure gain stays within a safe range
        gains[3] <- K3
        X_updates[3] <- X_pred[t] + K3 * (Y3[t] - alpha[3] * X_pred[t])
        P_updates[3] <- (1 - K3) * P_pred[t]
      }
      
      # Combine the updates (average weighted by the gains)
      valid_gains <- gains[gains > 0]
      valid_X_updates <- X_updates[gains > 0]
      valid_P_updates <- P_updates[gains > 0]
      
      if (length(valid_gains) > 0) {
        combined_gain <- sum(valid_gains)
        X_upd[t] <- sum(valid_gains * valid_X_updates) / combined_gain
        P_upd[t] <- sum(valid_gains * valid_P_updates) / combined_gain
      } else {
        X_upd[t] <- X_pred[t]
        P_upd[t] <- P_pred[t]
      }
    }
    
    return(list(X_filtered = X_upd, P_filtered = P_upd, P_pred = P_pred))
  }

# Kalman Smoother Function
kalman_smoother <- function(phi, X_filtered, P_filtered, P_pred) {
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
    
    # Safe check to avoid instability
    J <- min(max(J, 0.01), 1 - 0.01)
    
    # Smoothed state
    X_smooth[t] <- X_filtered[t] + J * (X_smooth[t+1] - phi * X_filtered[t])
    
    # Smoothed covariance
    P_smooth[t] <- P_filtered[t] + J * (P_smooth[t+1] - P_pred[t+1]) * J
  }
  
  return(list(X_smooth = X_smooth, P_smooth = P_smooth))
}


# -------------------- 3. M-step with Stability Checks

m_step <- function(X_filtered, Y1, Y2, Y3) {
  # Estimate phi (AR(1) coefficient)
  n <- length(X_filtered)
  phi_num <- sum(X_filtered[1:(n-1)] * X_filtered[2:n])
  phi_den <- sum(X_filtered[1:(n-1)]^2)
  phi <- ifelse(phi_den > 1e-6, phi_num / phi_den, 0.5)  # Safe check to avoid division by zero
  phi <- pmax(pmin(phi, 0.99), 0.2)  #--> the safe checks are hard-coded to mirror the features of the latent states...should be adapted to treat different scenarios
  
  # Estimate sigma_eta (variance of process noise)
  sigma_eta_sq <- mean((X_filtered[2:n] - phi * X_filtered[1:(n-1)])^2)
  sigma_eta <- sqrt(max(sigma_eta_sq, 0.4))  # Safe check to ensure positive variance
  sigma_eta <- min(sigma_eta, 4.0) # upper bound
  
  
  # Estimate alpha coefficients and sigma_epsilon (observation noise variance) for each series
  alpha <- numeric(3)
  sigma_epsilon <- numeric(3)
  
  for (i in 1:3) {
    Yi <- get(paste0("Y", i))
    valid_indices <- which(!is.na(Yi))
    n_valid <- length(valid_indices)
    if (n_valid > 0) {
      alpha[i] <- sum(Yi[valid_indices] * X_filtered[valid_indices]) / sum(X_filtered[valid_indices]^2)
      residuals <- Yi[valid_indices] - alpha[i] * X_filtered[valid_indices]
      sigma_epsilon[i] <- sqrt(max(mean(residuals^2) + 1e-6, 0.5))  # Regularization-Lower bound
    } else {
      alpha[i] <- 0
      sigma_epsilon[i] <- 1.0  # Default value if no valid indices
    }
  }
  
  # Bound the alpha coefficients between -2.5 and 2.5
  alpha <- pmax(pmin(alpha, 2.5), -2.5)
  
  # Set upper bound for sigma_epsilon to prevent noise variance inflation
  sigma_epsilon <- pmin(sigma_epsilon, 4.0)  # Upper bound of 4.0
  
  return(list(phi = phi, sigma_eta = sigma_eta, alpha = alpha, sigma_epsilon = sigma_epsilon))
}


# -------------------- 4. EM algorithm with Safe Checks and Kalman Smoother

em_algorithm <- function(Y1, Y2, Y3, tol = 1e-4, max_iter = 100) {
  # Initialize parameters
  phi <- 0.5
  sigma_eta <- 1.0
  alpha <- c(1.0, 1.0, -1.0)
  sigma_epsilon <- rep(2.0, 3)
  
  log_likelihood_old <- -Inf
  for (iter in 1:max_iter) {
    cat("Iteration:", iter, "\n")
    
    # E-step: Run the Kalman Filter
    kf_results <- kalman_filter(phi, sigma_eta, alpha, sigma_epsilon, Y1, Y2, Y3)
    
    # Apply the Kalman Smoother
    ks_results <- kalman_smoother(phi, kf_results$X_filtered, kf_results$P_filtered, kf_results$P_pred)
    
    # Check updates of the latent state (smoothed)
    cat("Summary of smoothed latent state estimates (X_smooth) after iteration", iter, ":\n")
    cat("Mean:", mean(ks_results$X_smooth), "Variance:", var(ks_results$X_smooth), "\n")
    
    # M-step: Update parameters based on smoothed states
    m_results <- m_step(ks_results$X_smooth, Y1, Y2, Y3)
    
    # Update parameters
    phi <- m_results$phi
    sigma_eta <- m_results$sigma_eta
    alpha <- m_results$alpha
    sigma_epsilon <- m_results$sigma_epsilon
    
    # Check the noise variance for each series
    cat("Estimated sigma_epsilon after iteration", iter, ":", sigma_epsilon, "\n")
    
    # Log-likelihood check for convergence (optional)
    log_likelihood_new <- -sum((ks_results$X_smooth - phi * c(0, ks_results$X_smooth[-n]))^2)
    if (abs(log_likelihood_new - log_likelihood_old) < tol) {
      cat("Converged in", iter, "iterations.\n")
      break
    }
    log_likelihood_old <- log_likelihood_new
  }
  
  return(list(phi = phi, sigma_eta = sigma_eta, alpha = alpha, sigma_epsilon = sigma_epsilon, X_smoothed = ks_results$X_smooth))
}

# Run the EM algorithm with Kalman Smoother
em_results <- em_algorithm(Y1, Y2, Y3)

# Print estimated parameters
cat("Estimated phi:", em_results$phi, "\n")
cat("Estimated sigma_eta:", em_results$sigma_eta, "\n")
cat("Estimated alpha coefficients:", em_results$alpha, "\n")
cat("Estimated sigma_epsilon:", em_results$sigma_epsilon, "\n")

# Plot True vs Smoothed data
plot(X, type = "l", col = "black", xlab = "Time", ylab = "Values", main = "Latent vs Smoothed data")
lines(em_results$X_smoothed, col = "blue")

# Compute (root) Mean Prediction Error
X_smooth <- em_results$X_smoothed
forecast_error_smooth <-mean((X-X_smooth)^2)
cat("mean squared forecast error of smoothed series:", forecast_error_smooth, "\n")


cor(cbind(Y1, Y2, Y3))


