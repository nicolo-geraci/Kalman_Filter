# Kalman_Filter
Simple scripts to estimate parameters of a linear State-Space Model with gaussian errors via Kalman Filter.


This repository contains a few scripts for R programming language to estimate the latent states and parameters of a State Space model with gaussian errors.
I wrote this code to better understand the Kalman Filter and Smoother in action. The code is not fully efficient as it does not use any matrix operation, however the scalar notation greatly facilitates the understanding of each passage in the algorithm.
In each script, the model assumes the following State Space structure:
 
# y_t = C * x_t     + eta_t,   eta_t~N(0,R)  (observation equation)
# x_t = A * x_{t-1} + eps_t,   eps_t~N(0,Q)  (transition equation)

where    y         = dxT matrix of T observations for d observed series
         R         = dxd observational error matrix
         Q         = kxk innovation error matrix, where k is the number of latent series assumed
         A         = transition matrix
         C         = measurement matrix

Note that only one latent series is assumed in these models (i.e. k=1). 

From the assumption of gaussian errors it is possible to derive a cumulative likelihood function. The scripts estimates the latent series and the parameters
in A, C, R and Q via Maximum Likelihood, either with direct optimization or implementing the Expectation -Maximization algorithm.
Codes are designed to handle any missing values in the observed series.



