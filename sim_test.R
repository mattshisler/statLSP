# this is a sim test
library(ggplot2)

rm(list=ls())

# load utility functions
source("utility_functions.R")

# true mean and variance for simulation
# sim_theta_mu  <- c(-1.80, 0.75, 120, 8, 290, 8) 
# sim_theta_var <- c( 0.01, 1e-3,  40, 1,  40, 1)
sim_theta_mu  <- c(-1.80, 0.75, 120, 8, 290, 8, -7.750)
sim_theta_var <- c( 0.01, 1e-3,  25, 1,  25, 1,  0.075)

npar    <- length(sim_theta_mu)
sigma   <- c(0.05)
nyear   <- 1000
min_obs <- 10
max_obs <- 10

# simulate the evi2 data
sim <- sim_evi2(sim_theta_mu, sim_theta_var, 
                sigma=sigma, nyear=nyear,
                min_obs = min_obs, max_obs = max_obs)
summary(sim$theta)

# plot a sample of the simulated simulated evi2 curves
plot_evi2(sim, "parallel", frac=0.1)
# plot_evi2(sim, "serial")

# create aggregated data set from 10% sample of simulated years
sampled_years <- sample(nyear,floor(0.1*nyear))
sim_agg <- do.call(rbind, sim$data[sampled_years])

# construct average model by combining data over years (frequentist non-lin)
# uses the nlsLM function of the minpack.lm package
avg_model_over_years <- avg_fit(sim_agg = sim_agg, npar = npar)
avg_model_theta <- coef(avg_model_over_years)
# avg_model_vcov  <- vcov(avg_model_over_years)
# avg_model_sigma <- sigma(avg_model_over_years)

# construct the basis functions from the 1st-order Taylor approx
B <- basis_functions(1:365, avg_model_theta)
matplot(B, type = "l")
matplot(scale(B), type = "l")

# numeric gradient check
numB <- B
for(k in 1:365){
  numB[k,] <- numDeriv::grad(double_logistic, x = avg_model_theta, t=k, method = "simple")
}
table(abs(B-numB)<1e-5)

# priors 
# In this idealized case we assume to know the residual variance and theta
# covariance.
V   <- diag(sim_theta_var)
Q   <- solve(V)
tau <- 1/sigma^2

# store results
theta_hat <- theta_hat_lower <- theta_hat_upper <- sim$theta

for (year in 1:nyear){
  
  # construct response and covariate matrix for bayesian regression
  DL_m  <- double_logistic(floor(sim$data[[year]]$t), avg_model_theta)
  X     <- B[floor(sim$data[[year]]$t),]
  y     <- sim$data[[year]]$y_noise1-DL_m

  # compure posterior mean and covariance
  V <- solve(tau*t(X)%*%X + Q)
  M <- tau*V%*%t(X)%*%(y)

  # estimate theta and 95% credible interval
  theta_hat[year,]   <- avg_model_theta + as.numeric(M)
  theta_hat_lower[year,] <- qnorm(.025, mean=theta_hat[year,], sd=sqrt(diag(V)))
  theta_hat_upper[year,] <- qnorm(.975, mean=theta_hat[year,], sd=sqrt(diag(V)))
}

# check if true theta lies within the credible interval
coverage <- theta_hat_lower<sim$theta & theta_hat_upper>sim$theta
colSums(coverage)/nyear


#### EXPERIMENTAL ####

# data <- data.frame(x=1:nyear, y_est=theta_hat[,3], lower = theta_lower[,3], upper = theta_upper[,3], y_true=sim$theta[,3])
# 
# ggplot(data=data,aes(x,y_est)) +
#   geom_point(aes(x,y_est)) +
#   geom_point(aes(x,y_true), color="red") +
#   geom_errorbar(aes(ymin=lower, ymax=upper))  
# 
# 
# ## ignore ##  
# widths <- theta_upper[,5]-theta_lower[,5]
# plot(sapply(sim$data, function(x) nrow(x)),widths)
# 
# 
# plot(1:365, double_logistic(1:365, sim$theta[year,]), type="l")
# lines(1:365, double_logistic(1:365, theta_hat[year,]))
# lines(1:365, double_logistic(1:365, theta_upper[year,]))
# lines(1:365, double_logistic(1:365, theta_lower[year,]))
# 
# 
# iters <- 5000
# keep_param <- matrix(NA, nrow=iters, ncol=npar+1)
# 
# DL_m  <- double_logistic(floor(sim$data[[year]]$t), avg_model_theta)
# X     <- B[floor(sim$data[[year]]$t),]
# y     <- sim$data[[year]]$y_noise1-DL_m
# 
# for (i in 1:iters){
#   # sample basis function coefficients
#   V <- solve(tau*t(X)%*%X + Q)
#   M <- tau*V%*%t(X)%*%(y)
#   b <- as.numeric(rmvnorm(1,M,V))
#   
#   # sample sigma
#   resid <- (y-DL_m-X%*%b)
#   T1    <- T0 + npar
#   D1    <- D0 + t(resid)%*%resid
#   tau <- rgamma(1,T1,D1)
#   
#   # store samples
#   keep_param[i,] <- c(avg_model_theta + b, sqrt(1/tau))
# }
# 
# 
# sim$theta[year,]
# colMeans(keep_param)
# apply(keep_param,2,function(x) quantile(x,0.025))
# apply(keep_param,2,function(x) quantile(x,0.975))
# 
# theta_hat[nyear,]
# theta_lower[nyear,]
# theta_upper[nyear,]
# 
# 
# 
# apply(keep_param, 2, function(x) quantile(x, c(0.975,0.5,0.025)))
# colMeans(keep_param)
# sim$theta[2,]
# 
# hist(keep_param[,1])
# 
# M+avg_model_theta
# 
# 
# diag((sigma^2/(1/8))*solve(t(X)%*%X))




