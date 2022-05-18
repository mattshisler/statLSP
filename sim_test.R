# this is a sim test

rm(list=ls())

# load utility functions
source("utility_functions.R")

# true mean and variance for simulation
sim_theta_mu6  <- c(-1.80, 0.75, 120, 8, 290, 8) 
sim_theta_var6 <- c( 0.01, 1e-3,  25, 1,  25, 1)
# sim_theta_mu7  <- c(-1.80, 0.75, 120, 8, 290, 8, -7.750) 
# sim_theta_var7 <- c( 0.01, 1e-3,  25, 1,  25, 1,  0.075)

npar    <- 6
sigma   <- 0.05
nyear   <- 20
min_obs <- 25
max_obs <- 25

# simulate the evi2 data
sim <- sim_evi2(sim_theta_mu6, sim_theta_var6, 
                sigma=sigma, nyear=nyear,
                min_obs = min_obs, max_obs = max_obs)

# plot the simulated evi2 curves
plot_evi2(sim, "parallel")
plot_evi2(sim, "serial")

# create aggregated data set
sim_agg <- do.call(rbind, sim$data)

# construct average model by combining data over years (frequentist non-lin)
# uses the nlsLM function of the minpack.lm package
avg_model_over_years <- avg_fit(sim_agg = sim_agg)
avg_model_theta <- coef(avg_model_over_years)
avg_model_vcov  <- vcov(avg_model_over_years)
avg_model_sigma <- sigma(avg_model_over_years)

# construct the basis functions from the 1st-order Taylor approx
B <- basis_functions(1:365, avg_model_theta)
matplot(B, type = "l")

# MCMC for the hierarchical model (IN-PROGRESS)

# priors
# Sigma <- diag(npar)
Sigma <- diag(sim_theta_var6)
# Sigma   <- 10*avg_model_vcov
Q     <- solve(Sigma)
tau   <- 1/sigma^2
T0    <- 1
D0    <- 0.1

year  <- 2
DL_m  <- double_logistic(floor(sim$data[[year]]$t), avg_model_theta)
X     <- B[floor(sim$data[[2]]$t),]
y     <- sim$data[[2]]$y_noise1-DL_m

iters <- 5000
keep_param <- matrix(NA, nrow=iters, ncol=npar+1)

for (i in 1:iters){
  # sample basis function coefficients
  V <- solve(tau*t(X)%*%X + Q)
  M <- tau*V%*%t(X)%*%(y-DL_m)
  b <- as.numeric(rmvnorm(1,M,V))
  
  # sample sigma
  resid <- (y-DL_m-X%*%b)
  T1    <- T0 + npar
  D1    <- D0 + t(resid)%*%resid
  tau <- rgamma(1,T1,D1)
    
  # compute theta sample
  theta_hat <- avg_model_theta + b
  
  # store samples
  keep_param[i,] <- c(theta_hat, sqrt(1/tau))
}

apply(keep_param, 2, function(x) quantile(x, c(0.975,0.025)))
colMeans(keep_param)
sim$theta[2,]

hist(keep_param[,1])









