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
nyear   <- 10
min_obs <- 5
max_obs <- 25

# simulate the evi2 data
sim <- sim_evi2(sim_theta_mu6, sim_theta_var6, 
                sigma=sigma, nyear=nyear,
                min_obs = min_obs, max_obs = max_obs)
summary(sim$theta)

# plot the simulated evi2 curves
plot_evi2(sim, "parallel")
# plot_evi2(sim, "serial")

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

theta_hat <- theta_lower <- theta_upper <- sim$theta

for (year in 1:nyear){

  DL_m  <- double_logistic(floor(sim$data[[year]]$t), avg_model_theta)
  X     <- B[floor(sim$data[[year]]$t),]
  y     <- sim$data[[year]]$y_noise1-DL_m

  V <- solve(tau*t(X)%*%X + Q)
  M <- tau*V%*%t(X)%*%(y)

  theta_hat[year,]   <- avg_model_theta + as.numeric(M)
  theta_lower[year,] <- qnorm(0.025, mean=theta_hat[year,], sd=sqrt(diag(V)))
  theta_upper[year,] <- qnorm(0.975, mean=theta_hat[year,], sd=sqrt(diag(V)))
}

coverage <- theta_lower<sim$theta & theta_upper>sim$theta
colSums(coverage)/nyear





plot(1:365, double_logistic(1:365, sim$theta[year,]), type="l")
lines(1:365, double_logistic(1:365, theta_hat[year,]))
lines(1:365, double_logistic(1:365, theta_upper[year,]))
lines(1:365, double_logistic(1:365, theta_lower[year,]))


iters <- 5000
keep_param <- matrix(NA, nrow=iters, ncol=npar+1)

DL_m  <- double_logistic(floor(sim$data[[year]]$t), avg_model_theta)
X     <- B[floor(sim$data[[year]]$t),]
y     <- sim$data[[year]]$y_noise1-DL_m

for (i in 1:iters){
  # sample basis function coefficients
  V <- solve(tau*t(X)%*%X + Q)
  M <- tau*V%*%t(X)%*%(y)
  b <- as.numeric(rmvnorm(1,M,V))
  
  # sample sigma
  resid <- (y-DL_m-X%*%b)
  T1    <- T0 + npar
  D1    <- D0 + t(resid)%*%resid
  tau <- rgamma(1,T1,D1)
  
  # store samples
  keep_param[i,] <- c(avg_model_theta + b, sqrt(1/tau))
}


sim$theta[year,]
colMeans(keep_param)
apply(keep_param,2,function(x) quantile(x,0.025))
apply(keep_param,2,function(x) quantile(x,0.975))

theta_hat[nyear,]
theta_lower[nyear,]
theta_upper[nyear,]



apply(keep_param, 2, function(x) quantile(x, c(0.975,0.5,0.025)))
colMeans(keep_param)
sim$theta[2,]

hist(keep_param[,1])

M+avg_model_theta







