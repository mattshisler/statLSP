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
avg_model_sigma <- sigma(avg_model_over_years)

# construct the basis functions from the 1st-order Taylor approx
B <- basis_functions(1:365, avg_model_theta)
matplot(B, type = "l")


# MCMC for the hierarchical model (IN-PROGRESS)

# priors
Sigma <- 100*diag(npar)
Q     <- solve(Sigma)



DL_m  <- double_logistic(floor(sim$data$year1$t), avg_model_theta)
tau   <- 1/sigma^2
X     <- B[floor(sim$data$year1$t),]
y     <- sim$data$year1$y_noise1

# sample basis function coefficients
cov_b     <- solve(tau*t(X)%*%X + Q)
mn_b      <- tau*cov_b%*%t(X)%*%(y-DL_m)
samps_b   <- t(chol(cov_b))%*%rnorm(npar)

# compute theta sample
theta_hat <- avg_model_theta + mn_b + samps_b


# sample from theta population distribution


# sample from 




r <- sim$data$year1$y_noise1-double_logistic(sim$data$year1$t,avg_model_theta)
X <- B[floor(sim$data$year1$t),]


b <- lm(r ~ X-1)
b$coefficients + avg_model_theta
sim$theta[1,]







