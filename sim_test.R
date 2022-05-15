# this is a sim test

rm(list=ls())

# load utility functions
source("utility_functions.R")

sim_theta_mu  <- c(-1.80, 0.75, 120, 8, 290, 8, -7.750) 
sim_theta_var <- c( 0.01, 1e-3,  25, 1,  25, 1,  0.075)
sigma         <- 0.05
nyear         <- 20
min_obs       <- 25
max_obs       <- 25

# simulate the evi2 data
sim <- sim_evi2(sim_theta_mu, sim_theta_var, 
                sigma=sigma, nyear=nyear,
                min_obs = min_obs, max_obs = max_obs)

# plot the simulated evi2 curves
plot_evi2(sim, "parallel")
plot_evi2(sim, "serial")
