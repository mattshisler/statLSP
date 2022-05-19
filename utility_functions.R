# equivalent to plogis with m=0 & s=1.
expit   <- function(x){1/(1+exp(-x))}

# first derivative of the expit function. equivalent to dlogis with m=0 & s=1.
expit_p <- function(x){expit(x)*(1-expit(x))}

# single year double logistic function. Used to model land surface phenology.
double_logistic <- function(t,theta){
  t   <- t%%365
  
  # start-of-season and end-of-season curves
  out <- expit((t - theta[3])/theta[4]) - 
         expit((t - theta[5])/theta[6])
 
  # returns 6 or 7 parameter version depending on the length
  # of the theta vector.
  if(length(theta)==6){
    # w/o greendown parameter
    out <- expit(theta[1]) + theta[2]*out
  } else if (length(theta)==7){
    # w/  greendown parameter
    out <- expit(theta[1]) + (theta[2]-expit(theta[7])*t)*out
    
  }
  
  return(out)
}

# construct the basis functions for the model. These are created using the
# gradient of the double logistic function with respect to the theta vector.
basis_functions <- function(t,theta){
  t   <- t%%365
  
  dl0 <- double_logistic(t,theta)-expit(theta[1])
  
  # returns 6 or 7 parameter version depending on the length
  # of the theta vector
  if(length(theta)==6){
    B1  <-  expit_p(theta[1])
    B2  <-  dl0/theta[2]
    B3  <- -(theta[2]*exp((theta[3] + t)/theta[4]))/
            (theta[4]*(exp(theta[3]/theta[4]) + exp(t/theta[4]))^2)
    B4  <-  (theta[2]*(theta[3]-t) * exp((theta[3]+t)/theta[4]))/
            ((theta[4]^2)*(exp(theta[3]/theta[4]) + exp(t/theta[4]))^2)
    B5  <-  (theta[2]*exp((theta[5]+t)/theta[6]))/
            (theta[6]*(exp(theta[5]/theta[6])+exp(t/theta[6]))^2)
    B6  <- -(theta[2]*(theta[5]-t)*exp((theta[5]+t)/theta[6]))/
            ((theta[6]^2)*(exp(theta[5]/theta[6])+exp(t/theta[6]))^2)
    B   <-  cbind(B1,B2,B3,B4,B5,B6)
  } else if (length(theta)==7){
    a <- expit_p(theta[7])
    theta[7] <- expit(theta[7])
    
    B1 <- expit_p(theta[1])
    B2 <- (1/(1+exp((theta[3]-t)/theta[4])))-(1/(1+exp((theta[5]-t)/theta[6])))
    B3 <- (theta[7]*t - theta[2])/(2*theta[4]*cosh((theta[3]-t)/theta[4])+2*theta[4])
    B4 <- ((theta[3] - t)*(theta[2]-theta[7]*t)*cosh((theta[3]-t)/(2*theta[4]))^(-2))/(4*theta[4]^2)
    B5 <- (theta[2] - theta[7]*t)/(2*theta[6]*cosh((theta[5]-t)/theta[6])+2*theta[6])
    B6 <- ((theta[5] - t)*(theta[7]*t-theta[2])*cosh((theta[5]-t)/(2*theta[6]))^(-2))/(4*theta[6]^2)
    B7 <- -t*a*((1/(1+exp((theta[3]-t)/theta[4])))-(1/(1+exp((theta[5]-t)/theta[6]))))
    B <- cbind(B1, B2, B3, B4, B5, B6, B7)
  }
  
  return(B)
}

sim_evi2 <- function(sim_theta_mu, sim_theta_var, sigma=0.05, nyear=20,
                     min_obs=25, max_obs=25){
  # number of parameters
  npar <- length(sim_theta_mu)
  nnoiselvl <- length(sigma)
  # nyear   <- 20
  # min_obs <- 25
  # max_obs <- 25
  
  # sample for number of observations per year
  if (min_obs==max_obs){
    obs_samp <- rep(min_obs,2)
  } else {
    obs_samp <- min_obs:max_obs
  }
  nobs <- sample(x=obs_samp, size=nyear, replace = T)
  
  # initialize data objects
  theta <- matrix(NA, nrow=nyear, ncol=npar)
  data  <- vector(mode="list", length=nyear)
  names(data) <- paste0("year",1:nyear)
  
  # sample from evi2 curve (double logistic model)
  for (i in 1:nyear){
    # allocate matrix object
    temp     <- data.frame(matrix(NA, nrow=nobs[i],ncol=2+nnoiselvl))
    
    # draw theta vector for year
    theta[i,] <- rnorm(npar,sim_theta_mu, sqrt(sim_theta_var))
    
    # draw observations and true evi2 value
    temp[,1] <- runif(nobs[i],1,365)
    temp[,2] <- double_logistic(temp[,1], theta[i,])
    
    # for each level of noise specified in sigma, generate
    # a new column with Gaussian error
    for (k in 1:nnoiselvl){
      temp[,(k+2)] <- pmax(rnorm(nobs[i],temp[,2],sigma[k]),0)
    }
    colnames(temp) <- c("t", "y_true", paste0("y_noise", 1:nnoiselvl))
    
    data[[i]] <- temp
  }
  
  return(list(theta=theta, data=data))
}

# plot evi2 in two modes:
# "parallel" - Year over year
# "serial"   - 
plot_evi2 <- function(sim, mode="parallel", frac=0.1){
  
  samps <- sample(1:length(sim$data), floor(frac*length(sim$data)))
  
  nyear <- floor(frac*length(sim$data))
  sim$data  <- sim$data[samps]
  sim$theta <- sim$theta[samps,]
  

  if (mode=="parallel"){
    plot(NA,xlim=c(1,365),ylim=c(0,1.1),xlab="t(DOY)",ylab="EVI2")
    for(k in 1:nyear){  
      points(sim$data[[k]]$t,sim$data[[k]][[3]], pch=19,col=k)
      lines(1:365, double_logistic(1:365,sim$theta[k,]), lwd=2, col=k)
    }
  } else if (mode=="serial"){
    plot(NA,xlim=c(1,365*nyear),ylim=c(0,1.1),xlab="t(DOY)",ylab="EVI2")
    for(k in 1:nyear){  
      points(sim$data[[k]]$t+(k-1)*365,sim$data[[k]][[3]], pch=19,col=k)
      lines((1:365)+(k-1)*365, double_logistic(1:365,sim$theta[k,]), lwd=2, col=k)
    }
  }
  
}

avg_fit <- function(sim_agg, npar){
  require(minpack.lm)
  
  if (npar==6){
    # six parameter model
    model_equ <- as.formula("y ~ 1/(1+exp(-theta1)) + (theta2) *
                            ((1 / (1 + exp((theta3 - t) / theta4))) -
                            (1 / (1 + exp((theta5 - t) / theta6))))")
  
    # fit using non-linear least squares
    fit <- nlsLM(model_equ,
                 data = list(y = sim_agg$y_noise1,
                             t = sim_agg$t),
                 weights = rep(1, nrow(sim_agg)),
                 start = list(theta1 = -2,
                              theta2 = 1,
                              theta3 = 120,
                              theta4 = 6,
                              theta5 = 290,
                              theta6 = 8),
                 lower = c(-10, 0.1, 1, 0, 185, 0),
                 upper = c(10, 100, 185, 100, 370, 100))
  } else if (npar==7){
    model_equ <- as.formula("y ~ 1/(1+exp(-theta1)) + 
                            (theta2-(1/(1+exp(-theta7)))*t) *
                            ((1 / (1 + exp((theta3 - t) / theta4))) -
                            (1 / (1 + exp((theta5 - t) / theta6))))")
    
    fit <- nlsLM(model_equ,
                 data = list(y = sim_agg$y_noise1,
                             t = sim_agg$t),
                 weights = rep(1, nrow(sim_agg)),
                 start = list(theta1 = -2,
                              theta2 = 1,
                              theta3 = 120,
                              theta4 = 6,
                              theta5 = 290,
                              theta6 = 8,
                              theta7 = -7),
                 lower = c(-10, 0.1, 1, 0, 185, 0, -10),
                 upper = c(10, 100, 185, 100, 370, 100, 10))
    
  }
  return(fit)
}


