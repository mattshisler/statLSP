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
    out <- expit(theta[1]) + (theta[2]-theta[7]*t)*out
    
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
    B1  <-  1
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
  }
  
  return(B)
}



