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





