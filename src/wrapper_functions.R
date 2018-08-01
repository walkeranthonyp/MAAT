################################
#
# Wrapper functions 
# 
# AWalker July 2018
#
################################



# MCMC function
################################

# calculate proposal likelihood 
# expects model output to be probability - as in the output from the mixture model
f_proposal_lklihood_log <- function(.) {
  log(.$dataf$out)
}
   
# standard error probability density function with iid error residuals
f_proposal_lklihood_ssquared <- function(.) {
  # number of measured data points
  measurement_num <- length(.$dataf$obs) 
  # calculate error residual 
  # - each chain is on rows of dataf$out take transpose
  error_residual_matrix <- t(.$dataf$out) - .$dataf$obs  
  # calculate sum of squared error
  # - error_residual_matrix now has chains on columns 
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2) )
  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix 
  -(measurement_num/2)*log(SSR)
}   

# standard error probability density function with iid error residuals
# this function incorporates measurement errors
f_proposal_lklihood_ssquared_se <- function(.) {
  # read in measurement error
  meas_sigma <- .$dataf$obsse
  # number of measured data points
  measurment_num <- length(.$dataf$obs)
  # calculate error residual
  # each chain is on rows of dataf$out, take transpose
  error_residual_matrix <- t(.$dataf$out) - .$dataf$obs
  # derive the log density
  if (dim(meas_sigma)[1] == 1) {
    # homoscedastic error
    # calculate sum of squared error
    SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2))
    log_density <- -(measurment_num/2)*log(2*pi) - measurment_num*log(meas_sigma) - (1/2)*meas_sigma^(-2)*SSR
  }
  else {
    # heteroscedastic error
    # calculate sum(log(sigma))
    sum_log_sigma <- apply(meas_sigma, 2, function(v) sum(log(v)))
    # calculate sum( (error/sigma)^2 )
    err_div_sig <- error_residual_matrix / meas_sigma
    sum_err_div_sig <- apply(error_div_sig, 2, function(v) sum(v^2))
    log_density <- -(measurment_num/2)*log(2*pi) - sum_log_sigma - (1/2)*sum_err_div_sig
  }  
  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  log_density
} 


### END ###
