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
  obs_n <- length(.$dataf$obs) 

  # calculate error residual 
  # - each chain is on rows of dataf$out take transpose
  error_residual_matrix <- t(.$dataf$out) - .$dataf$obs  

  # calculate sum of squared error
  # - error_residual_matrix now has chains on columns 
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2) )

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix 
  -(obs_n/2)*log(SSR)
}   


# standard error probability density function with iid error residuals
# this function incorporates measurement errors
f_proposal_lklihood_ssquared_se <- function(.) {
  # read in measurement error and remove zeros from measurement error  
  sspos <- which(.$dataf$obsse > 1e-9)

  # number of measured data points (that do not have zero uncertainty)
  obs_n <- length(sspos)

  # observed error
  obsse <- if(.$wpars$mcmc_homosced)   rep(mean(.$dataf$obsse[sspos]), obs_n)
           else                .$dataf$obsse[sspos]
 
  # calculate error residual (each chain is on rows of dataf$out, take transpose)
  error_residual_matrix <- ( t(.$dataf$out)[sspos,] - .$dataf$obs[sspos] ) / obsse
  
  # calculate sum of squared error
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2))
 
  # derive log density   
  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  -(obs_n/2)*log(2*pi) - sum(log(obsse)) - 0.5*SSR
} 



### END ###
