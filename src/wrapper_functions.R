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
 


### END ###
