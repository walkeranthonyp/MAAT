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
  # read in measurement error
  #meas_sigma <- .$dataf$obsse

  # remove zeros from measurement error  
  #subset1    <- sapply(meas_sigma, function(row) all(row > 0))
  # change this!
  #meas_sigma2 <- meas_sigma[subset1]
  sspos <- which(.$dataf$obsse > 1e-9)

  # number of measured data points (that do not have zero uncertainty)
  # obs_n <- length(.$dataf$obs)
  obs_n <- length(sspos)
  # calculate error residual
  # each chain is on rows of dataf$out, take transpose
  #error_residual_matrix <- t(.$dataf$out)[sspos,] - .$dataf$obs[sspos]
  # subset error_residual_matrix
  #error_residual_matrix2 <- error_residual_matrix[subset1, ]

  # derive the log density
  if (length(.$dataf$obsse) == 1) {
    # homoscedastic error - I think this is an assumption that the error is constant across all values of obs
    # calculate sum of squared error
    error_residual_matrix <- t(.$dataf$out)[sspos,] - .$dataf$obs[sspos]
    SSR                   <- apply(error_residual_matrix, 2, function(v) sum(v^2))
    log_density           <- -(obs_n/2)*log(2*pi) - obs_n*log(abs(.$dataf$obsse)) - 0.5*.$dataf$obsse^(-2)*SSR
  } else {
    # heteroscedastic error
    #err_div_sig     <- error_residual_matrix / as.vector(.$dataf$obsse[sspos])
    error_residual_matrix <- ( t(.$dataf$out)[sspos,] - .$dataf$obs[sspos] ) / .$dataf$obsse[sspos]
    #err_div_sig     <- error_residual_matrix / .$dataf$obsse[sspos]
    #sum_err_div_sig <- apply(err_div_sig, 2, function(v) sum(v^2))
    SSR                   <- apply(error_residual_matrix, 2, function(v) sum(v^2))
    # not sure if it is correct to take absolute value here??? - these are some metric of obs uncertainty so should all be > 0
    log_density           <- -(obs_n/2)*log(2*pi) - sum(log(abs(.$dataf$obsse[sspos]))) - 0.5*SSR
  }  

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  print(log_density)
  log_density
} 



### END ###
