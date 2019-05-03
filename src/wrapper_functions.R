################################
#
# Wrapper functions 
# 
# AWalker July 2018
#
################################



# MCMC functions
################################

# set parameter boundaries from prior distributions
boundary_handling_set <- function(.) {
  # number of samples to be used in boundary handling
  n <- 1000
  boundary_sample <- lapply(.$dynamic$pars_eval, function(cs) eval(parse(text=cs)) )
  .$mcmc$boundary_min <- unlist(lapply(boundary_sample,min))
  .$mcmc$boundary_max <- unlist(lapply(boundary_sample,max))
  rm(boundary_sample) 
}

# restrict parameter proposals that are beyond the boundary
boundary_handling <- function(., ii, jj ) {
  if      (.$dataf$pars[ii,jj] < .$mcmc$boundary_min[jj]) .$dataf$pars[ii,jj] <- .$mcmc$boundary_min[jj]
  else if (.$dataf$pars[ii,jj] > .$mcmc$boundary_max[jj]) .$dataf$pars[ii,jj] <- .$mcmc$boundary_max[jj]
}   
   
   
# generate proposal using DE-MC algorithm  
proposal_generate_demc <- function(.,j) {

  # number of data points to be used in boundary handling
  # this doesn't need to be done every iteration - can be done at the beginning of the MCMC and maxes and mins stored
  #n <- 1000
  #.$dynamic$pars_bndhndling <- lapply(.$dynamic$pars_eval, function(cs) eval(parse(text=cs)) )
  #minn <- unlist(lapply(.$dynamic$pars_bndhndling,min))
  #maxn <- unlist(lapply(.$dynamic$pars_bndhndling,max))
   
  # Metropolis sampling
  # scaling factor
  d <- ncol(.$dataf$pars)  
  gamma_star <- 2.38 / sqrt(d + d)
  # b-value should be small compared to width of target distribution
  # what is the name of this b parameter? 
  b  <- 0.01
  R1 <- 0
  R2 <- 0
  # evaluate for each chain
  for (ii in 1:.$dataf$lp) {
    # randomly select two different numbers R1 and R2 unequal to j
    # from a uniform distribution without replacement
    while ((R1 == 0) | (R1 == ii))               R1 <- ceiling(runif(1,min=0,max=1)*.$dataf$lp)
    while ((R2 == 0) | (R2 == ii) | (R2 == R1))  R2 <- ceiling(runif(1,min=0,max=1)*.$dataf$lp)
    # draw random number from uniform distribution on interval (-b,b)
    uniform_r <- runif(1,min=(-b),max=b)
    # evaluate for each parameter value
    for (jj in 1:d) {
      # generate proposal via Differential Evolution
      .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j] + gamma_star * (.$dataf$pars_array[R1,jj,j] - .$dataf$pars_array[R2,jj,j]) + uniform_r

      # call boundary handling function here
      boundary_handling(., ii, jj ) 
#      # boundary handling for minumum
#      if (.$dataf$pars[ii,jj] < minn[jj]) {
#        .$dataf$pars[ii,jj] <- minn[jj]
#      }
#      # boundary handling for maximum
#      if (.$dataf$pars[ii,jj] > maxn[jj]) {
#        .$dataf$pars[ii,jj] <- maxn[jj]
#      } 
    }
  }
} 


# calculate proposal acceptance using the Metropolis ratio (for DE-MC algorithm)
proposal_accept_demc <- function(., j, lklihood ) {

  # Metropolis ratio 
  metrop_ratio <- exp(lklihood - .$dataf$pars_lklihood[ ,j])
  alpha        <- pmin(1,metrop_ratio)
  for(kk in 1:.$dataf$lp) {
    # accept if Metropolis ratio > random number from uniform distribution on interval (0,1) 
    accept <- log(alpha[kk]) > log(runif(1,min=0,max=1))
    .$dataf$pars_array[kk,,j+1]   <- if(accept) .$dataf$pars[kk,] else .$dataf$pars_array[kk,,j]    
    .$dataf$pars_lklihood[kk,j+1] <- if(accept) lklihood[kk]      else .$dataf$pars_lklihood[kk,j]    
    #print(c(kk, accept))
    #print(c(.$dataf$pars_array[kk,,j+1], .$dataf$pars[kk,], .$dataf$pars_array[kk,,j] ))   
    out_n <- .$wpars$mcmc_maxiter/2       
    if (j > out_n) 
      .$dataf$out_mcmc[kk,,(j-out_n)] <- if(accept | j==out_n+1) .$dataf$out[kk,] else .$dataf$out_mcmc[kk,,(j-out_n-1)]    
#APW: check this once decided on j indexing - is j correctly specified given the many j+1 s?
#    if ((j+1) > out_n) 
#      .$dataf$out_mcmc[kk,,(j+1-out_n)] <- if(accept | j==out_n+1) .$dataf$out[kk,] else .$dataf$out_mcmc[kk,,(j-out_n-1)]    
  }
}   
   
   
######################################################################################################################################################

# static part of DREAM algorithm 
static_dream <- function(.) {
    
  # number of parameters being estimated
  .$mcmc$d <- ncol(.$dataf$pars)

  # preallocate memory space for algorithmic variables
  .$mcmc$J             <- numeric(.$mcmc$n_CR)
  .$mcmc$n_id          <- numeric(.$mcmc$n_CR)
  .$mcmc$CR            <- numeric(.$mcmc$n_CR)
  .$mcmc$p_CR          <- numeric(.$mcmc$n_CR)
  .$mcmc$R             <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$dataf$lp - 1)
  .$mcmc$current_state <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$mcmc$d)
  .$mcmc$p_state       <- numeric(.$dataf$lp)
  .$mcmc$std_state     <- numeric(.$mcmc$d)
  .$mcmc$jump          <- matrix(data=0,nrow=.$dataf$lp,ncol=.$mcmc$d)
  .$mcmc$draw          <- matrix(data=0,nrow=.$dataf$lp-1,ncol=.$dataf$lp)
  .$mcmc$lambda        <- matrix(data=0,nrow=.$dataf$lp,ncol=1)
  
  # index of chains for Differential Evolution
  for (kk in 1:.$dataf$lp) .$mcmc$R[kk, ] <- setdiff(1:.$dataf$lp, kk)

  # crossover values
  .$mcmc$CR[] <- 1:.$mcmc$n_CR / .$mcmc$n_CR

  # selection probability of crossover values
  .$mcmc$p_CR[] <- rep(1, .$mcmc$n_CR) / .$mcmc$n_CR
 
  # vector that stores how many times crossover value indices are used
  # initialized to 1's in order to avoid numeric issues
  .$mcmc$n_id[] <- rep(1, .$mcmc$n_CR)
  # originally initialized to 0's in Vrugt's algorithm
  #.$mcmc$n_id[] <- rep(0, .$mcmc$n_CR)
}

# generate proposal using DREAM algorithm
proposal_generate_dream <- function(., j) {
  
  # reset matrix of jump vectors to zero
  .$mcmc$jump[] <- matrix(data = 0)
  # .$mcmc$jump[] <- 0 # APW: I think you can just do this, the square brackets preserve the data structure 
  
  # current state ('mcmc_chains' number of samples of a d-variate distribution)
  .$mcmc$current_state[] <- matrix(.$dataf$pars_array[ , , j], nrow = .$dataf$lp, ncol = .$mcmc$d)

  # NOTE, this chunk is a repeat of DE-MC code
  # boundary handling
  # number of data points to be used in boundary handling
  #n <- 1000
  #.$dynamic$pars_bndhndling <- lapply(.$dynamic$pars_eval, function(cs) eval(parse(text = cs)) )
  #minn <- unlist(lapply(.$dynamic$pars_bndhndling,min))
  #maxn <- unlist(lapply(.$dynamic$pars_bndhndling,max))

  # permute [1,2,...,mcmc_chains-1] mcmc_chains number of times
  .$mcmc$draw[] <- apply(matrix(runif((.$dataf$lp - 1) * .$dataf$lp), .$dataf$lp - 1, .$dataf$lp), 2, function(v) sort(v, index.return = T)$ix)
  
  # create a .$dataf$lp x 1 matrix of continuous uniform random values between -c_rand and c_rand
  .$mcmc$lambda[] <- matrix(runif(.$dataf$lp * 1, -.$mcmc$c_rand, .$mcmc$c_rand), .$dataf$lp)

  # compute standard deviation of each dimension (ie,compute standard deviation of each column of current_state matrix)
  .$mcmc$std_state[] <- apply(.$mcmc$current_state, 2, sd)

  # NOTE vectorize this inner for-loop to improve computational efficiency, but this is non-trivial     
  # create proposals
  for (ii in 1:.$dataf$lp) {
    
    # select delta (equal selection probability) (ie, choose 1 value from the vector [1:delta] with replacement)
    D <- sample(1:.$mcmc$delta,1,replace=T)
 
    # extract vectors a and b not equal to ii
    a <- .$mcmc$R[ii, .$mcmc$draw[1:D, ii]]
    b <- .$mcmc$R[ii, .$mcmc$draw[(D + 1):(2 * D), ii]]
 
    # select index of crossover value (weighted sample with replacement)
    .$mcmc$id <- sample(1:.$mcmc$n_CR, 1, replace = T, prob = .$mcmc$p_CR)
 
    # draw d values from uniform distribution between 0 and 1
    zz <- runif(.$mcmc$d)
 
    # derive subset A of selected dimensions
    A <- which(zz < .$mcmc$CR[.$mcmc$id])

    #  how many dimensions are sampled
    d_star <- length(A)
 
    # make sure that A contains at least one value
    if (d_star == 0) A <- which.min(zz); d_star <- 1
 
    # calculate jump rate
    gamma_d <- 2.38 / sqrt(2 * D * d_star)

    # NOTE maybe there is a way to consolidate this  chunk of code more efficiently      
    # select gamma
    temp1 <- c(gamma_d, 1)
    temp2 <- c(1 - .$mcmc$p_gamma, .$mcmc$p_gamma)
    gamma <- sample(temp1, 1, replace = T, prob = temp2)

    # compute jump differential evolution of ii-th chain
    .$mcmc$jump[ii, A] <- .$mcmc$c_ergod * rnorm(d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[a, A] - .$mcmc$current_state[b, A]), dim = 1)
 
    # compute proposal of ii-th chain
    .$dataf$pars[ii,1:.$mcmc$d] <- .$mcmc$current_state[ii, 1:.$mcmc$d] + .$mcmc$jump[ii, 1:.$mcmc$d]

    # call boundary handling function here
    boundary_handling(., ii, jj ) 
    # NOTE, this chunk is a repeat of DE-MC code
    # APW: replace with boundary handling function 
    # more boundardy handling
    #for (jj in 1:.$mcmc$d) {
    #  # boundary handling for minumum
    #  if (.$dataf$pars[ii, jj] < minn[jj]) {
    #    .$dataf$pars[ii, jj] <- minn[jj]
    #  } 
    #  # boundary handling for maximum
    #  if (.$dataf$pars[ii, jj] > maxn[jj]) {
    #    .$dataf$pars[ii, jj] <- maxn[jj]
    #  } 
    #}
  }
} 

# proposal acceptance function for the DREAM algorithm
# in the future: could probably consolidate this with the acceptance function for the DE-MC algorithm
proposal_accept_dream <- function(., j, lklihood) {

  # likelihood of current state
  .$mcmc$p_state[] <- .$dataf$pars_lklihood[ ,j] 

  for (qq in 1:.$dataf$lp) {

    # compute Metropolis acceptance probability
    alpha <- min(1, exp(lklihood[qq] - .$mcmc$p_state[qq]))

    # in the future: figure out how to make this clunky block of code prettier      
    # determine if p_acc is larger than random number drawn from uniform distribution on interval [0,1]
    if (alpha > runif(1, min = 0, max = 1)) {
      # if true, accept the proposal
      accept <- TRUE
      .$mcmc$current_state[qq, 1:.$mcmc$d] <- .$dataf$pars[qq, 1:.$mcmc$d]
      .$mcmc$p_state[qq] <- lklihood[qq]
      # append accepted current_state and probability density to storage data frames
      .$dataf$pars_array[qq, 1:.$mcmc$d, j + 1] <- .$mcmc$current_state[qq, 1:.$mcmc$d]
      .$dataf$pars_lklihood[qq, j + 1] <- .$mcmc$p_state[qq]
	  # store generated proposal (regardless of whether accepted or not)
	  .$dataf$prop_storage[qq, 1:.$mcmc$d, j + 1] <- .$dataf$pars[qq, 1:.$mcmc$d]
    } else {
      accept <- FALSE
      # set jump back to zero for p_CR
      .$mcmc$jump[qq, 1:.$mcmc$d] <- 0
      # repeat previous current_state and probability density in storage data frames
      .$dataf$pars_array[qq, 1:.$mcmc$d, j + 1] <- .$dataf$pars_array[qq, 1:.$mcmc$d, j]
      .$dataf$pars_lklihood[qq, j + 1] <- .$dataf$pars_lklihood[qq, j]
	  # store generated proposal (regardless of whether accepted or not)
	  .$dataf$prop_storage[qq, 1:.$mcmc$d, j + 1] <- .$dataf$pars[qq, 1:.$mcmc$d]
    }
 
    # update jump distance crossover index
    .$mcmc$J[.$mcmc$id] <- .$mcmc$J[.$mcmc$id] + sum((.$mcmc$jump[qq, 1:.$mcmc$d] / .$mcmc$std_state)^2)

    # number of times index crossover is used
    .$mcmc$n_id[.$mcmc$id] <- .$mcmc$n_id[.$mcmc$id] + 1

    # NOTE this chunck is a reapeat of DE-MC code
    out_n <- .$wpars$mcmc_maxiter / 2
    if (j > out_n) 
      .$dataf$out_mcmc[qq, ,(j - out_n)] <- if(accept | j == out_n + 1) .$dataf$out[qq, ] else .$dataf$out_mcmc[qq, , (j - out_n - 1)]  

  }

  # update selection probability of crossover
  # altered original algorithm here to account for numerical issues
  if ((j < (.$wpars$mcmc_maxiter / 10)) & (sum(.$mcmc$J) > 0)) {
    .$mcmc$p_CR <- .$mcmc$J / .$mcmc$n_id
    .$mcmc$p_CR <- .$mcmc$p_CR / sum(.$mcmc$p_CR)
  }

}
 
# function for detection and correction of outlier chains
# outlier_check <- function(.) {
# }

# test for convergence using the R-statistic convergence diagnostic of Gelman and Rubin  
# chain_converge <- function(.) {
# }    

######################################################################################################################################################


# Likelihood functions
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
