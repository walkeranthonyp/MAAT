################################
#
# Wrapper functions for MCMC runs
#
# AJohnson, AWalker July 2019
#
################################



# MCMC functions
################################

# set parameter boundaries from prior distributions
boundary_handling_set <- function(.) {
  # number of samples to be used in boundary handling
  n <- 1e4
  boundary_sample     <- lapply(.$dynamic$pars_eval, function(cs) eval(parse(text=cs)) )
  .$mcmc$boundary_min <- unlist(lapply(boundary_sample,min))
  .$mcmc$boundary_max <- unlist(lapply(boundary_sample,max))
  rm(boundary_sample)
}


# restrict parameter proposals that are beyond the boundary
boundary_handling <- function(., ii, jj ) {
  if      (.$dataf$pars[ii,jj] < .$mcmc$boundary_min[jj]) .$dataf$pars[ii,jj] <- .$mcmc$boundary_min[jj]
  else if (.$dataf$pars[ii,jj] > .$mcmc$boundary_max[jj]) .$dataf$pars[ii,jj] <- .$mcmc$boundary_max[jj]
}



# DEMC functions
################################

# generate proposal using DE-MC algorithm
proposal_generate_mcmc_demc <- function(., j ) {

  # scaling factor
  d          <- ncol(.$dataf$pars)
  gamma_star <- 2.38 / sqrt(d + d)

  # b-value should be small compared to width of target distribution
  # b_rand specifies range for drawn "randomization" value
  b_rand  <- 0.01

  uniform_r <- runif(1, min = (-b_rand), max = b_rand)

  # evaluate for each chain
  for (ii in 1:.$dataf$lp) {

    # randomly select two different numbers R1 and R2 unequal to j, from uniform distribution without replacement
    R1 <- 0
    R2 <- 0
    while ((R1 == 0) | (R1 == ii))               R1 <- ceiling(runif(1, min = 0, max = 1) * .$dataf$lp)
    while ((R2 == 0) | (R2 == ii) | (R2 == R1))  R2 <- ceiling(runif(1, min = 0, max = 1) * .$dataf$lp)

    # evaluate for each parameter
    for (jj in 1:d) {

      # generate proposal via Differential Evolution
      .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r

      # call boundary handling function
      .$boundary_handling(ii = ii, jj = jj)

    }
  }
}


# calculate proposal acceptance using the Metropolis ratio (for DE-MC algorithm)
proposal_accept_mcmc_demc <- function(., j, lklihood ) {

  # Metropolis ratio
  metrop_ratio <- exp(lklihood - .$dataf$pars_lklihood[ ,j-1])

  alpha        <- pmin(1, metrop_ratio)

  # evaluate for each chain
  for(ii in 1:.$dataf$lp) {

    # accept if Metropolis ratio > random number from uniform distribution on interval (0,1)
    accept <- log(alpha[ii]) > log(runif(1, min = 0, max = 1))

    .$dataf$pars_array[ii,,j]   <- if(accept) .$dataf$pars[ii,] else .$dataf$pars_array[ii,,j-1]
    .$dataf$pars_lklihood[ii,j] <- if(accept) lklihood[ii]      else .$dataf$pars_lklihood[ii,j-1]

    # debug: store every model evaluation
    # out_n <- .$wpars$mcmc_maxiter / 2
    out_n <- 0
    if (j > out_n)
      .$dataf$out_mcmc[ii,,(j-out_n)] <- if(accept | j==out_n+1) .$dataf$out[ii,] else .$dataf$out_mcmc[ii,,(j-out_n-1)]
  }
}



# DREAM functions
################################

init_mcmc_demc <- function(.) NULL

# initialisation of DREAM algorithm
init_mcmc_dream <- function(.) {

  # number of parameters being estimated
  .$mcmc$d <- ncol(.$dataf$pars)

  # preallocate memory space for algorithmic variables
  .$mcmc$J             <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$n_id          <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$CR            <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$p_CR          <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$R             <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$dataf$lp - 1)
  .$mcmc$current_state <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$mcmc$d)
  .$mcmc$p_state       <- numeric(.$dataf$lp)
  .$mcmc$sd_state      <- numeric(.$mcmc$d)
  .$mcmc$jump          <- matrix(data=0, nrow=.$dataf$lp,   ncol=.$mcmc$d )
  .$mcmc$draw          <- matrix(data=0, nrow=.$dataf$lp - 1, ncol=.$dataf$lp )
  .$mcmc$lambda        <- matrix(data=0, nrow=.$dataf$lp,   ncol=1 )

  if(.$wpars$mcmc_debug) {
    .$mcmc$runif_seed     <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$wpars$mcmc_maxiter)
    .$mcmc$draw_seed      <- array(data = 0, dim = c(.$dataf$lp - 1, .$dataf$lp, .$wpars$mcmc_maxiter))
    .$mcmc$lambda_seed    <- array(data = 0, dim = c(.$dataf$lp, 1, .$wpars$mcmc_maxiter))
    .$mcmc$zz_seed        <- array(data = 0, dim = c(.$dataf$lp, .$mcmc$d, .$wpars$mcmc_maxiter))
    .$mcmc$prop_storage   <- array(data = 0, dim = c(dim(.$dataf$pars), .$wpars$mcmc_maxiter))
    .$mcmc$accept_storage <- array(data = 0, dim = c(.$dataf$lp, 1, .$wpars$mcmc_maxiter))
    .$mcmc$lklhd_storage  <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$wpars$mcmc_maxiter)
  }

  # index of chains for Differential Evolution
  for (ii in 1:.$dataf$lp) .$mcmc$R[ii, ] <- setdiff(1:.$dataf$lp, ii)

  # crossover values
  .$mcmc$CR[] <- 1:.$wpars$mcmc_n_CR / .$wpars$mcmc_n_CR

  # selection probability of crossover values
  .$mcmc$p_CR[] <- 1 / .$wpars$mcmc_n_CR

  # vector that stores how many times crossover value indices are used
  # ALJ: initialized to 1's in order to avoid numeric issues
  # debug: maybe noteworthy that this was originally initialized to 0's in Vrugt's algorithm
  # debug: play around with whether or not this makes a difference
  .$mcmc$n_id[] <- 1

  # debug: print crossover values and crossover probabilities at the end of each iteration
  print('CR = ')
  print(.$mcmc$CR)
  print('p_CR = ')
  print(.$mcmc$p_CR)
}


# generate proposal using DREAM algorithm
proposal_generate_mcmc_dream <- function(., j ) {

  # debug
  # print(paste0('iteration = ', j))

  # reset matrix of jump vectors to zero
  .$mcmc$jump[] <- 0

  # current state, 'mcmc_chains' number of samples of a d-variate distribution
  # future work: play around with whether the current_state and jump matrices are absolutely essential
  .$mcmc$current_state[] <- matrix(.$dataf$pars_array[ , , j-1], nrow = .$dataf$lp, ncol = .$mcmc$d)

  # debug: can make sure that this code is exactly analagous with Matlab version of "sort" function
  # permute [1,2,...,mcmc_chains-1] mcmc_chains number of times
  .$mcmc$draw[]                        <- apply(matrix(runif((.$dataf$lp - 1) * .$dataf$lp), .$dataf$lp - 1, .$dataf$lp), 2, function(v) sort(v, index.return = T)$ix)
  if(.$wpars$mcmc_debug) .$mcmc$draw[] <- apply(.$mcmc$draw_seed[ , , j], 2, function(v) sort(v, index.return = T)$ix)

  # debug: can make sure this code snipit is exactly analagous to Matlab version
  # create a .$dataf$lp x 1 matrix of continuous uniform random values between -c_rand and c_rand
  .$mcmc$lambda[]                        <- matrix(runif(.$dataf$lp * 1, -.$wpars$mcmc_c_rand, .$wpars$mcmc_c_rand), .$dataf$lp)
  if(.$wpars$mcmc_debug) .$mcmc$lambda[] <- .$mcmc$lambda_seed[ , , j]

  # compute standard deviation of each dimension (ie,compute standard deviation of each column of current_state matrix)
  .$mcmc$sd_state[]      <- apply(.$mcmc$current_state, 2, sd)
  # debug: replace any 0's in standard deviation array with 1e-9 to avoid division by 0
  idx                    <- which(.$mcmc$sd_state == 0)
  .$mcmc$sd_state[idx]   <- 1e-9

  # create proposals
  # future work: vectorize this for-loop to improve computational efficiency, but this is non-trivial
  for (ii in 1:.$dataf$lp) {

    # select delta (equal selection probability) (ie, choose 1 value from the vector [1:delta] with replacement)
    D <- sample(1:.$wpars$mcmc_delta, 1, replace = T)

    # debugging: can make sure this code snipit is exactly analagous to Matlab version (maybe there is a numerical issue here?)
    # extract vectors a and b not equal to ii
    a <- .$mcmc$R[ii, .$mcmc$draw[1:D, ii]]
    b <- .$mcmc$R[ii, .$mcmc$draw[(D + 1):(2 * D), ii]]

    # select index of crossover value (weighted sample with replacement)
    .$mcmc$id <- sample(1:.$wpars$mcmc_n_CR, 1, replace = T, prob = .$mcmc$p_CR)

    # draw d values from uniform distribution between 0 and 1
    zz                        <- runif(.$mcmc$d)
    if(.$wpars$mcmc_debug) zz <- .$mcmc$zz_seed[ii, 1:.$mcmc$d, j]

    # derive subset A of selected dimensions
    A  <- which(zz < .$mcmc$CR[.$mcmc$id])

    #  how many dimensions are sampled
    d_star <- length(A)

    # make sure that A contains at least one value
    if (d_star == 0) {
      A <- which.min(zz)
      d_star <- 1
    }

    # calculate jump rate
    gamma_d <- 2.38 / sqrt(2 * D * d_star)

    # select gamma
    gamma <- sample(c(gamma_d, 1), size = 1, replace = T, prob = c(1 - .$wpars$mcmc_p_gamma, .$wpars$mcmc_p_gamma))

    # compute jump differential evolution of ii-th chain
    .$mcmc$jump[ii, A] <- .$wpars$mcmc_c_ergod * rnorm(d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[a, A] - .$mcmc$current_state[b, A]), dim = 1)

    # compute proposal of ii-th chain
    .$dataf$pars[ii, 1:.$mcmc$d] <- .$mcmc$current_state[ii, 1:.$mcmc$d] + .$mcmc$jump[ii, 1:.$mcmc$d]

    # call boundary handling function
    for (jj in 1:.$mcmc$d) .$boundary_handling(ii = ii, jj = jj)

  }
}


# proposal acceptance function for the DREAM algorithm
proposal_accept_mcmc_dream <- function(., j, lklihood) {

  # future work: parallelize this part

  # debug: make sure this is being assigned at the right place in terms of function call order; also check function call order again
  # future work: play around with whether p_state assignment is absolutely essential
  # likelihood of current state
  .$mcmc$p_state[] <- .$dataf$pars_lklihood[ ,j-1]

  for (ii in 1:.$dataf$lp) {

    # debug: maybe try this part with the division and without the exponential to see it if makes a differnce?
    # compute Metropolis acceptance probability
    alpha <- min(1, exp(lklihood[ii] - .$mcmc$p_state[ii]))

    # future work: figure out how to make this clunky block of code "prettier"

    runif_val                         <- runif(1, min = 0, max = 1)
    if (.$wpars$mcmc_debug) runif_val <- .$mcmc$runif_seed[ii, j]

    # determine if p_acc is larger than random number drawn from uniform distribution on interval [0,1]
    if (alpha > runif_val) {

      # future work: may not need this accept vector since no longer using de-mc algorithm
      # accept the proposal
      accept <- TRUE

      .$mcmc$current_state[ii, 1:.$mcmc$d]   <- .$dataf$pars[ii, 1:.$mcmc$d]
      .$mcmc$p_state[ii]                     <- lklihood[ii]

      # append accepted current_state and probability density to storage data frames
      .$dataf$pars_array[ii, 1:.$mcmc$d, j]  <- .$mcmc$current_state[ii, 1:.$mcmc$d]
      .$dataf$pars_lklihood[ii, j]           <- .$mcmc$p_state[ii]

      # if debugging, store generated proposal (regardless of whether accepted or not)
      if (.$wpars$mcmc_debug) .$dataf$prop_storage[ii, 1:.$mcmc$d, j] <- .$dataf$pars[ii, 1:.$mcmc$d]

    } else {

      # reject the proposal
      accept <- FALSE

      # set jump back to zero for p_CR
      .$mcmc$jump[ii, 1:.$mcmc$d] <- 0

      # repeat previous current_state and probability density in storage data frames
      .$dataf$pars_array[ii, 1:.$mcmc$d, j]   <- .$dataf$pars_array[ii, 1:.$mcmc$d, j-1]
      .$dataf$pars_lklihood[ii, j]            <- .$dataf$pars_lklihood[ii, j-1]

      # if debugging, store generated proposal (regardless of whether accepted or not)
      if (.$wpars$mcmc_debug) .$dataf$prop_storage[ii, 1:.$mcmc$d, j] <- .$dataf$pars[ii, 1:.$mcmc$d]
    }

    # update jump distance crossover index
    .$mcmc$J[.$mcmc$id]    <- .$mcmc$J[.$mcmc$id] + sum((.$mcmc$jump[ii, 1:.$mcmc$d] / .$mcmc$sd_state)^2)

    # number of times index crossover is used
    .$mcmc$n_id[.$mcmc$id] <- .$mcmc$n_id[.$mcmc$id] + 1

    # future work: figure out how to align out_n with burn-in
    # out_n <- .$wpars$mcmc_maxiter / 2
    # APW: not a huge deal, but this is not quite right in the case where j==out_n+1 but not accepted as it adds the output from the rejected proposal
    out_n <- 0
    if (j > out_n)
      .$dataf$out_mcmc[ii,,(j-out_n)] <- if(accept | j == out_n + 1) .$dataf$out[ii, ] else .$dataf$out_mcmc[ii,,(j-out_n-1)]
  }

  # update selection probability of crossover in the first 10% of samples
  if ((j < (.$wpars$mcmc_maxiter / 10)) & (sum(.$mcmc$J) > 0)) {
    .$mcmc$p_CR <- .$mcmc$J    / .$mcmc$n_id
    .$mcmc$p_CR <- .$mcmc$p_CR / sum(.$mcmc$p_CR)
  }

  # debug: print crossover values and crossover probabilities at the end of each iteration
  print(paste0('iteration = ', j))
  print('CR = ')
  print(.$mcmc$CR)
  print('p_CR = ')
  print(.$mcmc$p_CR)

}

# function that generates/updates crossover values based on current probabilities
generate_CR <- function(.) {

}

# function that adapts crossover probabilities
adapt_CR <- function(.) {

}

# function that detects and corrects outlier Markov chains
outlier_handling <- function(.) {

}

# function that computes the R-statistic of Gelman and Rubin as a convergence diagnostic
Gelman_Rubin <- function(., j) {

  # TEMPORARY basing off of Dan's MATLAB code

  # need an R_stat storage array!

  # TEMPORARY dimensions
  # number of stored samples
  n   <- .$wpars$mcmc_maxiter
  # dimension / number of parameters
  nrY <- .$mcmc$d
  # number of Markov chains
  m   <- .$wpars$mcmc_chains

  # TEMPORARY reshape pars_array
  sequences <- aperm(.$dataf$pars_array, c(3, 2, 1))
  #print(dim(.$dataf$pars_array))
  #print(dim(sequences))

  if (n < 10) {

    # set R-statistic to large value
    R_stat <- -2 * rep(1, n)

  } else {

      # determin chain means
      meanSeq <- apply(sequences, 3, mean)
      print(length(meanSeq))


      R_stat <- 1

  }

  # within-chain variance (for each parameter)

  # between-chain variance (for each parameter)

  return(R_stat)
}



# Likelihood functions
################################

# expects model output to be probability - as in the output from the mixture model
f_proposal_lklihood_log <- function(.) {

  # derive log density
  log(.$dataf$out)

  # print(paste0('model likelihood = ', log(.$dataf$out)))

  return(log(.$dataf$out))
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
  # future work: screen out night-time values in a more efficient location
  sspos <- which(.$dataf$obsse > 1e-9)

  # number of measured data points (that do not have zero uncertainty)
  obs_n <- length(sspos)

  # observed error
  obsse <- if(.$wpars$mcmc_homosced)   rep(mean(.$dataf$obsse[sspos]), obs_n)
           else                .$dataf$obsse[sspos]

  # calculate error residual (each chain is on rows of dataf$out, take transpose)
  error_residual_matrix <- ( t(.$dataf$out)[sspos, ] - .$dataf$obs[sspos] ) / obsse

  # calculate sum of squared error
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2))

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  -(obs_n/2)*log(2*pi) - sum(log(obsse)) - 0.5*SSR
}



# Debuging functions
################################

# set seed function (to reproduce sequences of quasi-random numbers)
set_seed <- function(.) {

  # random number generation 1, set seed for draw matrix in proposal_generate_mcmc_dream
  set.seed(1703)
  .$mcmc$draw_seed[]   <- runif(((.$dataf$lp - 1) * .$dataf$lp * .$wpars$mcmc_maxiter), min = 0, max = 1)

  # random number generation 2, set seed for lambda matrix in proposal_generate_mcmc_dream
  set.seed(4050)
  .$mcmc$lambda_seed[] <- runif((.$dataf$lp * .$wpars$mcmc_maxiter), min = -.$wpars$mcmc_c_rand, max = .$wpars$mcmc_c_rand)

  # random number generation 3, set seed for zz vector in proposal_generate_mcmc_dream
  set.seed(1337)
  .$mcmc$zz_seed[]     <- runif((.$mcmc$d * .$dataf$lp * .$wpars$mcmc_maxiter), min = 0, max = 1)

  # random number generation 4, runif(1) value used in accept/reject step
  .$mcmc$runif_seed[]  <- runif((.$dataf$lp * .$wpars$mcmc_maxiter), min = 0, max = 1)

}

### END ###
