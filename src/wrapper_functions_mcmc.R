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

  # Metropolis sampling

  # scaling factor
  d <- ncol(.$dataf$pars)
  gamma_star <- 2.38 / sqrt(d + d)

  # b-value should be small compared to width of target distribution
  # b_rand specifies range for drawn "randomization" value
  b_rand  <- 0.01

  # debug: temporarily hardcode uniform_r
  # uniform_r <- rep(-0.000378, d)

  uniform_r <- runif(1, min = (-b_rand), max = b_rand)

  # debug: (1) set the seed for uniform_r generation
  # uniform_r <- .$mcmc$uniform_r_seed[j, ii]
  # uniform_r <- .$mcmc$uniform_r_seed[j, 1:.$dataf$lp]
  # uniform_r <- rep(.$mcmc$uniform_r_seed[j, ii], d)

  # print("uniform_r = ")
  # print(uniform_r)

  # evaluate for each chain
  for (ii in 1:.$dataf$lp) {

    # debug: (1) set the seed for uniform_r generation
    # uniform_r <- rep(.$mcmc$uniform_r_seed[j, ii], d)

    # debug: R1 and R2 iniitalization to 0 inside the for-loop
    # randomly select two different numbers R1 and R2 unequal to j
    # from a uniform distribution without replacement
    R1 <- 0
    R2 <- 0

    while ((R1 == 0) | (R1 == ii))               R1 <- ceiling(runif(1,min=0,max=1)*.$dataf$lp)
    while ((R2 == 0) | (R2 == ii) | (R2 == R1))  R2 <- ceiling(runif(1,min=0,max=1)*.$dataf$lp)

    # debug: temporarily hardcode R1 and R2
    # R1 <- 6
    # R2 <- 5

    # print(paste0('<<<< iteration = ', j, ', chain = ', ii, ' <<<<'))

    # print('.$dataf$pars_array = ')
    # print(.$dataf$pars_array)

    # debug: check dimensions
    # print(paste0('dim of .$dataf$pars_array = ', dim(.$dataf$pars_array)))
    # print(paste0('dim of .$dataf$pars = ', dim(.$dataf$pars)))

    # print(paste0("R1 = ", R1, ", R2 = ", R2))

    # debug: store R1 and R2 values that were generated
    # .$dataf$R1_R2_storage[ii, 1, j] <- R1
    # .$dataf$R1_R2_storage[ii, 2, j] <- R2

    # evaluate for each parameter
    for (jj in 1:d) {

      # generate proposal via Differential Evolution

      # debug: restructure to make setting the seed easier
      # .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r[j]

      # debug: note that this one below was the original at beginning of bug-fixing process
      # .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r[jj]

      # debug: create "scalar" to measure jump differential
      # jump <- gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r[jj]
      # .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + jump

      # store uniform_r value to check that set.seed() code is functioning
      # .$dataf$uniform_r_storage[ii, jj, j] <- uniform_r[jj]

      # print(paste0('ii = ', ii))
      # print(paste0('jj = ', jj))
      # print(paste0('j = ', j))

      # print(paste0('gamma_star = ', gamma_star))

      # print(paste0('uniform_r[jj] = ', uniform_r[jj]))

      # print(paste0('jump = ', jump))

      # print(paste0('.$dataf$pars_array[R1, jj, j-1] = ', .$dataf$pars_array[R1, jj, j-1]))
      # print(paste0('.$dataf$pars_array[R2, jj, j-1] = ', .$dataf$pars_array[R2, jj, j-1]))
      # print(paste0('.$dataf$pars_array[ii, jj, j-1] = ', .$dataf$pars_array[ii, jj, j-1]))
      # print(paste0('.$dataf$pars[ii, jj] = ', .$dataf$pars[ii, jj]))

      .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r

      # print(paste0('uniform_r = ', uniform_r[jj]))

      # debug: temporarily comment out boundary handling
      # call boundary handling function
      #boundary_handling(., ii, jj )
      .$boundary_handling(ii=ii,jj=jj)
    }

    # print('proposal generated = ')
    # print(.$dataf$pars[ii, ])

    # store proposals being generated (regardless of whether or not they are being accepted)
    # .$dataf$prop_storage[ii, 1:d, j] <- .$dataf$pars[ii, 1:d]

  }

  # debug: print statements at the end of each run
  # if (j == .$wpars$mcmc_maxiter) {
    # print('array of proposals that were generated = ')
    # print(.$dataf$prop_storage)
    # print('uniform_r values (randomly generated) = ')
    # print(.$dataf$uniform_r_storage)
    # print('R1 (1st col.) and R2 (2nd col.) values (randomly generated) = ')
    # print(.$dataf$R1_R2_storage)
  # }

}


# calculate proposal acceptance using the Metropolis ratio (for DE-MC algorithm)
proposal_accept_mcmc_demc <- function(., j, lklihood ) {

  # Metropolis ratio
  metrop_ratio <- exp(lklihood - .$dataf$pars_lklihood[ ,j-1])

  .$dataf$metrop_ratio_storage[1:.$wpars$mcmc_chains, 1, j] <- t(metrop_ratio)

  # print(paste0('likelihood of proposal = ', lklihood))

  # print(paste0('likelihood of current chain = ', .$dataf$pars_lklihood[ ,j-1]))

  # print(paste0('dim of .$dataf$pars_lklihood ' = dim(.$dataf$pars_lklihood)))

  # print(paste0('Metropolis ratio = ', metrop_ratio))
  # print(length(metrop_ratio))

  alpha        <- pmin(1, metrop_ratio)

  .$dataf$alpha_storage[1:.$wpars$mcmc_chains, 1, j] <- t(alpha)

  # debug: maybe try computing alpha this way (more numerically comprehensive/rigorous)
  # alpha <- numeric(.$dataf$lp)
  # for (q in 1:.$dataf$lp) {
  #   if (.$dataf$pars_lklihood[q, j-1] > 0) {
  #     alpha[q] <- min(1, metrop_ratio[q])
  #   } else {
  #     alpha[q] <- 1
  #   }
  # }

  # print('alpha = ')
  # print(alpha)

  # debug: (3) set the seed for runif(1) value chosen in accept/reject step
  # runif_val <- .$mcmc$runif_seed[j, ii]
  # runif_val <- .$mcmc$runif_seed[j, 1:.$dataf$lp]

  # print('runif_val = ')
  # print(runif_val)

  # debug: store runif(1) value
  # .$dataf$runif_val_storage[1:.$wpars$mcmc_chains, 1, j] <- t(runif_val)

  # print(paste0('length of alpha = ', length(alpha)))
  # print(paste0('type of alpha = ', typeof(alpha)))

  # print(paste0('length of runif_val = ', length(runif_val)))
  # print(paste0('type of runif_val = ', typeof(runif_val)))

  # evaluate for each chain
  for(ii in 1:.$dataf$lp) {

    # print(paste0('<<<< iteration = ', j, ', chain = ', ii, ' <<<<'))

    # accept if Metropolis ratio > random number from uniform distribution on interval (0,1)
    accept <- log(alpha[ii]) > log(runif(1, min = 0, max = 1))

    # debug: temporarily hard-code runif-value for debugging
    # accept <- log(alpha[ii]) > log(0.843332)

    # debug:(3) set the seed for runif(1) value chosen in accept/reject step
    accept <- log(alpha[ii]) > log(runif_val[ii])

    # debug: store accept value
    # .$dataf$accept_storage[ii, 1, j] <- accept

    # print(paste0('runif_val = ', runif_val[ii]))

    # print(paste0('length of accept = ', length(accept)))

    # print('accept = ')
    # print(accept)

    .$dataf$pars_array[ii,,j]   <- if(accept) .$dataf$pars[ii,] else .$dataf$pars_array[ii,,j-1]
    .$dataf$pars_lklihood[ii,j] <- if(accept) lklihood[ii]      else .$dataf$pars_lklihood[ii,j-1]


    # debug: store every model evaluation debugging purposes
    # out_n <- .$wpars$mcmc_maxiter / 2
    out_n <- 0
    if (j > out_n)
      .$dataf$out_mcmc[ii,,(j-out_n)] <- if(accept | j==out_n+1) .$dataf$out[ii,] else .$dataf$out_mcmc[ii,,(j-out_n-1)]
    # APW: not a huge deal, but this is not quite right in the case where j==out_n+1 but not accepted as it adds the output from the rejected proposal
    # if ((j+1) > out_n)
    #   .$dataf$out_mcmc[ii,,(j+1-out_n)] <- if(accept | j==out_n+1) .$dataf$out[ii,] else .$dataf$out_mcmc[ii,,(j-out_n-1)]
  }

  # debug: print statements
  # if (j == .$wpars$mcmc_maxiter) {
    # print('Metropolis ratio = ')
    # print(.$dataf$metrop_ratio_storage)
    # print('alpha = ')
    # print(.$dataf$alpha_storage)
    # print('randomly generated runif_val = ')
    # print(.$dataf$runif_val_storage)
    # print('accept/reject values (1 = TRUE and 0 = FALSE) = ')
    # print(.$dataf$accept_storage)
  # }

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
}


# generate proposal using DREAM algorithm
proposal_generate_mcmc_dream <- function(., j ) {

  # reset matrix of jump vectors to zero
  .$mcmc$jump[] <- 0

  # current state, 'mcmc_chains' number of samples of a d-variate distribution
  # APW: I think you can just use the pars_array here and later
  # ALJ: you could, but I am hesitant to change it
  #      becasue using the current_state and jump matrices circumvents "function call order" issue we ran into with the DE-MC algorithm
  #      i.e., I think these "placeholder" matrices  are relevant to parallelization and/or other forms of the DREAM alg
  #      I would like to leave them until we are absolutely certain that the DREAM algorithm is up and running perfectly in MAAT
  # future work: play around with whether the current_state and jump matrices are absolutely essential
  .$mcmc$current_state[] <- matrix(.$dataf$pars_array[ , , j-1], nrow = .$dataf$lp, ncol = .$mcmc$d)

  # permute [1,2,...,mcmc_chains-1] mcmc_chains number of times
  .$mcmc$draw[]          <- apply(matrix(runif((.$dataf$lp - 1) * .$dataf$lp), .$dataf$lp - 1, .$dataf$lp), 2, function(v) sort(v, index.return = T)$ix)

  # create a .$dataf$lp x 1 matrix of continuous uniform random values between -c_rand and c_rand
  .$mcmc$lambda[]        <- matrix(runif(.$dataf$lp * 1, -.$wpars$mcmc_c_rand, .$wpars$mcmc_c_rand), .$dataf$lp)

  # compute standard deviation of each dimension (ie,compute standard deviation of each column of current_state matrix)
  .$mcmc$sd_state[]     <- apply(.$mcmc$current_state, 2, sd)

  # create proposals
  # future work: vectorize this for-loop to improve computational efficiency, but this is non-trivial
  for (ii in 1:.$dataf$lp) {

    # select delta (equal selection probability) (ie, choose 1 value from the vector [1:delta] with replacement)
    D <- sample(1:.$wpars$mcmc_delta, 1, replace = T)
    
    # extract vectors a and b not equal to ii
    a <- .$mcmc$R[ii, .$mcmc$draw[1:D, ii]]
    b <- .$mcmc$R[ii, .$mcmc$draw[(D + 1):(2 * D), ii]]

    # select index of crossover value (weighted sample with replacement)
    .$mcmc$id <- sample(1:.$wpars$mcmc_n_CR, 1, replace = T, prob = .$mcmc$p_CR)

    # draw d values from uniform distribution between 0 and 1
    zz <- runif(.$mcmc$d)

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
    # future work: figure out a way to make this chunk of code prettier
    temp1 <- c(gamma_d, 1)
    temp2 <- c(1 - .$wpars$mcmc_p_gamma, .$wpars$mcmc_p_gamma)
    gamma <- sample(temp1, 1, replace = T, prob = temp2)

    # compute jump differential evolution of ii-th chain
    .$mcmc$jump[ii,A] <- .$wpars$mcmc_c_ergod * rnorm(d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[a,A] - .$mcmc$current_state[b,A]), dim = 1)

    # compute proposal of ii-th chain
    .$dataf$pars[ii,1:.$mcmc$d] <- .$mcmc$current_state[ii,1:.$mcmc$d] + .$mcmc$jump[ii,1:.$mcmc$d]

    # debug: temporarily comment out boundary handling
    # call boundary handling function
    #for (jj in 1:.$mcmc$d) boundary_handling(., ii, jj )
    for (jj in 1:.$mcmc$d) .$boundary_handling(ii=ii,jj=jj)
  }
}


# proposal acceptance function for the DREAM algorithm
proposal_accept_mcmc_dream <- function(., j, lklihood) {

  # likelihood of current state
  # APW: is this assignment totally necessary? Could we just use the pars_lklihood matrix?
  # ALJ: same reasoning as above
  # future work: play around with whether p_state is absolutely essential
  .$mcmc$p_state[] <- .$dataf$pars_lklihood[ ,j-1]

  for (ii in 1:.$dataf$lp) {

    # compute Metropolis acceptance probability
    alpha <- min(1, exp(lklihood[ii] - .$mcmc$p_state[ii]))

    # future work: figure out how to make this clunky block of code prettier

    # determine if p_acc is larger than random number drawn from uniform distribution on interval [0,1]
    if (alpha > runif(1, min = 0, max = 1)) {

      # accept the proposal
      accept <- TRUE

      # APW: are these assignments totally necessary? Could we just use the pars & lklihood matrix / vector?
      # ALJ: no, I don't think all of these assignments are totally necessary
      #      but I left them all in to remain as true as possible to Vrugt's original pseudocode
      #      and because they make debugging easier
      #      once DREAM is integrated in the new MAAT (and we're confident that it's working), I can alter some of the unnecessary assignments
      .$mcmc$current_state[ii, 1:.$mcmc$d]  <- .$dataf$pars[ii,1:.$mcmc$d]
      .$mcmc$p_state[ii]                    <- lklihood[ii]

      # append accepted current_state and probability density to storage data frames
      .$dataf$pars_array[ii,1:.$mcmc$d,j]   <- .$mcmc$current_state[ii,1:.$mcmc$d]
      .$dataf$pars_lklihood[ii,j]           <- .$mcmc$p_state[ii]

      # debug: store generated proposal (regardless of whether accepted or not)
      # .$dataf$prop_storage[ii,1:.$mcmc$d,j] <- .$dataf$pars[ii,1:.$mcmc$d]

    } else {

      # reject the proposal
      accept <- FALSE

      # set jump back to zero for p_CR
      .$mcmc$jump[ii, 1:.$mcmc$d] <- 0

      # repeat previous current_state and probability density in storage data frames
      .$dataf$pars_array[ii, 1:.$mcmc$d,j]   <- .$dataf$pars_array[ii,1:.$mcmc$d,j-1]
      .$dataf$pars_lklihood[ii,j]            <- .$dataf$pars_lklihood[ii,j-1]

      # debug: store generated proposal (regardless of whether accepted or not)
      # .$dataf$prop_storage[ii, 1:.$mcmc$d,j] <- .$dataf$pars[ii, 1:.$mcmc$d]
    }

    # update jump distance crossover index
    # APW: this step is introducing NaNs when a value in sd_state is 0
    # ALJ: I think this issue is now resolved?
    .$mcmc$J[.$mcmc$id]    <- .$mcmc$J[.$mcmc$id] + sum((.$mcmc$jump[ii, 1:.$mcmc$d] / .$mcmc$sd_state)^2)

    # number of times index crossover is used
    .$mcmc$n_id[.$mcmc$id] <- .$mcmc$n_id[.$mcmc$id] + 1

    # debug: see above comments in analagous DE-MC function
    # out_n <- .$wpars$mcmc_maxiter / 2
    out_n <- 0
    if (j > out_n)
      .$dataf$out_mcmc[ii,,(j-out_n)] <- if(accept | j == out_n + 1) .$dataf$out[ii, ] else .$dataf$out_mcmc[ii,,(j-out_n-1)]
  }

  # print(.$mcmc$id)

  # print(.$mcmc$n_id)

  # print(sum(.$mcmc$jump[ii,]))

  # print(.$mcmc$sd_state)

  # print(.$mcmc$J)

  # print(.$mcmc$p_CR)

  # update selection probability of crossover
  # ALJ: altered original algorithm here to account for numerical issues
  # APW: still something odd going on here.
  #      with the linear example mcmc$J has two zeros in vector elements 2 & 3
  #      and then at iteration 54 and after the first element is getting a NaN
  # APW: Actually the behaviour is weirder than that, I've seen NaNs come in anywhere from iteration 12 and up
  #      this only happens when .$mcmc$id is 1,
  #      and for any full MCMC run using linear test, .$mcmc$id is always the same value
  # APW: this is probably my fault from moving the variables from .$mcmc to .$wpars but I can't see where
  #      seems to be because mcmc$p_CR always contains just 1 and 0's
  # APW: update, p_CR coverges to a 1 and 0's acfter the second proposal, perhaps this is expected behaviour but worth more investigation
  # debug: test less than versus greater than
  # debug: testing below - to prevent immediate transition to 1/0 states for p_CR, seems to work
  #print(j)
  #print(.$wpars$mcmc_maxiter / 10)
  #print(.$mcmc$J)
  #print(sum(.$mcmc$J))

  if ((j > (.$wpars$mcmc_maxiter / 10)) & (sum(.$mcmc$J) > 0)) {
  # if ((j < (.$wpars$mcmc_maxiter / 10)) & (sum(.$mcmc$J) > 0)) {
    .$mcmc$p_CR <- .$mcmc$J    / .$mcmc$n_id
    .$mcmc$p_CR <- .$mcmc$p_CR / sum(.$mcmc$p_CR)
  }
}

# future work: function for detection and correction of outlier chains
# outlier_check <- function(.) {
# }

# future work: test for convergence using the R-statistic convergence diagnostic of Gelman and Rubin
# chain_converge <- function(.) {
# }



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
  #print(error_residual_matrix)

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

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  -(obs_n/2)*log(2*pi) - sum(log(obsse)) - 0.5*SSR
}



# Debuging functions
################################

#set seed functions (to reproduce sequences of quasi-random numbers)
# debug: function (1), set seed for uniform_r generation
# set_seed1 <- function(.) {
#  set.seed(1703)
#  .$mcmc$uniform_r_seed <- matrix(data = 0, nrow = .$wpars$mcmc_maxiter, ncol = .$dataf$lp)
#  .$mcmc$uniform_r_seed[] <- runif((.$wpars$mcmc_maxiter * .$dataf$lp), min = -0.01, max = 0.01)
# }


# debug: function (2), set seed for R1 and R2 general_functions
# set_seed2 <- function(.) {
#   set.seed(4050)
# }


# debug: function (3), set seed for runif(1) value chosen in accept/reject step
# set_seed3 <- function(.) {
#   set.seed(1337)
#   .$mcmc$runif_seed <- matrix(data = 0, nrow = .$wpars$mcmc_maxiter, ncol = .$dataf$lp)
#   .$mcmc$runif_seed[] <- runif((.$wpars$mcmc_maxiter * .$dataf$lp), min = 0, max = 1)
# }



### END ###
