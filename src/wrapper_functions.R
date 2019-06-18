################################
#
# Wrapper functions
#
# AWalker July 2018
#
################################



# set seed functions
################################

# setting the seed to reproduce sequences of quasi-random numbers

# (1) set the seed for uniform_r generation
set_seed1 <- function(.) {
  set.seed(1703)
  .$mcmc$uniform_r_seed <- matrix(data = 0, nrow = .$wpars$mcmc_maxiter, ncol = .$dataf$lp)
  .$mcmc$uniform_r_seed[] <- runif((.$wpars$mcmc_maxiter * .$dataf$lp), min = -0.01, max = 0.01)
}

# (2) set the seed for R1 and R2 general_functions
set_seed2 <- function(.) {
  set.seed(4050)
}


# (3) set the seed for runif(1) value chosen in accept/reject step
set_seed3 <- function(.) {
  set.seed(1337)
  .$mcmc$runif_seed <- matrix(data = 0, nrow = .$wpars$mcmc_maxiter, ncol = .$dataf$lp)
  .$mcmc$runif_seed[] <- runif((.$wpars$mcmc_maxiter * .$dataf$lp), min = 0, max = 1)
}

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



# DEMC functions
################################

# generate proposal using DE-MC algorithm
proposal_generate_demc <- function(., j ) {

  # Metropolis sampling
  # scaling factor
  d <- ncol(.$dataf$pars)
  gamma_star <- 2.38 / sqrt(d + d)
  # b-value should be small compared to width of target distribution
  # APW: does this b parameter have a name?
  # ALJ: regrettably, there is no name for "b" that I have been able to find
  # ALJ: specifies range for "randomization" value added to proposal being generated
  b_rand  <- 0.01
  # ALJ: made uniform_r a vector and pulled it outside for-loop
  # draw vector of random numbers from uniform distribution on interval (-b_rand, b_rand)
  # uniform_r <- runif(d,min=(-b_rand),max=b_rand)
  # ALJ: NEED TO TRY creating uniform_r as just a randomly drawn scalar value
  # temporarily hardcode uniform_r
  # uniform_r <- rep(-0.000378, d)
  # ALJ: like below is how uniform_r is in the original r-script
  # uniform_r <- runif(1,min=(-b_rand),max=b_rand)

  # (1) set the seed for uniform_r generation
  # uniform_r <- .$mcmc$uniform_r_seed[j, ii]
  # uniform_r <- .$mcmc$uniform_r_seed[j, 1:.$dataf$lp]
  # uniform_r <- rep(.$mcmc$uniform_r_seed[j, ii], d)

  # print("uniform_r = ")
  # print(uniform_r)

  # ALJ: NEED TO TRY moving R1 and R2 iniitalization to 0 inside the for-loop
  # ALJ: index for 1st randomly chosen chain used in proposal generation
  # R1 <- 0
  # ALJ: index for 2nd randomly chosen chain used in proposal generation
  # R2 <- 0

  # evaluate for each chain
  for (ii in 1:.$dataf$lp) {

    # (1) set the seed for uniform_r generation
    uniform_r <- rep(.$mcmc$uniform_r_seed[j, ii], d)

    # randomly select two different numbers R1 and R2 unequal to j
    # from a uniform distribution without replacement
    R1 <- 0
    R2 <- 0

    while ((R1 == 0) | (R1 == ii))               R1 <- ceiling(runif(1,min=0,max=1)*.$dataf$lp)
    while ((R2 == 0) | (R2 == ii) | (R2 == R1))  R2 <- ceiling(runif(1,min=0,max=1)*.$dataf$lp)

    # temporarily hardcode R1 and R2
    # R1 <- 6
    # R2 <- 5

    print(paste0('<<<< iteration = ', j, ', chain = ', ii, ' <<<<'))

    #print('.$dataf$pars_array = ')
    #print(.$dataf$pars_array)

    #print(paste0('dim of .$dataf$pars_array = ', dim(.$dataf$pars_array)))

    #print(paste0('dim of .$dataf$pars = ', dim(.$dataf$pars)))

    #print(paste0('j = ', j))

    print(paste0("R1 = ", R1, ", R2 = ", R2))

    # evaluate for each parameter
    for (jj in 1:d) {

      #print(paste0('iteration = ', j, ', chain = ', ii))
      #print(paste0("R1 = ", R1, ", R2 = ", R2))

      # generate proposal via Differential Evolution

      # restructure this a bit for setting the seed
      # .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r[j]

      # this one below is the original
      # .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r[jj]

      # print out everything for debugging
      jump <- gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r[jj]
      .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + jump

      print(paste0('ii = ', ii))
      print(paste0('jj = ', jj))
      print(paste0('j = ', j))
      print(paste0('gamma_star = ', gamma_star))
      print(paste0('uniform_r[jj] = ', uniform_r[jj]))
      print(paste0('jump = ', jump))
      print(paste0('.$dataf$pars_array[R1, jj, j-1] = ', .$dataf$pars_array[R1, jj, j-1]))
      print(paste0('.$dataf$pars_array[R2, jj, j-1] = ', .$dataf$pars_array[R2, jj, j-1]))
      print(paste0('.$dataf$pars_array[ii, jj, j-1] = ', .$dataf$pars_array[ii, jj, j-1]))
      print(paste0('.$dataf$pars[ii, jj] = ', .$dataf$pars[ii, jj]))


      # .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r

      #print(paste0('uniform_r = ', uniform_r[jj]))

      # call boundary handling function
      # boundary_handling(., ii, jj )
    }

    print('proposal generated = ')
    print(.$dataf$pars[ii, ])

  }
}


# calculate proposal acceptance using the Metropolis ratio (for DE-MC algorithm)
proposal_accept_demc <- function(., j, lklihood ) {

  # Metropolis ratio
  metrop_ratio <- exp(lklihood - .$dataf$pars_lklihood[ ,j-1])

  # print(paste0('likelihood of proposal = ', lklihood))
  # print(paste0('likelihood of current chain = ', .$dataf$pars_lklihood[ ,j-1]))

  # print(paste0('dim of .$dataf$pars_lklihood ' = dim(.$dataf$pars_lklihood)))

  # print(paste0('Metropolis ratio = ', metrop_ratio))
  # print(length(metrop_ratio))

  alpha        <- pmin(1, metrop_ratio)

  # ALJ: maybe try computing alpha this way (more numerically comprehensive)
  # ALJ: if this change needs to be made, make the code prettier
  #alpha <- numeric(.$dataf$lp)
  #for (q in 1:.$dataf$lp) {
  #  if (.$dataf$pars_lklihood[q, j-1] > 0) {
  #    alpha[q] <- min(1, metrop_ratio[q])
  #  } else {
  #    alpha[q] <- 1
  #  }
  # }

  # print('alpha = ')
  # print(alpha)

  # (3) set the seed for runif(1) value chosen in accept/reject step
  # runif_val <- .$mcmc$runif_seed[j, ii]
  runif_val <- .$mcmc$runif_seed[j, 1:.$dataf$lp]
  # print('runif_val = ')
  # print(runif_val)

  # print(paste0('length of alpha = ', length(alpha)))
  # print(paste0('type of alpha = ', typeof(alpha)))
  # print(paste0('length of runif_val = ', length(runif_val)))
  # print(paste0('type of runif_val = ', typeof(runif_val)))

  # evaluate for each chain
  # APW: change this iteration counter to ii for consistency
  # ALJ: changed kk to ii
  for(ii in 1:.$dataf$lp) {

    print(paste0('<<<< iteration = ', j, ', chain = ', ii, ' <<<<'))

    # ALJ: maybe try putting alpha here
    #if (.$dataf$pars_lklihood[ ,j-1] > 0) {
    #  alpha[ii] <- min(1, metrop_ratio)
    #} else {
    #  alpha[ii] <- 1
    #}

    # print(paste0(' new alpha = ', alpha))

    # accept if Metropolis ratio > random number from uniform distribution on interval (0,1)
    # APW: should this be inside or outside of the loop, ie could draw a random outside the loop
    # ALJ: put this outside for-loop and index accept
    # accept <- log(alpha[ii]) > log(runif(1,min=0,max=1))

    # temporarily hard-code runif-value for debugging
    # accept <- log(alpha[ii]) > log(0.843332)

    # (3) set the seed for runif(1) value chosen in accept/reject step
    accept <- log(alpha[ii]) > log(runif_val[ii])

    #print(paste0('runif_val = ', runif_val[ii]))

    # print(paste0('length of accept = ', length(accept)))
    # print('accept = ')
    # print(accept)

    # print(paste0('accept = ', accept))

    .$dataf$pars_array[ii,,j]   <- if(accept) .$dataf$pars[ii,] else .$dataf$pars_array[ii,,j-1]
    .$dataf$pars_lklihood[ii,j] <- if(accept) lklihood[ii]      else .$dataf$pars_lklihood[ii,j-1]
    #print(c(ii, accept))
    #print(c(.$dataf$pars_array[ii,,j], .$dataf$pars[ii,], .$dataf$pars_array[ii,,j-1] ))

    # ALJ: temporarily getting rid of burn-in for debugging purposes
    # out_n <- .$wpars$mcmc_maxiter/2
    out_n <- 0
    if (j > out_n)
      .$dataf$out_mcmc[ii,,(j-out_n)] <- if(accept | j==out_n+1) .$dataf$out[ii,] else .$dataf$out_mcmc[ii,,(j-out_n-1)]
# APW: not a huge deal, but this is not quite right in the case where j==out_n+1 but not accepted as it adds the output from the rejected proposal
#    if ((j+1) > out_n)
#      .$dataf$out_mcmc[ii,,(j+1-out_n)] <- if(accept | j==out_n+1) .$dataf$out[ii,] else .$dataf$out_mcmc[ii,,(j-out_n-1)]
  }
}



# DREAM functions
################################

# static part of DREAM algorithm
static_dream <- function(.) {

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
  # APW: lets call this variable sd_state
  .$mcmc$std_state     <- numeric(.$mcmc$d)
  .$mcmc$jump          <- matrix(data=0, nrow=.$dataf$lp,   ncol=.$mcmc$d )
  .$mcmc$draw          <- matrix(data=0, nrow=.$dataf$lp-1, ncol=.$dataf$lp )
  .$mcmc$lambda        <- matrix(data=0, nrow=.$dataf$lp,   ncol=1 )

  # index of chains for Differential Evolution
  # APW: change this iteration counter to ii for consistency
  for (kk in 1:.$dataf$lp) .$mcmc$R[kk, ] <- setdiff(1:.$dataf$lp, kk)

  # crossover values
  .$mcmc$CR[] <- 1:.$wpars$mcmc_n_CR / .$wpars$mcmc_n_CR

  # selection probability of crossover values
  .$mcmc$p_CR[] <- rep(1, .$wpars$mcmc_n_CR) / .$wpars$mcmc_n_CR
  # make changes to []
  #.$mcmc$p_CR[] <- 1 / .$wpars$mcmc_n_CR # APW: this should work, [] preserves data structure
  #print('static DREAM')
  #print(.$mcmc$p_CR)

  # vector that stores how many times crossover value indices are used
  # initialized to 1's in order to avoid numeric issues
  .$mcmc$n_id[] <- rep(1, .$wpars$mcmc_n_CR)
  #.$mcmc$n_id[] <- 1 # APW: as above
  # originally initialized to 0's in Vrugt's algorithm
  #.$mcmc$n_id[] <- rep(0, .$wpars$mcmc_n_CR)
}


# generate proposal using DREAM algorithm
proposal_generate_dream <- function(., j ) {

  # reset matrix of jump vectors to zero
  .$mcmc$jump[]          <- matrix(data = 0)
  # .$mcmc$jump[] <- 0 # APW: I think you can just do this, the square brackets preserve the data structure

  # current state ('mcmc_chains' number of samples of a d-variate distribution)
  # APW: I think you can just use the pars_array here and later
  .$mcmc$current_state[] <- matrix(.$dataf$pars_array[ , , j-1], nrow = .$dataf$lp, ncol = .$mcmc$d)

  # permute [1,2,...,mcmc_chains-1] mcmc_chains number of times
  .$mcmc$draw[]          <- apply(matrix(runif((.$dataf$lp - 1) * .$dataf$lp), .$dataf$lp - 1, .$dataf$lp), 2, function(v) sort(v, index.return = T)$ix)

  # create a .$dataf$lp x 1 matrix of continuous uniform random values between -c_rand and c_rand
  .$mcmc$lambda[]        <- matrix(runif(.$dataf$lp * 1, -.$wpars$mcmc_c_rand, .$wpars$mcmc_c_rand), .$dataf$lp)

  # compute standard deviation of each dimension (ie,compute standard deviation of each column of current_state matrix)
  .$mcmc$std_state[]     <- apply(.$mcmc$current_state, 2, sd)

  # create proposals
  # ALJ: could vectorize this inner for-loop to improve computational efficiency, but this is non-trivial
  for (ii in 1:.$dataf$lp) {

    # select delta (equal selection probability) (ie, choose 1 value from the vector [1:delta] with replacement)
    D <- sample(1:.$wpars$mcmc_delta,1,replace=T)

    # extract vectors a and b not equal to ii
    a <- .$mcmc$R[ii, .$mcmc$draw[1:D, ii]]
    b <- .$mcmc$R[ii, .$mcmc$draw[(D + 1):(2 * D), ii]]

    # select index of crossover value (weighted sample with replacement)
    .$mcmc$id <- sample(1:.$wpars$mcmc_n_CR, 1, replace = T, prob = .$mcmc$p_CR)
    #print(.$wpars$mcmc_n_CR)
    #print(.$mcmc$p_CR)
    #print(.$mcmc$id)

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
    # note maybe there is a way to consolidate this chunk of code more efficiently
    temp1 <- c(gamma_d, 1)
    temp2 <- c(1 - .$wpars$mcmc_p_gamma, .$wpars$mcmc_p_gamma)
    gamma <- sample(temp1, 1, replace = T, prob = temp2)

    # compute jump differential evolution of ii-th chain
    .$mcmc$jump[ii,A] <- .$wpars$mcmc_c_ergod * rnorm(d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[a,A] - .$mcmc$current_state[b,A]), dim = 1)

    # compute proposal of ii-th chain
    .$dataf$pars[ii,1:.$mcmc$d] <- .$mcmc$current_state[ii,1:.$mcmc$d] + .$mcmc$jump[ii,1:.$mcmc$d]

    # call boundary handling function
    for (jj in 1:.$mcmc$d) boundary_handling(., ii, jj )
  }
}


# proposal acceptance function for the DREAM algorithm
# in the future: could probably consolidate this with the acceptance function for the DE-MC algorithm
proposal_accept_dream <- function(., j, lklihood) {

  # likelihood of current state
  # APW: is this assignment totally necessary? Could we just use the pars_lklihood matrix?
  .$mcmc$p_state[] <- .$dataf$pars_lklihood[ ,j-1]

  # ALJ: make an accept vector of true/false values and pull it outside of for-loop (like done above with DE-MC)

  # APW: change this iteration counter to ii for consistency
  for (qq in 1:.$dataf$lp) {

    # APW: as with demc these two steps could be extracted to get a boolean accept vector
    # compute Metropolis acceptance probability
    alpha <- min(1, exp(lklihood[qq] - .$mcmc$p_state[qq]))

    # in the future: figure out how to make this clunky block of code prettier
    # determine if p_acc is larger than random number drawn from uniform distribution on interval [0,1]
    if (alpha > runif(1, min = 0, max = 1)) {

      # accept the proposal
      accept <- TRUE

      # APW: are these assignments totally necessary? Could we just use the pars & lklihood matrix / vector?
      .$mcmc$current_state[qq, 1:.$mcmc$d]  <- .$dataf$pars[qq,1:.$mcmc$d]
      .$mcmc$p_state[qq]                    <- lklihood[qq]

      # append accepted current_state and probability density to storage data frames
      .$dataf$pars_array[qq,1:.$mcmc$d,j]   <- .$mcmc$current_state[qq,1:.$mcmc$d]
      .$dataf$pars_lklihood[qq,j]           <- .$mcmc$p_state[qq]

      # store generated proposal (regardless of whether accepted or not)
      .$dataf$prop_storage[qq,1:.$mcmc$d,j] <- .$dataf$pars[qq,1:.$mcmc$d]

    } else {

      # reject the proposal
      accept <- FALSE

      # set jump back to zero for p_CR
      .$mcmc$jump[qq, 1:.$mcmc$d] <- 0

      # repeat previous current_state and probability density in storage data frames
      .$dataf$pars_array[qq, 1:.$mcmc$d,j]   <- .$dataf$pars_array[qq,1:.$mcmc$d,j-1]
      .$dataf$pars_lklihood[qq,j]            <- .$dataf$pars_lklihood[qq,j-1]

      # store generated proposal (regardless of whether accepted or not)
      .$dataf$prop_storage[qq, 1:.$mcmc$d,j] <- .$dataf$pars[qq, 1:.$mcmc$d]
    }

    # update jump distance crossover index
    # APW: this step is introducing NaNs when a value in std_state is 0
    .$mcmc$J[.$mcmc$id]    <- .$mcmc$J[.$mcmc$id] + sum((.$mcmc$jump[qq, 1:.$mcmc$d] / .$mcmc$std_state)^2)

    # number of times index crossover is used
    .$mcmc$n_id[.$mcmc$id] <- .$mcmc$n_id[.$mcmc$id] + 1

    # NOTE this chunk is a reapeat of DE-MC code
    # APW: can either move to mcmc function in wrapper (would require an accept vector)
    #      or, create a function that is called here
    out_n <- .$wpars$mcmc_maxiter / 2
    if (j > out_n)
      .$dataf$out_mcmc[qq,,(j-out_n)] <- if(accept | j == out_n + 1) .$dataf$out[qq, ] else .$dataf$out_mcmc[qq,,(j-out_n-1)]
  }

  print(.$mcmc$id)
  print(.$mcmc$n_id)
  print(sum(.$mcmc$jump[qq,]))
  print(.$mcmc$std_state)
  print(.$mcmc$J)
  print(.$mcmc$p_CR)

  # update selection probability of crossover
  # altered original algorithm here to account for numerical issues
  # APW: still something odd going on here.
  #      with the linear example mcmc$J has two zeros in vector elements 2 & 3
  #      and then at iteration 54 and after the first element is getting a NaN
  # APW: Actually the behaviour is weirder than that, I've seen NaNs come in anywhere from iteration 12 and up
  #      this only happens when .$mcmc$id is 1,
  #      and for any full MCMC run using linear test, .$mcmc$id is always the same value
  # APW: this is probably my fault from moving the variables from .$mcmc to .$wpars but I can't see where
  #      seems to be because mcmc$p_CR always contains just 1 and 0's
  # APW: update, p_CR coverges to a 1 and 0's acfter the second proposal, perhaps this is expected behaviour but worth more investigation
  #print(c(j,.$wpars$mcmc_maxiter,.$mcmc$J,sum(.$mcmc$J)))

  # test less than versus greater than
  # testing below - to prevent immediate transition to 1/0 states for p_CR, seems to work
  if ((j > (.$wpars$mcmc_maxiter / 10)) & (sum(.$mcmc$J) > 0)) {
  #if ((j < (.$wpars$mcmc_maxiter / 10)) & (sum(.$mcmc$J) > 0)) {
    .$mcmc$p_CR <- .$mcmc$J    / .$mcmc$n_id
    .$mcmc$p_CR <- .$mcmc$p_CR / sum(.$mcmc$p_CR)
  }
}

# function for detection and correction of outlier chains
#outlier_check <- function(.) {
#}

# test for convergence using the R-statistic convergence diagnostic of Gelman and Rubin
#chain_converge <- function(.) {
#}



# Likelihood functions
################################

# expects model output to be probability - as in the output from the mixture model
f_proposal_lklihood_log <- function(.) {
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
