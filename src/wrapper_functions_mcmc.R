################################
#
# Wrapper functions for MCMC runs
#
# AJohnson, AWalker July 2019
#
################################



# DEMC functions
################################

init_mcmc_demc <- function(.) NULL

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



# DREAM MCMC functions
################################

# initialisation of DREAM algorithm
init_mcmc_dream <- function(.) {

  # number of parameters being estimated
  .$mcmc$d <- ncol(.$dataf$pars)

  # preallocate memory space for algorithmic variables
  .$mcmc$R             <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$dataf$lp - 1)
  .$mcmc$current_state <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$mcmc$d)
  .$mcmc$p_state       <- numeric(.$dataf$lp)
  .$mcmc$jump          <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$mcmc$d)
  .$mcmc$draw          <- matrix(data = 0, nrow = .$dataf$lp - 1, ncol = .$dataf$lp)
  .$mcmc$lambda        <- matrix(data = 0, nrow = .$dataf$lp, ncol = 1)

  # ALJ: preallocate space for crossover variables (non-adaptive)
  # ALJ: may not need these any more AND/OR need to figure out how these relate to new CR vars (may be redundant)
  # ALJ: these are non-adaptive crossover vars from vrugt matlab paper
  .$mcmc$id            <- numeric(1)
  .$mcmc$J             <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$n_id          <- numeric(.$wpars$mcmc_n_CR)

  # ALJ: preallocate space for crossover variables (adaptive)
  # ALJ: some may be redundant with above non-adaptive crossover vars
  .$mcmc$t         <- numeric(1)
  .$mcmc$CR_burnin <- numeric(1)
  .$mcmc$L         <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$del       <- numeric(.$wpars$mcmc_n_CR)

  # ALJ: preallocate space for crossover variables (ones that are the same regardless of whther adaptive or not)
  # ALJ: standard deviation computed at different places in the alg depending whether adaptive CR or not
  .$mcmc$sd_state <- numeric(.$mcmc$d)
  .$mcmc$p_CR     <- numeric(.$wpars$mcmc_n_CR)

  # if user chooses adaptive crossover probabilties
  #if(.$wpars$mcmc_adapt_CR) {
    # ALJ: CR can either be a vector or scalar, and I need to pick one way to do it
    # ALJ: in vrugt matlab paper CR is a vector; in vrugt 2008 paper CR is a scalar
    #.$mcmc$CR     <- numeric(1)

    # initialize crossover variables
    #.$mcmc$t      <- 1
    #.$mcmc$L[]    <- 0

    # initial probability of each crossover value
    #.$mcmc$p_CR[] <- 1 / .$wpars$mcmc_n_CR

    # burn-in period for adapting crossover vals (compute number of iters)
    #.$mcmc$CR_burnin <- ceiling(.$wpars$mcmc_CR_burnin * .$wpars$mcmc_maxiter)

  #} else {
    # ALJ: CR can either be a vector or scalar, and I need to pick one way to do it
    # ALJ: in vrugt matlab paper CR is a vector; in vrugt 2008 paper CR is a scalar
    .$mcmc$CR   <- numeric(.$wpars$mcmc_n_CR)

    # crossover values
    .$mcmc$CR[] <- 1:.$wpars$mcmc_n_CR / .$wpars$mcmc_n_CR

    # selection probability of crossover values
    .$mcmc$p_CR[] <- 1 / .$wpars$mcmc_n_CR

    # vector that stores how many times crossover value indices are used
    # ALJ: initialized to 1's in order to avoid numeric issues
    # debug: maybe noteworthy that this was originally initialized to 0's in Vrugt's algorithm
    # debug: play around with whether or not this makes a difference and understand its purpose better
    .$mcmc$n_id[] <- 1
  #}

  # number of burn-in iterations
  .$mcmc$burnin <- ceiling(.$wpars$mcmc_burnin * .$wpars$mcmc_maxiter)

  # debug: print crossover values and crossover probabilities at the end of each iteration
  #print('check init_mcmc_dream part')
  #print(paste0('mcmc_adapt_CR is ', .$wpars$mcmc_adapt_CR))
  #print(paste0('CR_burnin = ', .$mcmc$CR_burnin))
  #print('CR = ')
  #print(.$mcmc$CR)
  #print('p_CR = ')
  #print(.$mcmc$p_CR)

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

}


# generate proposal using DREAM algorithm
proposal_generate_mcmc_dream <- function(., j ) {

  # debug
  print(paste0('iteration = ', j))

  # reset matrix of jump vectors to zero
  .$mcmc$jump[] <- 0

  # current state, 'mcmc_chains' number of samples of a d-variate distribution
  .$mcmc$current_state[] <- matrix(.$dataf$pars_array[ , , j-1], nrow = .$dataf$lp, ncol = .$mcmc$d)

  # debug: can make sure that this code is exactly analagous with Matlab version of "sort" function
  # permute [1,2,...,mcmc_chains-1] mcmc_chains number of times
  .$mcmc$draw[]                        <- apply(matrix(runif((.$dataf$lp - 1) * .$dataf$lp), .$dataf$lp - 1, .$dataf$lp), 2, function(v) sort(v, index.return = T)$ix)
  if(.$wpars$mcmc_debug) .$mcmc$draw[] <- apply(.$mcmc$draw_seed[ , , j], 2, function(v) sort(v, index.return = T)$ix)

  # debug: can make sure this code snipit is exactly analagous to Matlab version
  # create a .$dataf$lp x 1 matrix of continuous uniform random values between -c_rand and c_rand
  .$mcmc$lambda[]                        <- matrix(runif(.$dataf$lp * 1, -.$wpars$mcmc_c_rand, .$wpars$mcmc_c_rand), .$dataf$lp)
  if(.$wpars$mcmc_debug) .$mcmc$lambda[] <- .$mcmc$lambda_seed[ , , j]

  # ALJ: maybe can comment this out and write my own standard deviation calculation?
  # if not adapting crossover values, compute standard deviation of each parameter
  #if(!(.$wpars$mcmc_adapt_CR)) {
    # compute standard deviation of each dimension (ie,compute standard deviation of each column of current_state matrix)
    .$mcmc$sd_state[]    <- apply(.$mcmc$current_state, 2, sd)
    # bugfix: replace any 0's in standard deviation array with 1e-9 to avoid division by 0
    idx                  <- which(.$mcmc$sd_state == 0)
    .$mcmc$sd_state[idx] <- 1e-9
  #}

  # create proposals
  # future work: vectorize this for-loop to improve computational efficiency, but this is non-trivial
  for (ii in 1:.$dataf$lp) {

    # select delta (equal selection probability) (ie, choose 1 value from the vector [1:delta] with replacement)
    D <- sample(1:.$wpars$mcmc_delta, 1, replace = T)

    # debugging: can make sure this code snipit is exactly analagous to Matlab version (maybe there is a numerical issue here?)
    # extract vectors a and b not equal to ii
    a <- .$mcmc$R[ii, .$mcmc$draw[1:D, ii]]
    b <- .$mcmc$R[ii, .$mcmc$draw[(D + 1):(2 * D), ii]]

    # NEED TO FIGURE OUT A WAY TO SWITCH BACK AND FORTH BETWEEN ADAPTIVE AND NON-ADAPTIVE HERE
    # select index of crossover value (weighted sample with replacement)
    # ALJ: I think this should be being sampled from a multinomial distribution
    #      I don't think what I have right now is correct
    #      It should be M(CR, p_CR)
    .$mcmc$id <- sample(1:.$wpars$mcmc_n_CR, 1, replace = T, prob = .$mcmc$p_CR)
    # ALJ: I'm thinking maybe an if/else statement here and using mcmc$id as an adaptive crossover parm
    # INSERT GENERATE CROSSOVER VALUES FUNCTION CALL HERE
    # ALJ: maybe add more things to this subroutine once i figure things out better ?????
    #if (.$wpars$mcmc_adapt_CR) .$generate_CR()

    # draw d values from uniform distribution between 0 and 1
    zz                        <- runif(.$mcmc$d)
    if(.$wpars$mcmc_debug) zz <- .$mcmc$zz_seed[ii, 1:.$mcmc$d, j]

    # NEED TO FIGURE OUT A WAY TO SWITCH BACK AND FORTH BETWEEN ADAPTIVE AND NON-ADAPTIVE HERE
    # ALJ: I'm thinking maybe an if-else statement that computes A and d-star differenlty depending on adaptive criterion
    # derive subset A of selected dimensions
    A  <- which(zz < .$mcmc$CR[.$mcmc$id])
    # ALJ: try this instead? this seems like it would work better for adapt_CR <- F
    #      A  <- which(zz < (1 - .$mcmc$CR[.$mcmc$id]))
    #      need to double check -->
    #      A is subset of crossed-over dimensions or updated dimensions or dimensions not updating...

    # NEED TO FIGURE OUT A WAY TO SWITCH BACK AND FORTH BETWEEN ADAPTIVE AND NON-ADAPTIVE HERE
    # ALJ: is this instead of d' = d' - 1 ????? if not, how do i put d' in here ?????
    #  how many dimensions are sampled (i.e., how many dimenstions that will be updated jointly)
    d_star <- length(A)
    # numerical check: make sure that A contains at least one value
    if (d_star == 0) {
      A <- which.min(zz)
      d_star <- 1
    }

    # debug
    #print(paste0('.$mcmc$id = ', .$mcmc$id))
    #print('A = ')
    #print(A)
    #print(paste0('d_star = ', d_star))

    # calculate jump rate
    gamma_d <- 2.38 / sqrt(2 * D * d_star)

    # select gamma
    gamma <- sample(c(gamma_d, 1), size = 1, replace = T, prob = c(1 - .$wpars$mcmc_p_gamma, .$wpars$mcmc_p_gamma))

    # compute jump differential evolution of ii-th chain
    .$mcmc$jump[ii, A] <- .$wpars$mcmc_c_ergod * rnorm(d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[a, A] - .$mcmc$current_state[b, A]), dim = 1)

    # compute proposal of ii-th chain
    .$dataf$pars[ii, 1:.$mcmc$d] <- .$mcmc$current_state[ii, 1:.$mcmc$d] + .$mcmc$jump[ii, 1:.$mcmc$d]

    # call boundary handling function
    for (jj in 1:.$mcmc$d) .$mcmc_bdry_handling(j=j, ii = ii, jj = jj)

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

      # future work: make out_n dependent on whether or not mcmc_debug is true (ie, if it is, then store all model evaluations, otherwise...)

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

    #if (.$wpars$mcmc_adapt_CR & (.$mcmc$t < .$mcmc$CR_burnin) ) {

      # ALJ: probably better if I write my own sd code

      # standard deviation of columns of pars_array
      #.$mcmc$sd_state[]     <- apply(.$dataf$pars_array, 2, sd)
      #idx                   <- which(.$mcmc$sd_state == 0)
      #.$mcmc$sd_state[idx]  <- 1e-9

      # ALJ: make sure this math is correct; pretty sure i'm not using sum fxn correctly

      # compute squared normalized jumping distance
      #.$mcmc$del[.$mcmc$id] <- .$mcmc$del[.$mcmc$id] + sum(((.$dataf$pars_array[ii, 1:.$mcmc$d, .$mcmc$t] - .$dataf$pars_array[ii, 1:.$mcmc$d, .$mcmc$t-1]) / .mcmc$sd_state)^2)
      #temp <- numeric(.$mcmc$d)
      #temp <- rep(0, .$mcmc$d)
      #for (qq in 1:.$mcmc$d) {
      #  temp[qq] <- ((.$dataf$pars_array[ii, qq, .$mcmc$t] - .$dataf$pars_array[ii, qq, .$mcmc$t-1]) / .$mcmc$sd_state[qq])^2
      #}
      #.$mcmc$del[.$mcmc$id] <- .$mcmc$del[.$mcmc$id] + sum(temp)

      #print('.$mcmc$del = '); print(.$mcmc$del)

    #}

    # number of times index crossover is used
    .$mcmc$n_id[.$mcmc$id] <- .$mcmc$n_id[.$mcmc$id] + 1

    # future work: figure out how to align out_n with burn-in
    # out_n <- .$wpars$mcmc_maxiter / 2
    # APW: not a huge deal, but this is not quite right in the case where j==out_n+1 but not accepted as it adds the output from the rejected proposal
    out_n <- 0
    if (j > out_n)
      .$dataf$out_mcmc[ii,,(j-out_n)] <- if(accept | j == out_n + 1) .$dataf$out[ii, ] else .$dataf$out_mcmc[ii,,(j-out_n-1)]

  }

  # ADAPT CROSSOVER VALUES FUNCTION CALL HERE
  #if (.$wpars$mcmc_adapt_CR & (.$mcmc$t < .$mcmc$CR_burnin) ) .$adapt_CR()

  # ALJ: experiment with whether or not J is the same as delta and L is the same as n_id ....

  # update selection probability of crossover in the first 10% of samples
  if ((j < (.$wpars$mcmc_maxiter / 10)) & (sum(.$mcmc$J) > 0)) {
    .$mcmc$p_CR <- .$mcmc$J    / .$mcmc$n_id
    .$mcmc$p_CR <- .$mcmc$p_CR / sum(.$mcmc$p_CR)
  }

  # debug: print crossover values and crossover probabilities at the end of each iteration
  #print(paste0('iteration = ', j))
  #print('id = ')
  #print(.$mcmc$id)
  #print('del = ')
  #print(.$mcmc$del)
  #print('CR = ')
  #print(.$mcmc$CR)
  #print('p_CR = ')
  #print(.$mcmc$p_CR)
  #print('n_id =')
  #print(.$mcmc$n_id)
  #print('L = ')
  #print(.$mcmc$L)

}

# function that generates/updates crossover values based on current probabilities
generate_CR <- function(.) {

  .$mcmc$CR[.$mcmc$id]   <- .$mcmc$id / .$wpars$mcmc_n_CR

  # this might have duality with n_id?

  .$mcmc$L[.$mcmc$id]    <- .$mcmc$L[.$mcmc$id] + 1

}

# function that adapts crossover probabilities
adapt_CR <- function(.) {

  # update probability of different CR values
  for (m in 1:.$mcmc$n_CR) {
    p_CR[m]  <- (.$mcmc$t * .$wpars$mcmc_chains * (.$mcmc$del[m] / .$mcmc$L[m]) ) / sum(.$mcmc$del)
  }

  .$mcmc$t <- .$mcmc$t + 1

}


# boundary handling functions
################################

# set parameter boundaries from prior distributions
boundary_handling_set <- function(.) {
  # number of samples to be used in boundary handling
  n <- 1e4
  boundary_sample     <- lapply(.$dynamic$pars_eval, function(cs) eval(parse(text = cs)))
  .$mcmc$boundary_min <- unlist(lapply(boundary_sample, min))
  .$mcmc$boundary_max <- unlist(lapply(boundary_sample, max))
  rm(boundary_sample)

  # debugging/development
  print('max')
  print(.$mcmc$boundary_max)
  print('min')
  print(.$mcmc$boundary_min)
}

mcmc_bdry_handling_none <- function(., j, ii, jj) {
  if ((j == .$wpars$mcmc_maxiter) & (ii == .$wpars$mcmc_chains) & (jj == .$mcmc$d)) {
    print('No option was chosen for MCMC boundary handling.')
  }
}


# restrict parameter proposals that are beyond the boundary (set them to bound)
mcmc_bdry_handling_bound <- function(., j, ii, jj) {

  # if outside bound of parameter space, restrict proposal value to corresponding dimension minimum
  if      (.$dataf$pars[ii,jj] < .$mcmc$boundary_min[jj]) .$dataf$pars[ii,jj] <- .$mcmc$boundary_min[jj]

  # if outside bound of parameter space, restrict proposal to corresponding dimension maximum
  else if (.$dataf$pars[ii,jj] > .$mcmc$boundary_max[jj]) .$dataf$pars[ii,jj] <- .$mcmc$boundary_max[jj]
}



# convergence diagnostic functions
#####################################

# no testing for convergence
mcmc_converge_none <- function(., j) {
  if (j == .$wpars$mcmc_maxiter) print('No option was chosen to test for MCMC convergence.')
}


# compute the R-statistic of Gelman and Rubin as a convergence diagnostic
mcmc_converge_Gelman_Rubin <- function(., j) {

  print('Geman-Rubin convergence test being called')

  # TEMPORARY basing off of Dan's MATLAB code -- will make more MAAT-like later

  # need an R_stat storage array!

  # TEMPORARY dimensions (rename everything to make more MAAT-like later)
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
      #print(length(meanSeq))


      R_stat <- 1

  }

  # within-chain variance (for each parameter)

  # between-chain variance (for each parameter)

  return(R_stat)
}



# likelihood functions
################################

# expects model output to be probability - as in the output from the mixture model
f_proposal_lklihood_log <- function(.) {

  # derive log density
  log(.$dataf$out)

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
  -(obs_n / 2) * log(SSR)
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
  -(obs_n / 2) * log(2 * pi) - sum(log(obsse)) - 0.5 * SSR
}



# outlier handling functions
#####################################

# no outlier handling for Markov chians
mcmc_outlier_none <- function(., j) {
  if (j == .$wpars$mcmc_maxiter) print('No option was chosen to identify and correct outlier Markov chains during burn-in.')
}


# function that detects and corrects outlier Markov chains using the Inter Quartile-Range (IQR) statistic
mcmc_outlier_iqr <- function(., j) {

  # debug/development print statment
  print('IQR outlier test being called')

  counter <- j / .$wpars$mcmc_check_iter

  # extract last 50% of samples of each chain
  sbst <- .$dataf$pars_lklihood[1:.$wpars$mcmc_chains, (ceiling(j/2)):j]

  # take the mean of the log of the posterior densities and store in omega
  # IMPORANT: all current likelihood function options already return log-likelihood
  #           so it's not necessary to take the log of sbst components here
  #           but this may change in the future with different likelihood functions
  for (ii in 1:.$wpars$mcmc_chains) .$dataf$omega[ii, counter] <- mean(sbst[ii, ])

  # determine upper and lower quantiles of the N different chains
  q1 <- quantile(.$dataf$omega[1:.$wpars$mcmc_chains, counter], prob = 0.25, type = 1)
  q3 <- quantile(.$dataf$omega[1:.$wpars$mcmc_chains, counter], prob = 0.75, type = 1)

  # compute IQR statistic
  iqr <- q3 - q1
  # alternate way to compute IQR statistic
  # iqr <- IQR(.$dataf$omega[1:.$wpars$mcmc_chains, counter], type = 1)

  # determine which chains are outliers
  outliers <- which(.$dataf$omega[ , counter] < (q1 - 2 * iqr))

  # if outlier chains are detected
  if (length(outliers) > 0) {

    print(paste0('Outlier chain detected. Chain ', outliers, ' at iteration ', j))

    # replace outlier(s) by randomly choosing from the remaining chains
    replace_idx <- rep(0, length(outliers))
    for (qq in 1:length(outliers)) {
      # while ((replace_idx[qq] == 0) | (replace_idx[qq] == outliers[qq])) {
      # check: make sure no elements of replace_idx are equal to any other elements of outliers
      while ((replace_idx[qq] == 0) | (replace_idx[qq] %in% outliers)) {
        replace_idx[qq] <- ceiling(runif(1, min = 0, max = 1) * .$wpars$mcmc_chains)
      }
    }

    # check: make sure 2 or more outliers aren't being replaced by the same randomly-chosen chain
    while (any(duplicated(replace_idx))) {
      repeat_idx <- which(duplicated(replace_idx))
      for (qq in 1:length(repeat_idx)) {
          replace_idx[repeat_idx[qq]] <-  ceiling(runif(1, min = 0, max = 1) * .$wpars$mcmc_chains)
      }
    }

    # replace outlier chain(s)
    .$dataf$pars_array[outliers, 1, j] <- .$dataf$pars_array[replace_idx, 1, j]
    .$dataf$pars_lklihood[outliers, j] <- .$dataf$pars_lklihood[replace_idx, j]
  }

  # identifying and correcting outliers should only be done during burn-in
  # because it violates the balance of sampled chains and destroys reversibility
  # if outlier chain is detected, apply another burn-in period before generating posterior moments
  # ALJ: not sure how to implement this in code??? basically need to restart j???
}



# debuging functions
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
