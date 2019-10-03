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
  #d          <- ncol(.$dataf$pars)
  d          <- dim(.$dataf$pars)[1]
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
      #.$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r
      .$dataf$pars[jj,ii] <- .$dataf$pars_array[jj,ii,j-1] + gamma_star * (.$dataf$pars_array[jj,R1,j-1] - .$dataf$pars_array[jj,R2,j-1]) + uniform_r

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

    #.$dataf$pars_array[ii,,j]   <- if(accept) .$dataf$pars[ii,] else .$dataf$pars_array[ii,,j-1]
    .$dataf$pars_array[,ii,j]   <- if(accept) .$dataf$pars[,ii] else .$dataf$pars_array[,ii,j-1]
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
  #.$mcmc$d <- ncol(.$dataf$pars)
  .$mcmc$d <- dim(.$dataf$pars)[1]

  # preallocate memory space for algorithmic variables
  .$mcmc$p_state       <- numeric(.$dataf$lp)
  .$mcmc$R             <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$dataf$lp - 1)
  #.$mcmc$current_state <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$mcmc$d)
  .$mcmc$current_state <- matrix(data = 0, nrow = .$mcmc$d, ncol = .$dataf$lp)
  .$mcmc$draw          <- matrix(data = 0, nrow = .$dataf$lp - 1, ncol = .$dataf$lp)
  .$mcmc$lambda        <- matrix(data = 0, nrow = .$dataf$lp, ncol = 1)
  .$mcmc$jump          <- matrix(data = 0, nrow = .$mcmc$d, ncol = .$dataf$lp)

  # preallocate space for crossover variables
  # ALJ: these are non-adaptive/non-working crossover vars from vrugt matlab paper
  .$mcmc$id   <- numeric(1)
  .$mcmc$J    <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$n_id <- numeric(.$wpars$mcmc_n_CR)

  # preallocate space for crossover variables
  # ALJ: some may be redundant with above crossover vars?
  # ALJ: these are from vrugt 2009 paper and also dan lu's matlab code
  .$mcmc$t      <- numeric(1)
  .$mcmc$m      <- numeric(1)
  .$mcmc$d_star <- numeric(1)
  .$mcmc$L      <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$del    <- numeric(.$wpars$mcmc_n_CR)

  # preallocate space for crossover variables
  # ones that are the same regardless of whther adaptive or not
  .$mcmc$CR_burnin <- numeric(1)
  .$mcmc$sd_state  <- numeric(.$mcmc$d)
  .$mcmc$p_CR      <- numeric(.$wpars$mcmc_n_CR)

  # if user chooses adaptive crossover probabilties
  if(.$wpars$mcmc_adapt_CR) {

    # ALJ: in vrugt matlab paper CR is a vector; in vrugt 2008 paper CR is a scalar
    .$mcmc$CR <- numeric(1)

    # initialize crossover variables
    .$mcmc$t      <- 1
    .$mcmc$L[]    <- 0
    .$mcmc$d_star <- .$mcmc$d

    # initial probability of each crossover value
    .$mcmc$p_CR[] <- 1 / .$wpars$mcmc_n_CR

  } else {

    # ALJ: in this case, crossover probalities are computed and stored in a vector
    .$mcmc$CR   <- numeric(.$wpars$mcmc_n_CR)

    # crossover probabilities/values
    .$mcmc$CR[] <- 1:.$wpars$mcmc_n_CR / .$wpars$mcmc_n_CR

    # selection probability of crossover values
    .$mcmc$p_CR[] <- 1 / .$wpars$mcmc_n_CR

    # vector that stores how many times crossover value indices are used
    # ALJ: initialized to 1's in order to avoid numeric issues
    # ALJ: maybe noteworthy that this was originally initialized to 0's in Vrugt's algorithm
    # ALJ: try initializing to 0's and then index the divide-by-zero-problem-area by id
    .$mcmc$n_id[] <- 1
  }

  # burn-in period for adapting crossover vals (compute number of iters)
  # ALJ: need to make sure that CR_burnin aligns correclty with burnin (especially with the restarting part)
  .$mcmc$CR_burnin <- ceiling(.$wpars$mcmc_CR_burnin * .$wpars$mcmc_maxiter)

  # number of burn-in iterations
  # ALJ: not sure whether or not this is the best way to do burn-in?
  .$mcmc$burnin <- ceiling(.$wpars$mcmc_burnin * .$wpars$mcmc_maxiter)

  # index of chains for Differential Evolution
  for (ii in 1:.$dataf$lp) .$mcmc$R[ii, ] <- setdiff(1:.$dataf$lp, ii)

  # debug/dev
  print('check init_mcmc_dream part')
  print(paste0('mcmc_adapt_CR is ', .$wpars$mcmc_adapt_CR))
  print('CR = ')
  print(.$mcmc$CR)
  print('p_CR = ')
  print(.$mcmc$p_CR)

}


# generate proposal using DREAM algorithm
proposal_generate_mcmc_dream <- function(., j ) {

  # debug
  print(paste0('iteration = ', j))

  # reset matrix of jump vectors to zero
  .$mcmc$jump[] <- 0

  # current state, 'mcmc_chains' number of samples of a d-variate distribution
  .$mcmc$current_state[] <- matrix(.$dataf$pars_array[ , , j-1])

  # debug: can make sure that this code is exactly analagous with Matlab version of "sort" function
  # permute [1,2,...,mcmc_chains-1] mcmc_chains number of times
  .$mcmc$draw[] <- apply(matrix(runif((.$dataf$lp - 1) * .$dataf$lp), .$dataf$lp - 1, .$dataf$lp), 2, function(v) sort(v, index.return = T)$ix)

  # create a .$dataf$lp x 1 matrix of continuous uniform random values between -c_rand and c_rand
  .$mcmc$lambda[] <- matrix(runif(.$dataf$lp * 1, -.$wpars$mcmc_c_rand, .$wpars$mcmc_c_rand), .$dataf$lp)

  # ALJ: maybe can comment this out and write my own standard deviation calculation?
  # if not adapting crossover values, compute standard deviation of each dimension/parameter
  if(!(.$wpars$mcmc_adapt_CR)) {

    #.$mcmc$sd_state[]      <- apply(.$mcmc$current_state, 2, sd)
    .$mcmc$sd_state[]      <- apply(.$mcmc$current_state, 1, sd)

    # replace any 0's in standard deviation array with 1e-9 to avoid division by 0
    idx                    <- which(.$mcmc$sd_state == 0)
    .$mcmc$sd_state[idx]   <- 1e-9
  }

  # create proposals
  for (ii in 1:.$dataf$lp) {

    # select delta (equal selection probability) (ie, choose 1 value from the vector [1:delta] with replacement)
    D <- sample(1:.$wpars$mcmc_delta, 1, replace = T)

    # debugging: can make sure this code snipit is exactly analagous to Matlab version (maybe there is a numerical issue here?)
    # extract vectors a and b not equal to ii
    a <- .$mcmc$R[ii, .$mcmc$draw[1:D, ii]]
    b <- .$mcmc$R[ii, .$mcmc$draw[(D + 1):(2 * D), ii]]

    if (.$wpars$mcmc_adapt_CR) {

      # generate crossover probability
      .$generate_CR()

      # calculate jump rate (scaling factor)
      # ALJ: is it correct to use the d_star from the previous iteration?
      gamma_d <- 2.38 / sqrt(2 * D * .$mcmc$d_star)

      # select gamma
      # gamma = 1 approximately every 5 iterations (with default p_gamma = 0.2)
      # when gamma = 1, jump between different modes of the posterior
      gamma <- sample(c(gamma_d, 1), size = 1, replace = T, prob = c(1 - .$wpars$mcmc_p_gamma, .$wpars$mcmc_p_gamma))

      # compute jump differential evolution of ii-th chain
      # ALJ: NOT SURE THAT MY REARRANGING OF THE SUM(...) PART IS CORRECT, but i think it is
      # ALJ: would maybe be more efficient to structure the mcmc_adapt_CR <- T part with A as well in the future (might also need to do this for debugging)
      #.$mcmc$jump[1:.$mcmc$d, ii] <- .$wpars$mcmc_c_ergod * rnorm(d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[A, a] - .$mcmc$current_state[A, b]), dim = 1)
      .$mcmc$jump[1:.$mcmc$d, ii] <- .$wpars$mcmc_c_ergod * rnorm(.$mcmc$d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[1:.$mcmc$d, a] - .$mcmc$current_state[1:.$mcmc$d, b]), dim = 1)

      # debug/dev
      #print('mcmc$jump[1:.$mcmc$d, ii] = ')
      #print(.$mcmc$jump[1:.$mcmc$d, ii])

      # compute proposal of ii-th chain
      .$dataf$pars[1:.$mcmc$d, ii] <- .$mcmc$current_state[1:.$mcmc$d, ii] + .$mcmc$jump[1:.$mcmc$d, ii]

      # debug/dev
      #print('dataf$pars[1:.$mcmc$d, ii] =')
      #print(.$dataf$pars[1:.$mcmc$d, ii])

      # replace each element (jj = 1,...,d) of the proposal with the corresponding current_state element
      # using a binomial scheme with probability 1 - CR (where CR = crossover probability)
      # note: when CR = 1, all dimensions are updated jointly and d_star = d (ALJ: USE THIS AS A TEST)
      crossover <- logical(length = .$mcmc$d)
      for (jj in 1:.$mcmc$d) {
        if (runif(1, min = 0, max = 1) <= (1 - .$mcmc$CR)) {
          .$dataf$pars[jj, ii] <- .$mcmc$current_state[jj, ii]
          # ALJ: this is how it is in Vrugt 2009 paper, but this would go to 1 for an awfully small parameter dimension space
          #.$mcmc$d_star <- .$mcmc$d_star - 1
        } else {
          crossover <- T
        }
      }

      #print('new dataf$pars[1:.$mcmc$d, ii] =')
      #print(.$dataf$pars[1:.$mcmc$d, ii])

      # ALJ: not really sure about the d_star above (seem's like it would go to 1 very quickly?)
      # alternatively: d_star <- number of dimensions being updated (this is how it is in Vrugt MATLAB paper)
      .$mcmc$d_star <- length(crossover)

      #print(paste0('d_star = ', .$mcmc$d_star))

      # numerical check
      if (.$mcmc$d_star == 0) .$mcmc$d_star <- 1

      # the above code SHOULD effectively modify each dimension with probability CR each time a proposal vector is generated

    } else {

      # select index of crossover value (weighted sample with replacement)
      # this returns an integer
      # ALJ: I think this should be being sampled from a multinomial distribution, and I'm not convinced that's what this is
      .$mcmc$id <- sample(1:.$wpars$mcmc_n_CR, 1, replace = T, prob = .$mcmc$p_CR)

      # draw d values from uniform distribution between 0 and 1
      zz <- runif(.$mcmc$d)

      # derive subset A of selected dimensions
      A  <- which(zz < .$mcmc$CR[.$mcmc$id])

      #  how many dimensions are sampled (i.e., how many dimenstions that will be updated jointly)
      d_star <- length(A)

      # numerical check: make sure that A contains at least one value
      if (d_star == 0) {
        A <- which.min(zz)
        d_star <- 1
      }

      # calculate jump rate (scaling factor)
      gamma_d <- 2.38 / sqrt(2 * D * d_star)

      # select gamma
      # gamma = 1 approximately every 5 iterations (with default p_gamma = 0.2)
      # when gamma = 1, jump between different modes of the posterior
      gamma <- sample(c(gamma_d, 1), size = 1, replace = T, prob = c(1 - .$wpars$mcmc_p_gamma, .$wpars$mcmc_p_gamma))

      # compute jump differential evolution of ii-th chain
      .$mcmc$jump[A, ii] <- .$wpars$mcmc_c_ergod * rnorm(d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[A, a] - .$mcmc$current_state[A, b]), dim = 1)

      print('mcmc$jump[1:.$mcmc$d, ii] = ')
      print(.$mcmc$jump[1:.$mcmc$d, ii])

      # compute proposal of ii-th chain
      #.$dataf$pars[ii, 1:.$mcmc$d] <- .$mcmc$current_state[ii, 1:.$mcmc$d] + .$mcmc$jump[ii, 1:.$mcmc$d]
      .$dataf$pars[1:.$mcmc$d, ii] <- .$mcmc$current_state[1:.$mcmc$d, ii] + .$mcmc$jump[1:.$mcmc$d, ii]
      #.$dataf$pars[,ii] <- .$mcmc$current_state[,ii] + .$mcmc$jump[,ii]
    }

    # call boundary handling function
    for (jj in 1:.$mcmc$d) .$mcmc_bdry_handling(j=j, ii = ii, jj = jj)

  }
}


# proposal acceptance function for the DREAM algorithm
proposal_accept_mcmc_dream <- function(., j, lklihood) {

  # future work: parallelize this part maybe?

  # debug: make sure this is being assigned at the right place in terms of function call order; also check function call order again
  # likelihood of current state
  .$mcmc$p_state[] <- .$dataf$pars_lklihood[ ,j-1]

  for (ii in 1:.$dataf$lp) {

    # debug: maybe try this part with the division and without the exponential to see it if makes a differnce?
    # compute Metropolis acceptance probability
    alpha <- min(1, exp(lklihood[ii] - .$mcmc$p_state[ii]))

    # future work: figure out how to make this clunky block of code "prettier"

    # determine if p_acc is larger than random number drawn from uniform distribution on interval [0,1]
    if (alpha > runif(1, min = 0, max = 1)) {

      # future work: may not need this accept vector since no longer using de-mc algorithm
      # accept the proposal
      accept <- TRUE

      #.$mcmc$current_state[ii, 1:.$mcmc$d] <- .$dataf$pars[ii, 1:.$mcmc$d]
      .$mcmc$current_state[1:.$mcmc$d,ii] <- .$dataf$pars[1:.$mcmc$d,ii]
      #.$mcmc$current_state[,ii] <- .$dataf$pars[,ii]
      .$mcmc$p_state[ii]                   <- lklihood[ii]

      # append accepted current_state and probability density to storage data frames
      #.$dataf$pars_array[ii, 1:.$mcmc$d, j]  <- .$mcmc$current_state[ii, 1:.$mcmc$d]
      .$dataf$pars_array[1:.$mcmc$d,ii,j] <- .$mcmc$current_state[1:.$mcmc$d,ii]
      #.$dataf$pars_array[,ii,j]  <- .$mcmc$current_state[,ii]
      .$dataf$pars_lklihood[ii,j]        <- .$mcmc$p_state[ii]

      # if debugging, store generated proposal (regardless of whether accepted or not)
      #if (.$wpars$mcmc_debug) .$dataf$prop_storage[ii, 1:.$mcmc$d, j] <- .$dataf$pars[ii, 1:.$mcmc$d]
      if (.$wpars$mcmc_debug) .$dataf$prop_storage[1:.$mcmc$d,ii,j] <- .$dataf$pars[1:.$mcmc$d,ii]
      #if (.$wpars$mcmc_debug) .$dataf$prop_storage[,ii,j] <- .$dataf$pars[,ii]

      # future work: make out_n dependent on whether or not mcmc_debug is true (ie, if it is, then store all model evaluations, otherwise...)

    } else {

      # reject the proposal
      accept <- FALSE

      # set jump back to zero for p_CR
      .$mcmc$jump[1:.$mcmc$d,ii] <- 0
      #.$mcmc$jump[,ii] <- 0

      # repeat previous current_state and probability density in storage data frames
      #.$dataf$pars_array[ii, 1:.$mcmc$d, j]   <- .$dataf$pars_array[ii, 1:.$mcmc$d, j-1]
      .$dataf$pars_array[1:.$mcmc$d,ii,j]   <- .$dataf$pars_array[1:.$mcmc$d,ii,j-1]
      #.$dataf$pars_array[,ii,j]   <- .$dataf$pars_array[,ii,j-1]
      .$dataf$pars_lklihood[ii, j]          <- .$dataf$pars_lklihood[ii, j-1]

      # if debugging, store generated proposal (regardless of whether accepted or not)
      #if (.$wpars$mcmc_debug) .$dataf$prop_storage[ii, 1:.$mcmc$d, j] <- .$dataf$pars[ii, 1:.$mcmc$d]
      if (.$wpars$mcmc_debug) .$dataf$prop_storage[1:.$mcmc$d,ii,j] <- .$dataf$pars[1:.$mcmc$d,ii]
      #if (.$wpars$mcmc_debug) .$dataf$prop_storage[,ii,j] <- .$dataf$pars[,ii]
    }

    if (.$wpars$mcmc_adapt_CR & (.$mcmc$t < .$mcmc$CR_burnin) ) {

      # compute the squared normalized jumping distance
      .$calc_del(j = j, ii = ii)

    } else if (!.$wpars$mcmc_adapt_CR) {

      # PRETTY SURE THIS IS WHERE THE PROBLEM IS!
      # update jump distance crossover index
      # ALJ / debug: why is J computed after jump is reset back to 0 if proposal is rejected
      #.$mcmc$J[.$mcmc$id]    <- .$mcmc$J[.$mcmc$id] + sum((.$mcmc$jump[ii, 1:.$mcmc$d] / .$mcmc$sd_state)^2)
      .$mcmc$J[.$mcmc$id]    <- .$mcmc$J[.$mcmc$id] + sum((.$mcmc$jump[1:.$mcmc$d,ii] / .$mcmc$sd_state)^2)
      #.$mcmc$J[.$mcmc$id]    <- .$mcmc$J[.$mcmc$id] + sum((.$mcmc$jump[,ii] / .$mcmc$sd_state)^2)

      # number of times index crossover is used
      .$mcmc$n_id[.$mcmc$id] <- .$mcmc$n_id[.$mcmc$id] + 1

    }

    # future work: figure out how to align out_n with burn-in
    # out_n <- .$wpars$mcmc_maxiter / 2
    # APW: not a huge deal, but this is not quite right in the case where j==out_n+1 but not accepted as it adds the output from the rejected proposal
    out_n <- 0
    if (j > out_n)
      .$dataf$out_mcmc[ii,,(j-out_n)] <- if(accept | j == out_n + 1) .$dataf$out[ii, ] else .$dataf$out_mcmc[ii,,(j-out_n-1)]

  }

  # ADAPT CROSSOVER VALUES FUNCTION CALL HERE
  if (.$wpars$mcmc_adapt_CR & (.$mcmc$t < .$mcmc$CR_burnin)) {

    # make sure to fix function names in wrapper object
    # update the probability of the differnt CR values
    .$adapt_pCR()

    .$mcmc$t <- .$mcmc$t + 1

  } else if ((.$wpars$mcmc_adapt_CR) & .$mcmc$t >= .$mcmc$burnin) {

    .$mcmc$t <- .$mcmc$t + 1

  } else if (!(.$wpars$mcmc_adapt_CR) & (j < .$mcmc$CR_burnin)) {
  #if (j < (.$wpars$mcmc_maxiter / 10)) {
  #if ((j < (.$wpars$mcmc_maxiter / 10)) & (sum(.$mcmc$J) > 0)) {

    # PRETTY SURE THIS IS WHERE THE PROBLEM IS!
    # ALJ: make sure this is component-wise division
    # update selection probability of crossover if during crossover burn-in <-- not working
    .$mcmc$p_CR <- .$mcmc$J    / .$mcmc$n_id
    .$mcmc$p_CR <- .$mcmc$p_CR / sum(.$mcmc$p_CR)
  }

  # debug: print crossover values and crossover probabilities at the end of each iteration
  #print(paste0('iteration = ', j))
  #print('id = ')
  #print(.$mcmc$id)
  #print('del = ')
  #print(.$mcmc$del)
  print('CR = ')
  print(.$mcmc$CR)
  print('p_CR = ')
  print(.$mcmc$p_CR)
  #print('n_id =')
  #print(.$mcmc$n_id)
  #print(paste0('m = ', .$mcmc$m))
  #print('L = ')
  #print(.$mcmc$L)

}



# crossover adaption functions
################################

# function that generates/updates crossover values based on current probabilities
generate_CR <- function(.) {

  # sample m from numbers 1,...,n_CR using multinomial distribution (ie, probabilities p_CR)
  .$mcmc$m <- sample(1:.$wpars$mcmc_n_CR, size = 1, replace = T, prob = .$mcmc$p_CR)

  # set crossover probability/value
  .$mcmc$CR <- .$mcmc$m / .$wpars$mcmc_n_CR

  # index of which crosover probabilities/values are selected
  .$mcmc$L[.$mcmc$m] <- .$mcmc$L[.$mcmc$m] + 1

}


# function to compute the squared normalized jumping distance
calc_del <- function(., j, ii) {

  # compute standard deviation of each dimension/parameter
  .$mcmc$sd_state[] <- apply(.$dataf$pars_array, 1, sd)

  # replace any 0's in sd vector with 1e-9 to avoid division by 0
  idx                  <- which(.$mcmc$sd_state == 0)
  .$mcmc$sd_state[idx] <- 1e-9

  # debug: can check summation
  #summation <- sum(((.$dataf$pars_array[ , ii, .$mcmc$t] - .$dataf$pars_array[ , ii, (.$mcmc$t - 1)]) / .$mcmc$sd_state)^2)
  #summation <- sum(((.$dataf$pars_array[1:.$mcmc$d, ii, .$mcmc$t] - .$dataf$pars_array[1:.$mcmc$d, ii, .$mcmc$t - 1]) / .$mcmc$sd_state)^2)
  summation <- sum(((.$dataf$pars_array[1:.$mcmc$d, ii, j] - .$dataf$pars_array[1:.$mcmc$d, ii, j - 1]) / .$mcmc$sd_state)^2)

  .$mcmc$del[.$mcmc$m] <- .$mcmc$del[.$mcmc$m] + summation

  # PROBLEM 1: del vector is zeros <- fixed it!
  # PROBLEM 2: need to index m due to parallelization (ie, need to make it a vector of length chains)

  #print(paste0('.$mcmc$t = ', .$mcmc$t))
  #print(.$dataf$pars_array[ , , .$mcmc$t])
  # problem is because t - 1 = 0 for j = 2
  #print(paste0('j = ', j))
  #print(.$dataf$pars_array[ , , j])
  #print(paste0('j - 1 = ', j - 1))
  #print(.$dataf$pars_array[ , , j - 1])

  #print(paste0('chain = ', ii))
  #print('standard deviation = ')
  #print(.$mcmc$sd_state)
  #print(paste0('summation = ', summation))
  #print(paste0('m = ', .$mcmc$m))
  #print('.$mcmc$del = ')
  #print(.$mcmc$del)

}


# function that adapts crossover probabilities
adapt_pCR <- function(.) {

  # note: dan's code doesn't have the multiplication by t???
  # update the probability of the different CR values being selected
  for (qq in 1:.$wpars$mcmc_n_CR) {
    .$mcmc$p_CR[qq] <- .$mcmc$t * .$wpars$mcmc_chains * (.$mcmc$del[.$mcmc$m] / .$mcmc$L[.$mcmc$m]) / sum(.$mcmc$del)
  }

}



# prior distribution functions
################################

# initialize chains with uniform distributions
mcmc_prior_uniform <- function(.) {

  # number of Markov chains
  n <- .$wpars$mcmc_chains

  # number of parameters (dimensionality of parameter space)
  d <- length(.$dynamic$pars)

  # determine minimums and maximums for parameter values
  min_vals <- rep(0, d); max_vals <- rep(0, d)
  for(i in 1:d) {
    min_vals[i] <- .$dynamic$pars[[names(.$dynamic$pars)[i]]][['min']]
    max_vals[i] <- .$dynamic$pars[[names(.$dynamic$pars)[i]]][['max']]
  }

  # extract parameter names
  pars_names <- names(.$dynamic$pars)

  # draw priors from uniform distribution
  # create pars / proposal matrix
  #.$dataf$pars <- matrix(0, nrow = n, ncol = d)
  for (i in 1:d) {
    #.$dataf$pars[1:n, i] <- runif(n, min = min_vals[i], max = max_vals[i])
    nam <- paste('prior', i, sep = '')
    assign(nam, runif(n, min = min_vals[i], max = max_vals[i]))
  }

  .$dataf$pars <- cbind(prior1, prior2, prior3, prior4)
  print(.$dataf$pars)

  # assign parameter names
  #colnames(.$dataf$pars) <- pars_names
}


# initialize chains with Latin hypercube sampling
mcmc_prior_latin <- function(.) {

}


# initialize chains with multi-normal distribution (multivariate normal distribution)
mcmc_prior_normal <- function(.) {

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
}

# no boundary handling: a good choice when the search space is not "physically" restricted
mcmc_bdry_handling_none <- function(., j, ii, jj) {

  if ((j == .$wpars$mcmc_maxiter) & (ii == .$wpars$mcmc_chains) & (jj == .$mcmc$d)) {
    print('No option was chosen for MCMC boundary handling.')
  }
}


# restrict parameter proposals that are beyond the boundary
# reset them back to the determined bound value
mcmc_bdry_handling_bound <- function(., j, ii, jj) {

  # if outside bound of parameter space, restrict proposal value to corresponding dimension min
  if (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj]) .$dataf$pars[jj, ii] <- .$mcmc$boundary_min[jj]

  # if outside bound of parameter space, restrict proposal to corresponding dimension max
  else if (.$dataf$pars[jj, ii] > .$mcmc$boundary_max[jj]) .$dataf$pars[jj, ii] <- .$mcmc$boundary_max[jj]
}


# restrict parameter proposals that are beyond the boundary
# reflect them back accros the boundary
mcmc_bdry_handling_reflect <- function(., j, ii, jj) {

  if (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj]) {
    .$dataf$pars[jj, ii] <- 2 * .$mcmc$boundary_min[jj] - .$dataf$pars[jj, ii]
  } else if (.$dataf$pars[ii, jj] > .$mcmc$boundary_max[jj]) {
    .$dataf$pars[jj, ii] <- 2 * .$mcmc$boundary_max[jj] - .$dataf$pars[jj, ii]
  }

  # double check: if new reflected proposal value is out of bounds
  if ((.$dataf$pars[jj, ii] > .$mcmc$boundary_max[jj]) | (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj])) {
    .$dataf$pars[jj, ii] <- .$mcmc$boundary_min[jj] + runif(1, min = 0, max = 1) * (.$mcmc$boundary_max[jj] - .$mcmc$boundary_min[jj])
  }
}


# restrict parameter proposals that are beyond the boundary
# "fold" them back across the boundary
# note: this boundary handling approach maintains detailed MCMC balance
mcmc_bdry_handling_fold <- function(., j, ii, jj) {

  if (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj]) {
    .$dataf$pars[jj, ii] <- .$mcmc$boundary_max[jj] - (.$mcmc$boundary_min[jj] - .$dataf$pars[jj, ii])
  } else if (.$dataf$pars[jj, ii] > .$mcmc$boundary_max[jj]) {
    .$dataf$pars[jj, ii] <- .$mcmc$boundary_min[jj] + (.$dataf$pars[jj, ii] - .$mcmc$boundary_max[jj])
  }

  # double check: if new reflected proposal value is out of bounds
  if ((.$dataf$pars[jj, ii] > .$mcmc$boundary_max[jj]) | (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj])) {
    .$dataf$pars[jj, ii] <- .$mcmc$boundary_min[jj] + runif(1, min = 0, max = 1) * (.$mcmc$boundary_max[jj] - .$mcmc$boundary_min[jj])
  }
}



# convergence diagnostic functions
#####################################

# no testing for convergence
mcmc_converge_none <- function(., j) {
  # ALJ: I don't think this print statement will work every time...
  if (j == .$mcmc$burnin) print('No option was chosen to test for MCMC convergence.')
}


# compute the R-statistic of Gelman and Rubin as a convergence diagnostic
mcmc_converge_Gelman_Rubin <- function(., j) {

  #print('Geman-Rubin convergence test being called')

  # need an R_stat storage array! (.$dataf$R_hat)
  # need to add R_hat to output system so I can graph it


  # within-chain variance (for each parameter)

  # between-chain variance (for each parameter)

  #return(R_stat)
}


# auto- correlation function

# Geweke diagnostic (within chain)

# Raftery and Lewis diagnostic (within chain)

# post processing and plotting functions
#########################################


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

  # debugging
  #print('lenth of sspos')
  #print(length(sspos))
  #print('dim of .$dataf$obs')
  #print(dim(.$dataf$obs))
  #print('length of .$dataf$obs')
  #print(length(.$dataf$obs))
  #print('dim of .$dataf$out')
  #print(dim(.$dataf$out))
  #print('dim of transpose .$dataf$out')
  #print(dim(t(.$dataf$out)))
  #print('length of obsse')
  #print(length(obsse))
  #print('type of dataf$out')
  #print(typeof(.$dataf$out))
  #print('type of dataf$obs')
  #print(typeof(.$dataf$obs))
  #print(.$dataf$out)

  # calculate error residual (each chain is on rows of dataf$out, take transpose)
  error_residual_matrix <- ( t(.$dataf$out)[sspos, ] - .$dataf$obs[sspos] ) / obsse
  #error_residual_matrix <- ( .$dataf$out[ , sspos] - .$dataf$obs[sspos] ) / obsse

  #print("dim of error residual matrix")
  #print(dim(error_residual_matrix))

  # calculate sum of squared error
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2))

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  -(obs_n / 2) * log(2 * pi) - sum(log(obsse)) - 0.5 * SSR
}



# outlier handling functions
#####################################

# no outlier handling for Markov chians
mcmc_outlier_none <- function(., j) {
  # NEED TO FIX THIS; PRINT STATEMENT NOT BEING CALLED AT CORRECT TIME
  if (j == .$wpars$mcmc_maxiter) print('No option was chosen to identify and correct outlier Markov chains during burn-in.')
}


# function that detects and corrects outlier Markov chains using the Inter Quartile-Range (IQR) statistic
mcmc_outlier_iqr <- function(., j) {

  # debug/development print statment
  #print('IQR outlier test being called')

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

    # IMPORTANT QUESTOIN: do i replace just the current iteration, or do i replace the whole sampled history???

    # replace outlier chain(s)
    #.$dataf$pars_array[outliers, 1, j] <- .$dataf$pars_array[replace_idx, 1, j]
    .$dataf$pars_array[1:.$mcmc$d, outliers, j] <- .$dataf$pars_array[1:.$mcmc$d, replace_idx, j]
    .$dataf$pars_lklihood[outliers, j] <- .$dataf$pars_lklihood[replace_idx, j]
  }

  # identifying and correcting outliers should only be done during burn-in
  # because it violates the balance of sampled chains and destroys reversibility
  # if outlier chain is detected, apply another burn-in period before generating posterior moments
  # ALJ: not sure how to implement this in code??? basically need to restart j???
}


mcmc_outlier_grubbs <- function(., j) {}


mcmc_outlier_pierce <- function(., j) {}


mcmc_outlier_chauvenet <- function(., j) {}



### END ###
